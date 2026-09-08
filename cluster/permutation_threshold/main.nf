#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Genome-wide significance by phenotype permutation.
 *
 *   nextflow run main.nf -profile hoffman2 \
 *     --pheno traits/2023_pos1_association_traits.csv \
 *     --traits vst_ctrl_pos-1_T2 \
 *     --n_perm 1000 --name pos1_perm
 *
 * WHY. The eigen threshold (Li & Ji) and Bonferroni both approximate the
 * multiple-testing burden analytically. Permutation measures it in THIS panel,
 * with its actual linkage disequilibrium and relatedness, and needs no
 * assumption about how many independent tests there are.
 *
 * HOW. Genotypes are fixed; the phenotype is shuffled across strains. Each
 * shuffle is mapped genome-wide with the same leave-one-chromosome-out kinship
 * matrices as the real scan, and the genome-wide MAXIMUM -log10 p is recorded.
 * The (1 - alpha) quantile of those maxima is the threshold. The observed
 * phenotype is carried as permutation 0, so its own maximum and empirical
 * p-value come out of the same machinery.
 *
 * THE ONE THING THAT MAKES THIS AFFORDABLE. Shuffling the phenotype does not
 * change the genotypes, so the kinship matrices are computed ONCE and reused by
 * every permutation. Recomputing them per permutation would multiply the cost
 * by the number of permutations for no change in the result.
 *
 * WHAT PERMUTING A LABEL DOES AND DOES NOT DO -- read before quoting the number.
 * Shuffling phenotype labels destroys the relatedness structure the mixed model
 * is fitted to. The null it samples is therefore "no association AND no
 * population structure", while the model assumes structure is present and
 * corrects for it. The resulting threshold is the standard one in the C. elegans
 * GWAS literature and is what cegwas/NemaScan report, but it is not exact: it
 * tends to be slightly ANTI-conservative where structure inflates the real
 * scan, because the permuted scans have no structure to inflate them. The
 * defensible reading is "a threshold calibrated to this panel's LD and marker
 * density", not "an exact family-wise error rate". A structure-preserving
 * alternative is to permute in the space rotated by the eigenvectors of the
 * kinship matrix; that is not implemented here and would be the thing to do if
 * a reviewer presses.
 *
 * COST. Jobs are (n_perm / perm_batch) x 6 chromosomes, plus 6 GRMs and the
 * conversions. At n_perm 1000 and perm_batch 25 that is 240 mapping jobs, each
 * running 25 GEMMA calls in series. Batching exists because 6000 one-call jobs
 * would spend most of the wall clock in the SGE queue rather than in GEMMA.
 */

def required(name, value) {
    if (value == null) {
        log.error "--${name} is required"
        System.exit(1)
    }
    return value
}

/* Param defaults are declared in nextflow.config, NOT here. A ${params.x}
 * interpolated inside a process or profile block in the config can only resolve
 * against params defined in the config file itself. Declaring them in both
 * places invites the two copies to drift.
 */

/* ------------------------------------------------------------------ */

/* The panel: strains carrying a value for the requested trait(s). MAF must be
 * computed among THESE strains, not all 540 -- a marker at 5% in the full
 * collection can be below 5% in the phenotyped subset and vice versa. The first
 * version filtered on all 540 and produced 519,341 markers against the scan's
 * 464,045, which is why its observed maximum came out at 8.69 instead of 8.84. */
process PREP_PANEL {
    publishDir "${params.outdir}/panel", mode: 'copy'
    input:  path phenofile
    output: path 'panel.txt', emit: keep
            path 'panel_summary.txt'
    script:
    """
    Rscript ${projectDir}/bin/prep_panel.R --pheno ${phenofile} \\
        --traits '${params.traits}'
    """
}

/* MARKER SET. Rather than re-derive the scan's filter chain -- which is not
 * recoverable from the archived output and which the first attempt got wrong --
 * the permutation scan tests EXACTLY the markers the real scan tested, supplied
 * as an id list. The multiple-testing burden then matches by construction
 * instead of by luck, and the threshold provably applies to that scan.
 *
 * No --maf or --geno here for the same reason: the extract list already defines
 * the set, and any further filter would silently shrink it. The assertion below
 * fails the run immediately if the count does not match, rather than after a
 * thousand permutations. */
process PLINK_CONVERT {
    publishDir "${params.outdir}/plink", mode: 'copy', pattern: '*.log'
    input:  path vcf
            path keep
            path markers
    output: tuple path('all.bed'), path('all.bim'), path('all.fam'), emit: bed
            path 'all.log'
    script:
    """
    set -euo pipefail
    zcat -f ${markers} > markers.txt
    ${params.plink} --vcf ${vcf} --allow-extra-chr --set-missing-var-ids '@:#' \\
        --keep ${keep} --extract markers.txt \\
        --make-bed --out all --threads ${task.cpus} --memory 6000

    n_mk=\$(wc -l < all.bim)
    n_id=\$(wc -l < all.fam)
    echo "markers retained: \$n_mk   strains retained: \$n_id"
    if [ "${params.expect_markers}" != "0" ] && [ "\$n_mk" != "${params.expect_markers}" ]; then
      echo "ERROR: \$n_mk markers, expected ${params.expect_markers}." >&2
      echo "  The permutation scan must test the same markers as the scan it" >&2
      echo "  thresholds. Check --markers, or set --expect_markers 0 to skip." >&2
      exit 1
    fi
    """
}

/* Per-chromosome BIMBAM for mapping, and the complement for the LOCO kinship.
 * BIMBAM is built from plink's A-transpose (.traw) rather than --recode bimbam,
 * because .traw's column layout is unambiguous: CHR SNP CM POS COUNTED ALT then
 * one dosage column per sample, in .fam order. */
process BUILD_CHROM {
    tag "${chrom}"
    input:  tuple path(bed), path(bim), path(fam)
            each chrom
    output: tuple val(chrom), path("geno_${chrom}.bimbam"), path("anno_${chrom}.txt"),
                  path("notchr_${chrom}.bimbam"), emit: sets
    script:
    """
    set -euo pipefail
    ${params.plink} --bfile all --allow-extra-chr --chr ${chrom} \\
        --recode A-transpose --out chr_${chrom} --threads ${task.cpus} --memory 6000
    ${params.plink} --bfile all --allow-extra-chr --not-chr ${chrom} \\
        --recode A-transpose --out not_${chrom} --threads ${task.cpus} --memory 6000

    # .traw -> BIMBAM geno: "snp, minor, major, dosage..."  (dosage is ALT count)
    # .traw is: CHR SNP (C)M POS COUNTED ALT <one dosage column per sample>,
    # and each dosage is the count of the COUNTED allele (\$5). BIMBAM's dosages
    # count the allele listed FIRST after the marker id, so \$5 must come before
    # \$6. Reversing them flips the allele coding, which changes the sign of beta;
    # it leaves p_wald untouched, so it would not have broken this threshold --
    # but it would quietly corrupt any effect size taken from these files.
    awk 'NR>1 {printf "%s, %s, %s", \$2, \$5, \$6; for(i=7;i<=NF;i++) printf ", %s", \$i; printf "\\n"}' \\
        chr_${chrom}.traw > geno_${chrom}.bimbam
    awk 'NR>1 {printf "%s, %s, %s", \$2, \$5, \$6; for(i=7;i<=NF;i++) printf ", %s", \$i; printf "\\n"}' \\
        not_${chrom}.traw > notchr_${chrom}.bimbam

    # annotation: snp, position, chromosome
    awk 'NR>1 {print \$2 ", " \$4 ", " \$1}' chr_${chrom}.traw > anno_${chrom}.txt

    test -s geno_${chrom}.bimbam
    test -s notchr_${chrom}.bimbam
    """
}

/* One kinship matrix per chromosome, from every marker NOT on it. Computed once
 * and reused by every permutation -- see the header.
 *
 * -gk 2, NOT -gk 1. GEMMA's -gk 1 is the centered relatedness matrix and -gk 2
 * the standardized one, where each marker is divided by its own standard
 * deviation before the cross-product. They are different matrices, they give
 * different p-values, and the scan being thresholded here used -gk 2 (the lab
 * gemma_nf pipeline's GEMMA_GRM process, and its archived gemmeGRM.*.sXX.txt
 * output). Using -gk 1 put the observed maximum at 8.5700 against the scan's
 * 8.8361. Nothing about the marker set or the panel was wrong at that point --
 * both matched exactly -- so this is the whole of the remaining discrepancy.
 * The output filename follows the flag: -gk 1 writes .cXX.txt, -gk 2 .sXX.txt. */
process GEMMA_GRM {
    tag "${chrom}"
    publishDir "${params.outdir}/kinship", mode: 'copy'
    input:  tuple val(chrom), path(geno), path(anno), path(notchr)
            path pheno_placeholder
    output: tuple val(chrom), path("kin_${chrom}.sXX.txt"), emit: kin
            path "kin_${chrom}.log.txt"
    script:
    """
    set -euo pipefail

    # The split must be exact: every marker belongs to this chromosome or to the
    # kinship set, never to both and never to neither. A silent drop here would
    # change the kinship without changing anything the marker-count assertion in
    # PLINK_CONVERT can see.
    n_chr=\$(wc -l < ${geno})
    n_not=\$(wc -l < ${notchr})
    echo "chr ${chrom}: \$n_chr markers tested, \$n_not in the kinship"
    if [ "${params.expect_markers}" != "0" ] \\
       && [ \$(( n_chr + n_not )) != "${params.expect_markers}" ]; then
      echo "ERROR: \$n_chr + \$n_not != ${params.expect_markers}" >&2
      exit 1
    fi

    ${params.gemma} -g ${notchr} -p ${pheno_placeholder} -gk 2 -o kin_${chrom}
    mv output/kin_${chrom}.sXX.txt .
    # The log carries GEMMA's version banner and the counts it actually analysed,
    # which is the evidence for whether the kinship was built from the markers
    # intended. Published, because a kinship is not self-describing.
    mv output/kin_${chrom}.log.txt .
    grep -E "Version|analyzed|total SNPs" kin_${chrom}.log.txt || true
    """
}

/* Permuted phenotype columns, in .fam order. Column 1 is the OBSERVED
 * phenotype, so the real scan's genome-wide maximum comes from the same code
 * path as the permutations and cannot drift from them. */
process MAKE_PERMS {
    tag "${trait}"
    publishDir "${params.outdir}/permutations", mode: 'copy', pattern: '*.tsv'
    input:  tuple path(bed), path(bim), path(fam)
            path phenofile
            each trait
    output: tuple val(trait), path("perm_${trait}_*.txt"), path("strains_${trait}.tsv"), emit: perms
    script:
    """
    Rscript ${projectDir}/bin/make_permutations.R \\
        --fam all.fam --pheno ${phenofile} --trait '${trait}' \\
        --n_perm ${params.n_perm} --batch ${params.perm_batch} --seed ${params.seed}
    """
}

/* One job per (trait, chromosome, batch). Runs GEMMA once per column in the
 * batch, keeping only each column's maximum -- the per-marker output of a
 * permutation is never needed and would be terabytes. */
process GEMMA_PERM {
    tag "${trait}:${chrom}:${batch.baseName}"
    /* The observed scan (permutation 0) is published in full. Everything else
     * keeps only its genome-wide maximum, which is all a threshold needs -- but
     * when the observed maximum does not match the scan being thresholded, a
     * single number gives nothing to diagnose with. The per-marker table can be
     * differenced against the shipped scan directly, which turns "it is 0.04
     * off" into a statement about which markers and by how much. Perm 0 is one
     * scan per chromosome, so this costs nothing. */
    publishDir "${params.outdir}/observed_scan", mode: 'copy',
               pattern: 'observed_*.assoc.txt.gz'
    publishDir "${params.outdir}/observed_scan", mode: 'copy',
               pattern: 'gemma_*.log.txt'
    input:  tuple val(trait), val(chrom), path(geno), path(anno), path(kin), path(batch)
    output: path "maxima_${trait}_${chrom}_${batch.baseName}.tsv", emit: maxima
            path "observed_${trait}_${chrom}.assoc.txt.gz", optional: true
            path "gemma_${trait}_${chrom}.log.txt", optional: true
    script:
    """
    set -euo pipefail
    ncol=\$(awk 'NR==1{print NF}' ${batch})
    offset=\$(echo ${batch.baseName} | sed 's/.*_b//')
    : > maxima_${trait}_${chrom}_${batch.baseName}.tsv
    for k in \$(seq 1 \$ncol); do
      # -lmm 1 is the Wald test ALONE, matching the shipped scan's p_wald column.
      # -lmm 4 would emit p_wald, p_lrt and p_score, so a positional \$NF would
      # silently read p_score -- a different statistic from the one the
      # threshold is meant to apply to.
      ${params.gemma} -g ${geno} -p ${batch} -a ${anno} -k ${kin} \\
          -lmm 1 -n \$k -o run_\$k > /dev/null 2>&1
      # p_wald is located BY NAME from the header, not by position: GEMMA's
      # column layout differs between -lmm modes and versions.
      mx=\$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="p_wald") c=i; next}
             c && \$c!="" {p=\$c+0; if(p>0 && (m==""||p<m)) m=p}
             END{if(m=="") print "NA"; else printf "%.6f", -log(m)/log(10)}' \\
             output/run_\$k.assoc.txt)
      pid=\$(( offset + k - 1 ))
      printf "%s\\t%s\\t%s\\t%s\\n" "${trait}" "${chrom}" "\$pid" "\$mx" \\
          >> maxima_${trait}_${chrom}_${batch.baseName}.tsv
      # permutation 0 IS the observed phenotype -- keep its scan and its GEMMA
      # log, the latter because the version banner and the analyzed marker and
      # individual counts are what a mismatch is diagnosed from
      if [ "\$pid" -eq 0 ]; then
        gzip -c output/run_\$k.assoc.txt > observed_${trait}_${chrom}.assoc.txt.gz
        cp output/run_\$k.log.txt gemma_${trait}_${chrom}.log.txt
      fi
      rm -f output/run_\$k.assoc.txt output/run_\$k.log.txt
    done
    """
}

process COLLECT_THRESHOLD {
    publishDir "${params.outdir}", mode: 'copy'
    input:  path maxima
    output: path 'permutation_maxima.tsv'
            path 'permutation_thresholds.tsv'
            path 'permutation_threshold.pdf'
            path 'permutation_threshold.png'
    script:
    """
    cat ${maxima} > all_maxima.tsv
    Rscript ${projectDir}/bin/collect_threshold.R \\
        --maxima all_maxima.tsv --alpha '${params.alpha.join(",")}' \\
        --n_perm ${params.n_perm} \\
        --expect_observed_max ${params.expect_observed_max} \\
        --observed_tol ${params.observed_tol}
    """
}

workflow {
    required('pheno', params.pheno)
    phenofile = file(params.pheno, checkIfExists: true)
    vcf       = file(params.vcf,   checkIfExists: true)
    traits    = params.traits.tokenize(',')*.trim()

    markers = file(params.markers, checkIfExists: true)
    PREP_PANEL(phenofile)
    PLINK_CONVERT(vcf, PREP_PANEL.out.keep, markers)
    BUILD_CHROM(PLINK_CONVERT.out.bed, params.chromosomes)

    // GEMMA needs SOME phenotype column to compute a kinship matrix; the values
    // are irrelevant to -gk, so a column of the observed trait is passed. The
    // matrix depends on genotypes only.
    MAKE_PERMS(PLINK_CONVERT.out.bed, phenofile, traits)

    // the first batch of the first trait doubles as the -gk placeholder
    grm_pheno = MAKE_PERMS.out.perms
        .map { trait, batches, strains -> (batches instanceof List ? batches[0] : batches) }
        .first()
    GEMMA_GRM(BUILD_CHROM.out.sets, grm_pheno)

    // (trait, chrom, geno, anno, kin, batch)
    per_chrom = BUILD_CHROM.out.sets
        .map { chrom, geno, anno, notchr -> tuple(chrom, geno, anno) }
        .join(GEMMA_GRM.out.kin)

    jobs = MAKE_PERMS.out.perms
        .flatMap { trait, batches, strains ->
            (batches instanceof List ? batches : [batches]).collect { b -> tuple(trait, b) } }
        .combine(per_chrom)
        .map { trait, batch, chrom, geno, anno, kin ->
            tuple(trait, chrom, geno, anno, kin, batch) }

    GEMMA_PERM(jobs)
    COLLECT_THRESHOLD(GEMMA_PERM.out.maxima.collect())
}
