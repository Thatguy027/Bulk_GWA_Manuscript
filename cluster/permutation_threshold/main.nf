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

params.vcf         = '/u/project/kruglyak/thatguy0/genomics/vcf/WI.20210121.hard-filter.isotype_with_cM.vcf.gz'
params.pheno       = null
params.traits      = 'vst_ctrl_pos-1_T2'   // comma-separated; one threshold per trait
params.maf         = 0.05
params.chromosomes = ['I', 'II', 'III', 'IV', 'V', 'X']
params.n_perm      = 1000
params.perm_batch  = 25
params.seed        = 1
params.alpha       = [0.05, 0.10, 0.01]
params.name        = null
params.outdir      = "results_perm_${new Date().format('yyyyMMdd')}${params.name ? '_' + params.name : ''}"
params.plink       = '/u/project/kruglyak/thatguy0/bin/plink'
params.gemma       = '/u/project/kruglyak/thatguy0/bin/gemma'
params.r_env_bin   = '/u/project/kruglyak/thatguy0/conda/envs/gemma_plots/bin'
params.conda_bin   = '/u/project/kruglyak/thatguy0/conda/bin'

/* ------------------------------------------------------------------ */

process PLINK_CONVERT {
    publishDir "${params.outdir}/plink", mode: 'copy', pattern: '*.log'
    input:  path vcf
    output: tuple path('all.bed'), path('all.bim'), path('all.fam'), emit: bed
            path 'all.log'
    script:
    """
    ${params.plink} --vcf ${vcf} --allow-extra-chr --set-missing-var-ids '@:#' \\
        --snps-only --biallelic-only strict --maf ${params.maf} --geno 0.10 \\
        --make-bed --out all --threads ${task.cpus}
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
        --recode A-transpose --out chr_${chrom} --threads ${task.cpus}
    ${params.plink} --bfile all --allow-extra-chr --not-chr ${chrom} \\
        --recode A-transpose --out not_${chrom} --threads ${task.cpus}

    # .traw -> BIMBAM geno: "snp, minor, major, dosage..."  (dosage is ALT count)
    awk 'NR>1 {printf "%s, %s, %s", \$2, \$6, \$5; for(i=7;i<=NF;i++) printf ", %s", \$i; printf "\\n"}' \\
        chr_${chrom}.traw > geno_${chrom}.bimbam
    awk 'NR>1 {printf "%s, %s, %s", \$2, \$6, \$5; for(i=7;i<=NF;i++) printf ", %s", \$i; printf "\\n"}' \\
        not_${chrom}.traw > notchr_${chrom}.bimbam

    # annotation: snp, position, chromosome
    awk 'NR>1 {print \$2 ", " \$4 ", " \$1}' chr_${chrom}.traw > anno_${chrom}.txt

    test -s geno_${chrom}.bimbam
    test -s notchr_${chrom}.bimbam
    """
}

/* One kinship matrix per chromosome, from every marker NOT on it. Computed once
 * and reused by every permutation -- see the header. */
process GEMMA_GRM {
    tag "${chrom}"
    publishDir "${params.outdir}/kinship", mode: 'copy'
    input:  tuple val(chrom), path(geno), path(anno), path(notchr)
            path pheno_placeholder
    output: tuple val(chrom), path("kin_${chrom}.cXX.txt"), emit: kin
    script:
    """
    set -euo pipefail
    ${params.gemma} -g ${notchr} -p ${pheno_placeholder} -gk 1 -o kin_${chrom}
    mv output/kin_${chrom}.cXX.txt .
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
    input:  tuple val(trait), val(chrom), path(geno), path(anno), path(kin), path(batch)
    output: path "maxima_${trait}_${chrom}_${batch.baseName}.tsv", emit: maxima
    script:
    """
    set -euo pipefail
    ncol=\$(awk 'NR==1{print NF}' ${batch})
    offset=\$(echo ${batch.baseName} | sed 's/.*_b//')
    : > maxima_${trait}_${chrom}_${batch.baseName}.tsv
    for k in \$(seq 1 \$ncol); do
      ${params.gemma} -g ${geno} -p ${batch} -a ${anno} -k ${kin} \\
          -lmm 4 -n \$k -o run_\$k > /dev/null 2>&1
      # smallest p_wald in the chromosome -> largest -log10 p
      mx=\$(awk 'NR>1 && \$NF!="" {p=\$NF; if(p>0 && (m==""||p<m)) m=p} END{
              if(m=="") print "NA"; else printf "%.6f", -log(m)/log(10)}' \\
              output/run_\$k.assoc.txt)
      pid=\$(( offset + k - 1 ))
      printf "%s\\t%s\\t%s\\t%s\\n" "${trait}" "${chrom}" "\$pid" "\$mx" \\
          >> maxima_${trait}_${chrom}_${batch.baseName}.tsv
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
        --n_perm ${params.n_perm}
    """
}

workflow {
    required('pheno', params.pheno)
    phenofile = file(params.pheno, checkIfExists: true)
    vcf       = file(params.vcf,   checkIfExists: true)
    traits    = params.traits.tokenize(',')*.trim()

    PLINK_CONVERT(vcf)
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
