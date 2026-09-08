#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Genome-wide significance by phenotype permutation, reconstructing the 2023
 * pos-1 scan exactly.
 *
 *   nextflow run main.nf -profile hoffman2 \
 *     --pheno traits/2023_pos1_association_traits.csv \
 *     --traits vst_ctrl_pos-1_T2 \
 *     --n_perm 1000 --name pos1_perm
 *
 * WHY. The eigen threshold (Li & Ji) and Bonferroni both approximate the
 * multiple-testing burden analytically. Permutation measures it in THIS panel,
 * with its actual linkage disequilibrium and relatedness.
 *
 * HOW. Genotypes are fixed; the phenotype is shuffled among the strains that
 * have one. Each shuffle is mapped genome-wide with the same leave-one-
 * chromosome-out kinship matrices as the real scan, and the genome-wide MAXIMUM
 * -log10 p is recorded. The (1 - alpha) quantile of those maxima is the
 * threshold. The observed phenotype is carried as permutation 0, so its own
 * maximum comes out of the same code path.
 *
 * WHY THE STEPS BELOW LOOK THE WAY THEY DO -- read before changing any flag.
 *
 * A threshold is only meaningful for the scan it is computed against, so this
 * pipeline reproduces that scan rather than doing the same thing in spirit.
 * Three earlier attempts did the latter and each returned a different observed
 * maximum against the scan's 8.8361:
 *
 *   8.6894  MAF computed on all 540 strains, so 519,341 markers were tested
 *   8.5700  kinship built with -gk 1 (centered) where the scan used -gk 2
 *   8.8759  markers correct, panel correct, kinship type correct -- and still
 *           34,316 markers short, because of the last item below
 *
 * The remaining difference was MISSING GENOTYPES, and it is the reason this
 * pipeline now goes through plink's oxford format rather than .traw. The scan's
 * plink step removed 339,834 variants at >10% missingness and kept 464,209 that
 * carry up to 10% each. GEMMA then dropped NONE of them -- but its -miss
 * default is 0.05, which those markers plainly exceed. The explanation is the
 * encoding: `--recode oxford` writes a missing call as "0 0 0", and the dosage
 * expression below turns that into 0, a homozygous call for the second allele.
 * GEMMA therefore never sees a missing value. A .traw route passes "NA" through
 * instead, GEMMA counts it missing, and 34,316 markers fall out.
 *
 * So the oxford route is kept DELIBERATELY, and it is worth being clear that it
 * is not the better choice on its own terms -- it substitutes a fabricated
 * genotype for an absent one, and does so most often on the chromosome arms
 * where the hyper-divergent regions are and where calls fail. It is here
 * because the published scan was built that way and a threshold has to match
 * the scan it thresholds. Anything drawn from a properly-missing scan needs its
 * own threshold, computed by changing this one expression.
 *
 * THE ONE THING THAT MAKES THIS AFFORDABLE. Shuffling the phenotype does not
 * change the genotypes, so the kinship matrices are computed ONCE and reused by
 * every permutation.
 *
 * WHAT PERMUTING A LABEL DOES AND DOES NOT DO. Shuffling phenotype labels
 * destroys the relatedness structure the mixed model is fitted to. The null it
 * samples is "no association AND no population structure", while the model
 * assumes structure is present and corrects for it. This is the standard
 * threshold in the C. elegans GWAS literature and is what cegwas/NemaScan
 * report, but it is not exact: it tends to be slightly ANTI-conservative where
 * structure inflates the real scan. The defensible reading is "a threshold
 * calibrated to this panel's LD and marker density", not "an exact family-wise
 * error rate".
 *
 * COST. Jobs are (n_perm / perm_batch) x 6 chromosomes, plus 6 GRMs and the
 * conversions. At n_perm 1000 and perm_batch 25 that is 240 mapping jobs, each
 * running 25 GEMMA calls in series.
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
 * against params defined in the config file itself. */

/* ------------------------------------------------------------------ */

/* THE PANEL IS EVERY STRAIN IN THE PHENOTYPE FILE, not only the strains with a
 * value for the trait being mapped. That is what the scan did: its plink step
 * kept all 366, so --maf 0.05 and --geno were computed across 366, and a marker
 * carried by two of the 231 phenotyped strains can still be in the scan. GEMMA
 * then drops the strains with no value for the trait, analysing 231.
 *
 * Getting this wrong is the first bug in the list above: filtering on the 231
 * gave 519,341 markers instead of 464,209.
 *
 * The order comes from the scan's own traits.sample rather than from anything
 * recomputed here. BIMBAM has no sample IDs -- dosage columns are positional --
 * so this file is the definition of which column is which strain, and it is
 * checked against plink's own output below rather than trusted. */
process PREP_ORDER {
    publishDir "${params.outdir}/panel", mode: 'copy'
    input:  path phenofile
            path sample_order
    output: path 'keep.txt',         emit: keep
            path 'strain_order.txt', emit: order
            path 'panel_summary.txt'
    script:
    """
    set -euo pipefail
    # oxford .sample carries two header lines before the samples
    sed 1,2d ${sample_order} | awk '{print \$1}' > strain_order.txt
    awk '{print \$1, \$1}' strain_order.txt > keep.txt

    Rscript ${projectDir}/bin/prep_order.R --pheno ${phenofile} \\
        --order strain_order.txt --traits '${params.traits}' \\
        --expect_individuals ${params.expect_individuals}
    """
}

/* plink, with the scan's flags verbatim. Deviations that were tried and are
 * wrong are noted rather than removed, because each one changed the answer:
 *
 *   --snps-only --biallelic-only  restricts to 2,919,294 of the VCF's 3,504,135
 *                                records. Without them, --set-missing-var-ids
 *                                can give an indel at the same position the
 *                                same chr:pos id as the SNP.
 *   --maf 0.05 --geno            computed across all 366 kept strains. --geno
 *                                with no argument is 0.1.
 *   --recode oxford              the missing-genotype encoding described in the
 *                                header. This is the load-bearing flag.
 *   --output-missing-genotype 9  carried over from the scan. It governs text
 *                                formats such as .ped, not the .gen
 *                                probabilities, so it has no effect here -- it
 *                                is kept only so this command is the scan's.
 *
 * --pheno is NOT passed. The scan passed it and plink reported "0 phenotype
 * values present", because the file is strain/trait columns where plink expects
 * FID IID value. It filtered nothing, so omitting it changes no genotype. */
process PLINK_CONVERT {
    publishDir "${params.outdir}/plink", mode: 'copy', pattern: '*.log'
    input:  path vcf
            path keep
    output: tuple path('traits.gen'), path('traits.sample'), emit: gen
            path 'traits.log'
    script:
    """
    set -euo pipefail
    ${params.plink} --vcf ${vcf} \\
        --snps-only --biallelic-only \\
        --maf ${params.maf} --geno \\
        --set-missing-var-ids '@:#' \\
        --keep ${keep} \\
        --output-missing-genotype 9 \\
        --recode oxford \\
        --out traits \\
        --allow-extra-chr \\
        --threads ${task.cpus} --memory 6000

    n_var=\$(wc -l < traits.gen)
    n_ind=\$(( \$(wc -l < traits.sample) - 2 ))
    echo "variants: \$n_var   individuals: \$n_ind"
    if [ "${params.expect_variants}" != "0" ] && [ "\$n_var" != "${params.expect_variants}" ]; then
      echo "ERROR: \$n_var variants, expected ${params.expect_variants}." >&2
      echo "  The scan's plink step retained ${params.expect_variants}; see" >&2
      echo "  snps/pos1_2023_plink_traits.log for its exact counts. A different" >&2
      echo "  number means a different VCF release or a changed filter, and the" >&2
      echo "  threshold would not apply to the scan. --expect_variants 0 skips." >&2
      exit 1
    fi
    if [ "\$n_ind" != "${params.expect_individuals}" ]; then
      echo "ERROR: \$n_ind individuals, expected ${params.expect_individuals}." >&2
      exit 1
    fi
    """
}

/* .gen -> BIMBAM, and the annotation.
 *
 * THE DOSAGE EXPRESSION IS THE SCAN'S, character for character. .gen carries
 * chr, rs, pos, alleleA, alleleB then three genotype probabilities per sample,
 * so sample i occupies columns 3i+3, 3i+4, 3i+5 and 2*P(AA) + P(AB) is its
 * count of allele A. BIMBAM's dosages count the allele listed first after the
 * marker id, which is why alleleA (\$4) precedes alleleB (\$5).
 *
 * A missing call is "0 0 0", so this yields 0 -- see the header. Changing it to
 * emit NA is the one-line change that makes missingness honest, and it produces
 * a DIFFERENT scan needing its own threshold.
 *
 * The annotation deliberately does NOT reproduce the scan's `awk 'NR!=1'`. A
 * .gen file has no header, so that expression silently dropped the scan's first
 * variant: 464,209 retained by plink, 464,208 in the marker list it handed
 * GEMMA. Rather than replicate the bug, the shipped marker list -- which is
 * that 464,208, exactly what GEMMA received -- is passed as -snps, and the
 * annotation is built complete. GEMMA tests the intersection, which is the
 * scan's 464,045 once MtDNA is excluded by -loco. */
process BUILD_BIMBAM {
    publishDir "${params.outdir}/bimbam", mode: 'copy', pattern: '*.tsv'
    input:  tuple path(gen), path(sample)
    output: path 'traits.csv',                  emit: geno
            path 'traits_gemmaAnnotation.tsv',  emit: anno
    script:
    """
    set -euo pipefail
    nsamp=\$(( \$(wc -l < ${sample}) - 2 ))
    awk -v s=\$nsamp '{ printf \$2","\$4","\$5;
                        for (i = 1; i <= s; i++) printf ","\$(i*3+3)*2+\$(i*3+4);
                        printf "\\n" }' ${gen} > traits.csv

    # ANNOTATION CHROMOSOME NAMES MUST MATCH WHAT -loco IS GIVEN.
    #
    # plink's oxford export writes chromosomes by its own numeric codes, and X
    # is a name plink knows: it comes out as 23. The roman numerals I-V are not
    # names plink knows, so --allow-extra-chr passes them through untouched.
    # The result is an annotation reading I, II, III, IV, V, 23, MtDNA.
    #
    # GEMMA's -loco X then matches NOTHING, and it does not fail -- it silently
    # tests every marker in the file against the chromosome-X-excluded kinship.
    # That is how an earlier run reported a genome-wide maximum of 9.0220: the
    # chromosome III peak was re-tested inside the -loco X job, where the
    # kinship still contains chromosome III, and 8.6837 inflated to 9.0220. The
    # six per-chromosome scans were individually correct the whole time.
    #
    # So plink's codes are mapped back to the VCF's names, and then every
    # chromosome the workflow will ask for is asserted to be present. The
    # assertion is the part that matters: a silent fallback that inflates the
    # answer is worth failing the run over.
    cut -f-3 -d' ' ${gen} \\
      | awk 'BEGIN{OFS="\\t"; m["23"]="X"; m["24"]="Y"; m["25"]="XY"; m["26"]="MT"}
             {c=\$1; if (c in m) c=m[c]; print \$2, \$3, c}' \\
      > traits_gemmaAnnotation.tsv

    test "\$(wc -l < traits.csv)" = "\$(wc -l < ${gen})"
    test "\$(wc -l < traits_gemmaAnnotation.tsv)" = "\$(wc -l < ${gen})"

    cut -f3 traits_gemmaAnnotation.tsv | sort -u > .chroms
    echo "chromosomes in the annotation: \$(tr '\\n' ' ' < .chroms)"
    for c in ${params.chromosomes.join(' ')}; do
      if ! grep -qx "\$c" .chroms; then
        echo "ERROR: chromosome '\$c' is not in the annotation." >&2
        echo "  GEMMA's -loco would match nothing and would then test EVERY" >&2
        echo "  marker against that chromosome's kinship instead of failing." >&2
        echo "  Annotation has: \$(tr '\\n' ' ' < .chroms)" >&2
        exit 1
      fi
    done
    """
}

/* One kinship matrix per chromosome, from every marker NOT on it, computed once
 * and reused by every permutation.
 *
 * -gk 2, NOT -gk 1. GEMMA's -gk 1 is the centered relatedness matrix and -gk 2
 * the standardized one, where each marker is divided by its own standard
 * deviation before the cross-product. They are different matrices and give
 * different p-values; the scan used -gk 2. Using -gk 1 put the observed maximum
 * at 8.5700 against 8.8361 with everything else already correct.
 *
 * -loco does the chromosome selection, rather than splitting the genotype file
 * per chromosome as an earlier version did. Same reasoning as everywhere else
 * here: it is what the scan did.
 *
 * GEMMA needs A phenotype column to decide which individuals to include, and
 * uses column 1. The permutation matrix's first batch is passed, whose column 1
 * is the OBSERVED phenotype -- the same 231 strains the scan's column 1 had, so
 * the matrix is identical. The values are irrelevant to -gk; only the
 * individual set is. */
process GEMMA_GRM {
    tag "${chrom}"
    publishDir "${params.outdir}/kinship", mode: 'copy'
    input:  each chrom
            path geno
            path anno
            path snps
            path pheno_placeholder
    output: tuple val(chrom), path("gemmeGRM.${chrom}.sXX.txt"), emit: kin
            path "gemmeGRM.${chrom}.log.txt"
    script:
    """
    set -euo pipefail
    zcat -f ${snps} > snps.txt
    ${params.gemma} -g ${geno} -p ${pheno_placeholder} -gk 2 -loco ${chrom} \\
        -a ${anno} -snps snps.txt -o gemmeGRM.${chrom} -outdir .
    grep -E "Version|analyzed|total SNPs" gemmeGRM.${chrom}.log.txt || true
    """
}

/* Permuted phenotype columns, in strain_order (which is .sample order, which is
 * BIMBAM column order). Column 1 is the OBSERVED phenotype, so the real scan's
 * genome-wide maximum comes from the same code path as the permutations and
 * cannot drift from them. Strains with no value stay NA and GEMMA drops them,
 * exactly as in the scan. */
process MAKE_PERMS {
    tag "${trait}"
    publishDir "${params.outdir}/permutations", mode: 'copy', pattern: '*.tsv'
    input:  path order
            path phenofile
            each trait
    output: tuple val(trait), path("perm_${trait}_b*.txt"), path("strains_${trait}.tsv"), emit: perms
    script:
    """
    awk '{print \$1, \$1}' ${order} > order.fam
    Rscript ${projectDir}/bin/make_permutations.R \\
        --fam order.fam --pheno ${phenofile} --trait '${trait}' \\
        --n_perm ${params.n_perm} --batch ${params.perm_batch} --seed ${params.seed}
    """
}

process GEMMA_PERM {
    tag "${trait}:${chrom}:${batch.baseName}"
    /* The observed scan (permutation 0) is published in full, with its GEMMA
     * log. Everything else keeps only its genome-wide maximum, which is all a
     * threshold needs -- but when the observed maximum does not match the scan
     * being thresholded, a single number gives nothing to diagnose with. The
     * per-marker table can be differenced against the shipped scan by
     * scripts/compare_observed_scan.R, which is how the missing-genotype
     * encoding was found. Perm 0 is one scan per chromosome, so it costs
     * nothing. */
    publishDir "${params.outdir}/observed_scan", mode: 'copy',
               pattern: 'observed_*.assoc.txt.gz'
    publishDir "${params.outdir}/observed_scan", mode: 'copy',
               pattern: 'gemma_*.log.txt'
    input:  tuple val(trait), val(chrom), path(kin), path(batch)
            path geno
            path anno
            path snps
    output: path "maxima_${trait}_${chrom}_${batch.baseName}.tsv", emit: maxima
            path "observed_${trait}_${chrom}.assoc.txt.gz", optional: true
            path "gemma_${trait}_${chrom}.log.txt", optional: true
    script:
    """
    set -euo pipefail
    zcat -f ${snps} > snps.txt
    ncol=\$(awk 'NR==1{print NF}' ${batch})
    offset=\$(echo ${batch.baseName} | sed 's/.*_b//')
    : > maxima_${trait}_${chrom}_${batch.baseName}.tsv
    for k in \$(seq 1 \$ncol); do
      # -lmm 1 is the Wald test ALONE, matching the scan's p_wald column. -lmm 4
      # would emit p_wald, p_lrt and p_score, so a positional \$NF would silently
      # read p_score -- a different statistic from the one being thresholded.
      ${params.gemma} -g ${geno} -p ${batch} -n \$k -lmm 1 -k ${kin} \\
          -loco ${chrom} -a ${anno} -snps snps.txt -o run_\$k -outdir . \\
          > /dev/null 2>&1
      # p_wald is located BY NAME from the header, not by position: GEMMA's
      # column layout differs between -lmm modes and versions.
      # A -loco run must test ONLY its own chromosome. If -loco silently
      # matched nothing, this file would carry the whole genome and its maximum
      # would be a marker scored against the wrong kinship -- the 9.0220 bug.
      # Checked on every call, because it costs one awk pass and the failure it
      # catches produces a plausible number rather than an error.
      n_other=\$(awk -F'\t' -v want="${chrom}" \\
        'NR==1{for(i=1;i<=NF;i++) if(\$i=="chr") c=i; next} c && \$c!=want {n++}
         END{print n+0}' run_\$k.assoc.txt)
      if [ "\$n_other" -ne 0 ]; then
        echo "ERROR: the -loco ${chrom} scan tested \$n_other markers that are" >&2
        echo "  not on ${chrom}. GEMMA's -loco matched no chromosome and fell" >&2
        echo "  back to the whole genome; check the annotation's chromosome" >&2
        echo "  names against ${chrom}." >&2
        exit 1
      fi
      mx=\$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="p_wald") c=i; next}
             c && \$c!="" {p=\$c+0; if(p>0 && (m==""||p<m)) m=p}
             END{if(m=="") print "NA"; else printf "%.6f", -log(m)/log(10)}' \\
             run_\$k.assoc.txt)
      pid=\$(( offset + k - 1 ))
      printf "%s\\t%s\\t%s\\t%s\\n" "${trait}" "${chrom}" "\$pid" "\$mx" \\
          >> maxima_${trait}_${chrom}_${batch.baseName}.tsv
      if [ "\$pid" -eq 0 ]; then
        gzip -c run_\$k.assoc.txt > observed_${trait}_${chrom}.assoc.txt.gz
        cp run_\$k.log.txt gemma_${trait}_${chrom}.log.txt
      fi
      rm -f run_\$k.assoc.txt run_\$k.log.txt
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
    snps      = file(params.snps,  checkIfExists: true)
    order     = file(params.sample_order, checkIfExists: true)
    traits    = params.traits.tokenize(',')*.trim()

    PREP_ORDER(phenofile, order)
    PLINK_CONVERT(vcf, PREP_ORDER.out.keep)
    BUILD_BIMBAM(PLINK_CONVERT.out.gen)

    MAKE_PERMS(PREP_ORDER.out.order, phenofile, traits)

    // the first batch of the first trait doubles as the -gk placeholder
    grm_pheno = MAKE_PERMS.out.perms
        .map { trait, batches, strains -> (batches instanceof List ? batches[0] : batches) }
        .first()

    GEMMA_GRM(params.chromosomes, BUILD_BIMBAM.out.geno,
              BUILD_BIMBAM.out.anno, snps, grm_pheno)

    jobs = MAKE_PERMS.out.perms
        .flatMap { trait, batches, strains ->
            (batches instanceof List ? batches : [batches]).collect { b -> tuple(trait, b) } }
        .combine(GEMMA_GRM.out.kin)
        .map { trait, batch, chrom, kin -> tuple(trait, chrom, kin, batch) }

    GEMMA_PERM(jobs, BUILD_BIMBAM.out.geno, BUILD_BIMBAM.out.anno, snps)
    COLLECT_THRESHOLD(GEMMA_PERM.out.maxima.collect())
}
