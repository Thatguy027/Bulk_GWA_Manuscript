# calculate LD with focal sid-2 variant for all of chrIII
plink --vcf bcsq.vcf.gz \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --geno 0.1 \
  --chr III \
  --r2 \
  --ld-snp III:13680248 \
  --ld-window-kb 9999999 \
  --ld-window 9999999 \
  --ld-window-r2 0 \
  --out focal_variant_ld_chrIII_filtered


# Step 1: Filter for genotyping rate ≥90% and create binary files
plink --vcf bcsq.vcf.gz \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --geno 0.1 \
  --make-bed \
  --out filtered_geno10

# Step 2: Calculate inter-chromosomal LD (one-by-all across ALL chromosomes)
plink --bfile filtered_geno10 \
  --r2 inter-chr \
  --allow-extra-chr \
  --ld-snp III:13680248 \
  --ld-window-r2 0 \
  --out focal_variant_ld_genome_wide