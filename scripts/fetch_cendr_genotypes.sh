#!/usr/bin/env bash
#
# Rebuild the CaeNDR-derived genotype inputs from the public releases -------
#
#   scripts/fetch_cendr_genotypes.sh plink     the 20210121 PLINK panel
#   scripts/fetch_cendr_genotypes.sh bcsq      the 20231213 BCSQ-annotated VCF
#   scripts/fetch_cendr_genotypes.sh all       both
#
# WHY THIS SCRIPT EXISTS
# The analysis reads two large genotype resources that are NOT this study's
# data: a PLINK conversion of the CaeNDR 20210121 release (540 isotypes) and a
# bcftools-csq annotation of the CaeNDR 20231213 release. Together they are
# about 25 GB. They are third-party data, so the deposit cites the releases and
# ships this script rather than redistributing them, and a reader rebuilds them
# byte-for-byte from the public files.
#
# Neither conversion was recorded anywhere when it was first run. Both were
# recovered from provenance the files carry themselves: the PLINK command from
# the .log beside each fileset, and the csq command from the VCF header
# (`bcftools view -h ... | grep csqCommand`). Nothing here is reconstructed
# from memory.
#
# WHAT IT COSTS
#   plink  8.3 GB download, ~0.5 GB written (or ~6.5 GB with --with-text)
#   bcsq   8.3 GB download plus a 29 MB reference and an 8 MB GFF3,
#          ~10 GB written; csq takes hours
# The downloads land in $WORK (default $TMPDIR), NOT in the repository, and are
# kept so a re-run resumes rather than re-downloads. Delete $WORK when done.
#
# THE .ped/.map QUESTION
# The original run passed `--recode 01`, which wrote .ped and .map beside every
# fileset: 5.9 GB of the 6.5 GB directory. NOTHING IN THE ANALYSIS READS THEM
# -- verified by grep over every script -- so they are off by default. Pass
# --with-text to reproduce the original directory exactly.
#
# URLs verified 2026-09-22. If a host has moved, the release pages are
# https://caendr.org/data/data-release/c-elegans/20210121 and /20231213.
#
# THE csq STEP IS VERIFIED, not just transcribed. The release VCF is indexed on
# its host, so a 20 kb window around sid-2 was pulled over HTTP, annotated with
# the command below, and compared with the same window of the bcsq.vcf.gz the
# analysis used. The BCSQ strings are identical, including the two-record
# haplotype form at III:13,680,412 that residue 151 turns on:
#
#   missense|sid-2|ZK520.2.1|protein_coding|+|151A>151I|13680412G>A+13680413C>T,
#   missense|sid-2|ZK520.2.1|protein_coding|+|151A>151T|13680412G>A,@13680236
#
# To repeat that check without downloading 8.3 GB:
#   bcftools view -r III:13670000-13690000 -Oz -o w.vcf.gz "$VCF_2023"
#   bcftools csq -Oz --fasta-ref ref.fa --gff-annot annot.gff3 --ncsq 32 \
#       --phase a -o w.bcsq.vcf.gz w.vcf.gz && bcftools index -t w.bcsq.vcf.gz
#   bcftools query -r III:13680412 -f '%POS\t%INFO/BCSQ\n' w.bcsq.vcf.gz
#
# The PLINK step is NOT verified end to end -- that needs the whole 8.3 GB
# release -- so instead it asserts the per-chromosome variant and sample counts
# of the panel the paper used, and stops on any mismatch.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${WORK:-${TMPDIR:-/tmp}/cendr_fetch}"
PLINK_OUT="${PLINK_OUT:-$ROOT/data/genotypes/CeNDR20210121_Plink}"
BCSQ_OUT="${BCSQ_OUT:-$ROOT/data/genotypes/cendr_20231213}"
WITH_TEXT=0

S3="https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_elegans"
VCF_2021="$S3/20210121/variation/WI.20210121.hard-filter.isotype.vcf.gz"
VCF_2023="$S3/20231213/variation/WI.20231213.hard-filter.isotype.vcf.gz"
# WormBase's own host refuses these paths (403); the EBI mirror serves them.
REF_FA="https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS276/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS276.genomic.fa.gz"
GFF3="https://ftp.ensembl.org/pub/release-112/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.112.gff3.gz"

CHROMS=(I II III IV V X)
# Variants each chromosome must end with, from the .bim files of the original
# conversion; the .log beside each one reports the same number. A mismatch
# means the release or the tool changed, and the script says so rather than
# writing a panel that silently differs from the one the paper used.
declare -a WANT=(358837 463722 341971 435133 909818 408516)
N_SAMPLES=540

msg() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

need() {
  command -v "$1" >/dev/null 2>&1 || die "$1 not on PATH. $2"
}

fetch() {  # url dest -- resumable, skipped when already complete
  local url="$1" dest="$2" want
  want="$(curl -sIL --max-time 60 "$url" | tr -d '\r' \
          | awk 'tolower($1)=="content-length:"{l=$2} END{print l+0}')"
  if [[ -f "$dest" ]]; then
    local have; have="$(wc -c < "$dest" | tr -d ' ')"
    if [[ "$want" -gt 0 && "$have" == "$want" ]]; then
      msg "have $(basename "$dest") ($have bytes)"; return
    fi
    msg "resuming $(basename "$dest") ($have of $want bytes)"
  else
    msg "downloading $(basename "$dest") ($want bytes)"
  fi
  curl -fL --retry 3 --retry-delay 5 -C - -o "$dest" "$url"
  if [[ "$want" -gt 0 ]]; then
    local got; got="$(wc -c < "$dest" | tr -d ' ')"
    [[ "$got" == "$want" ]] || die "$(basename "$dest"): got $got bytes, expected $want"
  fi
}

do_plink() {
  need plink "Needs PLINK 1.9 (the original used v1.90b6.21); plink2 will not do -- it has neither --biallelic-only nor --recode 01."
  need bcftools "Needed to split the release VCF per chromosome."
  local ver; ver="$(plink --version 2>/dev/null | head -1)"
  [[ "$ver" == *"v1.9"* ]] || die "found '$ver'; this conversion needs PLINK 1.9."
  msg "using $ver"

  mkdir -p "$WORK" "$PLINK_OUT"
  fetch "$VCF_2021" "$WORK/WI.20210121.hard-filter.isotype.vcf.gz"
  fetch "$VCF_2021.tbi" "$WORK/WI.20210121.hard-filter.isotype.vcf.gz.tbi"

  local i=0 c
  for c in "${CHROMS[@]}"; do
    local part="$WORK/$c.vcf.gz"
    if [[ ! -s "$part" ]]; then
      msg "splitting chromosome $c"
      bcftools view -r "$c" -Oz -o "$part" \
        "$WORK/WI.20210121.hard-filter.isotype.vcf.gz"
    fi
    if [[ ! -s "$PLINK_OUT/$c.bed" ]]; then
      msg "converting $c"
      # exactly the options the original .log records, in its own order
      local args=(--vcf "$part" --allow-extra-chr --biallelic-only --snps-only
                  --set-missing-var-ids '@:#' --output-missing-genotype 9
                  --make-bed --out "$PLINK_OUT/$c")
      [[ "$WITH_TEXT" == 1 ]] && args+=(--recode 01)
      plink "${args[@]}" >/dev/null
    fi
    local got; got="$(wc -l < "$PLINK_OUT/$c.bim" | tr -d ' ')"
    local fam; fam="$(wc -l < "$PLINK_OUT/$c.fam" | tr -d ' ')"
    [[ "$got" == "${WANT[$i]}" ]] \
      || die "$c: $got variants, the paper's panel has ${WANT[$i]}"
    [[ "$fam" == "$N_SAMPLES" ]] \
      || die "$c: $fam samples, expected $N_SAMPLES"
    msg "  $c ok -- $got variants, $fam isotypes"
    i=$((i + 1))
  done
  msg "PLINK panel complete in $PLINK_OUT"
  [[ "$WITH_TEXT" == 1 ]] || msg "  (.ped/.map omitted; nothing in the analysis reads them -- pass --with-text for them)"
}

do_bcsq() {
  need bcftools "The original used bcftools 1.11."
  mkdir -p "$WORK" "$BCSQ_OUT"
  fetch "$VCF_2023" "$WORK/WI.20231213.hard-filter.isotype.vcf.gz"
  fetch "$REF_FA" "$WORK/c_elegans.PRJNA13758.WS276.genomic.fa.gz"
  fetch "$GFF3" "$WORK/Caenorhabditis_elegans.WBcel235.112.gff3.gz"

  # csq wants an uncompressed or bgzipped fasta it can index; the WormBase
  # file is gzip, not bgzip, so it is expanded rather than re-compressed.
  [[ -s "$WORK/ref.fa" ]] || {
    msg "expanding the reference"
    gunzip -c "$WORK/c_elegans.PRJNA13758.WS276.genomic.fa.gz" > "$WORK/ref.fa"
  }
  [[ -s "$WORK/annot.gff3" ]] || gunzip -c "$WORK/Caenorhabditis_elegans.WBcel235.112.gff3.gz" > "$WORK/annot.gff3"

  if [[ ! -s "$BCSQ_OUT/bcsq.vcf.gz" ]]; then
    msg "running bcftools csq -- this takes hours"
    # verbatim from the VCF header of the file the analysis used:
    #   csq -Oz --fasta-ref c_elegans.PRJNA13758.WS276.genomic.fa \
    #       --gff-annot Caenorhabditis_elegans.WBcel235.112.gff3 \
    #       --ncsq 32 --phase a -o bcsq.vcf.gz WI.20231213.hard-filter.isotype.vcf.gz
    bcftools csq -Oz --fasta-ref "$WORK/ref.fa" --gff-annot "$WORK/annot.gff3" \
      --ncsq 32 --phase a -o "$BCSQ_OUT/bcsq.vcf.gz" \
      "$WORK/WI.20231213.hard-filter.isotype.vcf.gz"
    bcftools index -t "$BCSQ_OUT/bcsq.vcf.gz"
  fi
  msg "BCSQ VCF complete: $BCSQ_OUT/bcsq.vcf.gz"
  msg "  scripts that read it take its path from \$CENDR_BCSQ; export it:"
  msg "  export CENDR_BCSQ=$BCSQ_OUT/bcsq.vcf.gz"
}

TARGET=""
for a in "$@"; do
  case "$a" in
    plink|bcsq|all) TARGET="$a" ;;
    --with-text)    WITH_TEXT=1 ;;
    -h|--help)      sed -n '2,45p' "$0"; exit 0 ;;
    *)              die "unknown argument '$a' (expected plink, bcsq, all, --with-text)" ;;
  esac
done
[[ -n "$TARGET" ]] || { sed -n '2,45p' "$0"; exit 1; }

need curl "Needed to fetch the public releases."
msg "work directory: $WORK"
case "$TARGET" in
  plink) do_plink ;;
  bcsq)  do_bcsq ;;
  all)   do_plink; do_bcsq ;;
esac
msg "done. The downloads in $WORK can be deleted."
