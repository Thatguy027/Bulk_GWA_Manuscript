#!/usr/bin/env python3
"""Build the tracked-changes view of the Results section.

    python3 claude_science/results_review/make_tracked_changes.py

Reads results_draft_v1.txt and applies the edit list below. Every `old` string
must occur exactly once in the draft -- the script asserts this, so an edit
can never silently apply to the wrong sentence or fail to apply at all.

Writes, all inside claude_science/results_review/:
    results_tracked_changes.md   full text, ~~deleted~~ -> **[inserted]**
    results_corrected.txt        clean text with every edit applied
    edits_table.tsv              find/replace pairs for the Google Doc

Evidence for each edit is in claims_verified.tsv, keyed by the same claim id.
"""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
DRAFT = HERE / "results_draft_v1.txt"

# (claim_id, severity, old, new, why)
# severity: ERROR   the number or word contradicts the deposit
#           JUDGE   defensible but the repository's own record words it differently
#           STYLE   typo or wording, no numeric consequence
EDITS = [
    ("S3", "ERROR",
     "This simulation framework revealed that NNLS regression can accurately infer strain frequencies with as little as 1X sequencing depth (Figure S1).",
     "This simulation framework revealed that NNLS regression recovers strain frequencies with a mean r\u00b2 of 0.76 at 1X sequencing depth (range 0.51-0.91 across the seven traits), rising above 0.95 from 10X (Figure S1).",
     "simulation_seeded_r2.tsv: mean r2 per trait at 1x is 0.51-0.91, median 0.76; the mean first reaches 0.95 at 10x. SUPP_FIG_XX_simulation_depth.R prints this under the heading 'accuracy at 1x, the depth the text claims'."),

    ("D1", "ERROR",
     "two pooled populations composed of ~48 strains each",
     "two pooled populations composed of 46 and 40 strains",
     "dilution_strain_sets.tsv: 174 isolates in four sets, A=46, B=46, C=40, D=42. The titration is set B against set C, so the two populations are 46 and 40, not ~48 each."),

    ("P2a", "ERROR",
     "We grew these populations for two generations in these conditions across two replicates,",
     "We grew these populations for two generations in these conditions across four pos-1 replicate pools and two control pools,",
     "pos1_2023_sample_frequencies.csv.gz: pos-1 T2 has replicates 1,2,3,4 and ctrl T2 has A,B. Four pos-1 pools is also what makes the six pairwise comparisons cited in the next sentence possible - two replicates give only one pair."),

    ("C1", "ERROR",
     "against each of nine target genes",
     "against each of ten target genes",
     "pooled_vst_traits.csv.gz and bundle $pheno both carry ten RNAi targets at 93 strains each: fog-2, mig-6, pos-1, rde-3, ric-3, rpn-12, spe-19, spe-43, unc-39, vha-5."),

    ("N4", "ERROR",
     "raised its hatching rate from 4.5% to 18.4%",
     "raised its hatching rate from 5.4% to 18.4%",
     "ju_allele_swaps_hatching.csv, pos-1 condition: JU2466_A[96K] is 11 hatched of 204 = 5.4%. Figure4_sid2.R prints 0.054. The 18.4% for wSZ206 is correct. Looks like a digit transposition."),

    ("P4b", "JUDGE",
     "the center of chromosome X that passed the bonferroni significance threshold",
     "the left arm of chromosome X (X:4,875,969) that passed the Bonferroni significance threshold",
     "The single chromosome X marker clearing Bonferroni is at 4,875,969, which is 28% along a 17.72 Mb chromosome. Also fixes the lower-case 'bonferroni'."),

    ("P4d", "JUDGE",
     "and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C).",
     "and a third QTL on the right arm of chromosome III above the eigen-decomposition threshold (Figure 1C). A tenth Bonferroni-passing marker, on chromosome III at 5.97 Mb, has no supporting marker within 100 kb of it and we do not carry it forward.",
     "FIGURE_REPORT.md 'GWAS interval admission' settles this: the 5.966 Mb marker clears Bonferroni at 8.68 with zero eigen-passing neighbours in 100 kb, while the 12.70-12.80 Mb cluster peaks below Bonferroni at 6.31 with 14. The draft's framing follows the repository's admission rule, but a reader comparing against Figure 1C will see a red marker on chromosome III that the text does not explain."),

    ("N6", "JUDGE",
     "we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4% (Figure 4B)",
     "we found that the 96K allele lowered the N2 hatching rate from 32.3% to 4.0% (two independently edited lines, 4.4% and 3.7%; Figure 4B)",
     "n2_allele_swaps_hatching.tsv at the 25% dose: N2[96T] 32.3% (n=220), and the pooled 96K value 4.0% is two independent lines wSZ203 4.4% (n=273) and wSZ204 3.7% (n=295); 273+295=568, the pooled n. The pooled figure is right; naming the two lines matches the caption and is the stronger claim."),

    ("M1", "JUDGE",
     "(Spearman\u2019s \u03c1 = 0.97, n = 98 strains, p < 1e-4)",
     "(Spearman\u2019s \u03c1 = 0.974, n = 98 strains)",
     "Recomputed rho is 0.974 over 98 strains. Figure 1A deliberately reports a bootstrap interval rather than a p value - Figure1_common.R notes that a p value against rho = 0 is not the question. Quoting p < 1e-4 in the text reintroduces what the figure dropped."),

    ("T1", "STYLE",
     "it is unclear that these QTL are pos-1-specifc QTL",
     "it is unclear whether these QTL are pos-1-specific QTL",
     "typo: specifc -> specific; 'unclear that' -> 'unclear whether'."),

    ("T2", "STYLE",
     "we identified RNA-specific QTL",
     "we identified RNAi-specific QTL",
     "RNA-specific should be RNAi-specific; the contrast is between RNAi treatments."),

    ("T3", "STYLE",
     "lowered the strains hatching rate",
     "lowered the strain\u2019s hatching rate",
     "possessive apostrophe."),

    ("T4", "STYLE",
     "increasing the net charge of it\u2019s local environment",
     "increasing the net charge of its local environment",
     "it\u2019s -> its."),

    ("T5", "STYLE",
     "a single-pass transmembrane protein localized apical membrane of intestinal cells",
     "a single-pass transmembrane protein localized to the apical membrane of intestinal cells",
     "missing preposition."),

    ("T6", "STYLE",
     "localize the QTL to a 37 kb interval spanning 13.658 - 13.695.",
     "localize the QTL to a 37 kb interval spanning 13.658-13.695 Mb.",
     "missing units. The wSZ191 introgression in nil_introgression_ranges.bed is 13,657,700-13,695,000 = 37.3 kb, so both the width and the bounds are right."),

    ("T7", "STYLE",
     "Additionally, we identified  QTL that are shared between conditions",
     "Additionally, we identified QTL that are shared between conditions",
     "double space."),

    ("T8", "STYLE",
     "we calculated the per-reside local net charge",
     "we calculated the per-residue local net charge",
     "typo: reside -> residue."),
]


def main():
    if not DRAFT.exists():
        sys.exit(f"draft not found: {DRAFT}")
    text = DRAFT.read_text()

    problems = []
    for cid, sev, old, new, why in EDITS:
        n = text.count(old)
        if n != 1:
            problems.append(f"  {cid}: matches {n} times, expected 1 -> {old[:70]!r}")
    if problems:
        sys.exit("edit strings did not match uniquely:\n" + "\n".join(problems))

    corrected = text
    tracked = text
    for cid, sev, old, new, why in EDITS:
        corrected = corrected.replace(old, new, 1)
        tracked = tracked.replace(old, f"~~{old}~~ **[{new}]**", 1)

    (HERE / "results_corrected.txt").write_text(corrected)

    md = ["# Results section, tracked changes", "",
          "Struck text is the current draft; **[bracketed]** text is the proposed",
          "replacement. Every change is keyed to a claim id in `claims_verified.tsv`,",
          "which carries the recomputed value and the file it came from.", "",
          f"{sum(1 for e in EDITS if e[1] == 'ERROR')} corrections of fact, "
          f"{sum(1 for e in EDITS if e[1] == 'JUDGE')} wording calls for the author, "
          f"{sum(1 for e in EDITS if e[1] == 'STYLE')} typos.", "",
          "---", "", tracked, "", "---", "", "## Why each change", ""]
    for cid, sev, old, new, why in EDITS:
        md += [f"**{cid} - {sev}**", "", f"- current: {old}", f"- proposed: {new}",
               f"- reason: {why}", ""]
    (HERE / "results_tracked_changes.md").write_text("\n".join(md))

    with open(HERE / "edits_table.tsv", "w") as fh:
        fh.write("id\tseverity\tfind\treplace\treason\n")
        for cid, sev, old, new, why in EDITS:
            f = old.replace("\n", " ").replace("\t", " ")
            r = new.replace("\n", " ").replace("\t", " ")
            fh.write(f"{cid}\t{sev}\t{f}\t{r}\t{why}\n")

    print(f"{len(EDITS)} edits applied, all matched uniquely")
    for sev in ("ERROR", "JUDGE", "STYLE"):
        ids = [e[0] for e in EDITS if e[1] == sev]
        print(f"  {sev:6} {len(ids):2}  {', '.join(ids)}")
    print("\nwrote results_tracked_changes.md, results_corrected.txt, edits_table.tsv")


if __name__ == "__main__":
    main()
