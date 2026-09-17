#!/usr/bin/env python3
"""Assemble the figures the manuscript actually cites into manuscript_figures/.

The draft refers to figures by manuscript number (Figure 1, Figure S1, ...)
while the repo names them by content (Figure1_pos1, SUPP_FIG_XX_simulation_depth,
...). This script is the one place that mapping lives. It

  1. resolves each manuscript number to its source .pdf/.png in plots/,
  2. pulls that figure's caption out of FIGURE_CAPTIONS.txt by block header
     rather than copying the text, so captions stay live,
  3. writes manuscript_figures/ with numbered copies, FIGURES.md (shareable,
     PNGs inline) and MANIFEST.tsv,
  4. names every figure in plots/ that the draft does NOT cite, so a figure
     cannot drop out of the package silently.

Re-runnable: manuscript_figures/ is rebuilt from scratch each time.
Nothing outside manuscript_figures/ is written.

Usage:  python3 scripts/finalize_manuscript_figures.py [--check]
        --check verifies sources and captions resolve, writes nothing.
"""
from __future__ import annotations
import argparse
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PLOTS = ROOT / "plots"
CAPTIONS = ROOT / "FIGURE_CAPTIONS.txt"      # working record, with caveats
CLEAN = ROOT / "MANUSCRIPT_CAPTIONS.txt"     # publication text, preferred
OUT = ROOT / "manuscript_figures"

# label, source basename in plots/, string that uniquely identifies the caption
# block header, and where the draft cites it.
MANIFEST = [
    ("Figure 1",  "Figure1_pos1",                                  "Figure1_pos1.pdf",
     "Fig 1A MIP-seq vs NNLS; 1B pos-1 response across 231 isolates; 1C pilot GWAS"),
    ("Figure 2",  "Figure2",                                       "plots/Figure2.pdf",
     "nine-target pooled GWAS with cross-QTL arrowheads"),
    ("Figure 3",  "Figure3_quad",                                  "Figure3_quad.pdf",
     "Fig 3A cross parents in the pooled assay; 3B chrIII allele frequency in the F2 pool; 3C NIL series"),
    ("Figure 4",  "Figure4_sid2",                                  "Figure4_sid2.pdf",
     "Fig 4A JU reciprocal edits; 4B N2 edits at 25% pos-1; 4C ectodomain local net charge"),
    ("Figure S1", "SUPP_FIG_XX_simulation_depth",                  "SUPP_FIG_XX_simulation_depth",
     "NNLS recovery against sequencing depth, seven traits"),
    ("Figure S2", "SUPP_FIG_XX_dilution_validation",               "SUPP_FIG_XX_dilution_validation",
     "inferred against designed pool ratios, r = 0.997"),
    ("Figure S3", "SUPP_FIG_XX_downsample_per_sample",             "SUPP_FIG_XX_downsample_per_sample",
     "agreement with MIP-seq under read downsampling"),
    ("Figure S4", "SUPP_FIG_XX_original_pos1_dfreq_rep_correlation","SUPP_FIG_XX_original_pos1_dfreq_rep_correlation",
     "replicate agreement of the 2023 pooled pos-1 response"),
    ("Figure S5", "SUPP_FIG_plate_vs_paaby_vs_pos1original",       "SUPP_FIG_plate_vs_paaby_vs_pos1original",
     "S5A manual plate scores against pooled; S5B manual against Paaby et al. 2015"),
    ("Figure S6", "SUPP_FIG_XX_pooled_phenotype_ranks",            "SUPP_FIG_XX_pooled_phenotype_ranks",
     "ranked pooled mig-6 and pos-1 responses with the cross parents marked"),
    # S7 is the cross contrasts because the draft has cited it as S7 since before
    # the NIL figure existed, and because Figure 2's section cites it well ahead
    # of the NIL work in the SID-2 section -- this list is in citation order.
    ("Figure S7", "SUPP_FIG_XX_cross_contrast_panels",             "SUPP_FIG_XX_cross_contrast_panels",
     "S7B-C the two cross contrast panels"),
    ("Figure S8", "SUPP_FIG_XX_nil_hatching_full",                 "SUPP_FIG_XX_nil_hatching_full",
     "the full NIL hatching experiment, all ten strains on both food conditions"),
    ("Figure S9", "SUPP_FIG_XX_n2_swap_dose",                      "SUPP_FIG_XX_n2_swap_dose",
     "N2 residue-96 swap across the pos-1 dose series"),
    # One figure, two panels: (A) the cross-species alignment of the sequon
    # window, (B) the N94A / 96K editing series. Cited as S10A and S10B.
    ("Figure S10","SUPP_FIG_XX_sid2_ortholog_conservation",        "SUPP_FIG_XX_sid2_ortholog_conservation",
     "S10A the N94-C95-T96 sequon across Caenorhabditis; S10B the N94A editing series on pos-1 RNAi"),
    ("Figure S11","SUPP_FIG_XX_sid2_local_charge",                 "SUPP_FIG_XX_sid2_local_charge",
     "local net charge percentile across the ectodomain"),
]

DASHES = re.compile(r"^-{40,}\s*$")

# Captions are copied VERBATIM -- this script never rewrites the author's prose.
# But FIGURE_CAPTIONS.txt is a working file, so some captions orient themselves
# against neighbouring blocks or against repo paths, which is meaningless once
# the caption is packaged on its own. These are reported, not fixed.
INTERNAL_REF = [
    (re.compile(r"\bthe block above\b", re.I),      "refers to a neighbouring caption block"),
    (re.compile(r"\bAs Figure \d"),                  "defers to another figure's caption"),
    (re.compile(r"\bSUPP_FIG[A-Za-z0-9_]*"),         "names a repo figure slug, not a manuscript number"),
    (re.compile(r"\b(plots|scripts)/[A-Za-z0-9_.]+"),"names a repo path"),
    (re.compile(r"\bsuperseded\b", re.I),           "mentions a superseded version"),
]


def lint_caption(label: str, cap: str) -> list[str]:
    out = []
    for ln in cap.splitlines():
        for pat, why in INTERNAL_REF:
            m = pat.search(ln)
            if m:
                out.append(f"{label}: {why} -- {m.group(0)!r} in: {ln.strip()[:88]}")
                break
    return out


def caption_blocks(path: Path) -> list[tuple[str, str]]:
    """Split FIGURE_CAPTIONS.txt into (header, body) pairs on its rule lines."""
    lines = path.read_text().splitlines()
    rules = [i for i, ln in enumerate(lines) if DASHES.match(ln)]
    chunks = []
    for a, b in zip(rules, rules[1:]):
        chunk = "\n".join(lines[a + 1:b]).strip("\n")
        chunks.append(chunk)
    out = []
    for i, chunk in enumerate(chunks):
        first = chunk.strip().splitlines()[0] if chunk.strip() else ""
        if first.startswith(("FIGURE ", "SUPP_FIG")) and i + 1 < len(chunks):
            out.append((chunk.strip(), chunks[i + 1].strip("\n")))
    return out


def clean_captions(path: Path) -> dict[str, str]:
    """Parse MANUSCRIPT_CAPTIONS.txt on its '## <label>' headings."""
    if not path.exists():
        return {}
    out, label, buf = {}, None, []
    for ln in path.read_text().splitlines():
        if ln.startswith("## "):
            if label:
                out[label] = "\n".join(buf).strip("\n")
            label, buf = ln[3:].strip(), []
        elif label is not None:
            buf.append(ln)
        # lines before the first heading are the file's own header comment
    if label:
        out[label] = "\n".join(buf).strip("\n")
    return {k: v.rstrip() for k, v in out.items()}


def find_caption(blocks, key: str, label: str) -> str:
    hits = [(h, b) for h, b in blocks if key in h]
    if len(hits) != 1:
        sys.exit(f"ERROR {label}: caption key {key!r} matched {len(hits)} blocks "
                 f"(need exactly 1). Fix MANIFEST or FIGURE_CAPTIONS.txt.")
    return hits[0][1].rstrip()


def sort_key(label: str):
    m = re.match(r"Figure (S?)(\d+)([A-Z]?)$", label)
    return (1 if m.group(1) else 0, int(m.group(2)), m.group(3))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="verify sources and captions resolve, write nothing")
    args = ap.parse_args()

    if not CAPTIONS.exists():
        sys.exit(f"ERROR: {CAPTIONS} not found")
    blocks = caption_blocks(CAPTIONS)
    clean = clean_captions(CLEAN)

    resolved, problems, fellback = [], [], []
    for label, base, key, cited in MANIFEST:
        png, pdf = PLOTS / f"{base}.png", PLOTS / f"{base}.pdf"
        for f in (png, pdf):
            if not f.exists():
                problems.append(f"{label}: missing {f.relative_to(ROOT)}")
        if label in clean:
            cap = clean[label]
        else:
            cap = find_caption(blocks, key, label)
            fellback.append(label)
        resolved.append((label, base, png, pdf, cap, cited))
    extra = sorted(set(clean) - {m[0] for m in MANIFEST})
    if extra:
        problems.append(f"{CLEAN.name} has captions with no MANIFEST entry: "
                        + ", ".join(extra))
    if problems:
        sys.exit("ERROR:\n  " + "\n  ".join(problems))

    used = {base for _, base, _, _, _, _ in resolved}
    everything = sorted(p.stem for p in PLOTS.glob("*.pdf"))
    unused = [b for b in everything if b not in used]

    src = f"{CLEAN.name} ({len(resolved) - len(fellback)})"
    if fellback:
        src += f" + {CAPTIONS.name} ({len(fellback)}: {', '.join(fellback)})"
    print(f"resolved {len(resolved)} figures, all sources present")
    print(f"captions from {src}")
    if unused:
        print(f"\nin plots/ but NOT cited by the draft ({len(unused)}) -- not packaged:")
        for b in unused:
            print(f"  {b}")
    if args.check:
        print("\n--check: nothing written")
        return 0

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()

    rows, md = [], []
    md.append("# Manuscript figures\n")
    md.append("Figures cited by the current draft, in citation order.\n")
    md.append("PNGs are shown inline; a print-resolution PDF sits beside each one.\n")
    md.append("Generated by `scripts/finalize_manuscript_figures.py` from "
              "`MANUSCRIPT_CAPTIONS.txt` -- edit the captions or the figure "
              "scripts and re-run rather than editing this file.\n")
    md.append("\n---\n")

    for label, base, png, pdf, cap, cited in sorted(resolved, key=lambda r: sort_key(r[0])):
        slug = label.replace("Figure ", "Figure_").replace(" ", "_")
        shutil.copy2(png, OUT / f"{slug}.png")
        shutil.copy2(pdf, OUT / f"{slug}.pdf")
        rows.append((label, base, f"{slug}.png", f"{slug}.pdf", cited))
        md.append(f"\n## {label}\n")
        md.append(f"![{label}]({slug}.png)\n")
        md.append(f"*[Print-resolution PDF]({slug}.pdf)*\n")
        md.append(f"\n{cap}\n")
        md.append("\n---\n")

    if unused:
        md.append("\n## Not included\n")
        md.append("\nThese figures exist in `plots/` but the current draft does not "
                  "cite them, so they are not packaged here:\n\n")
        for b in unused:
            md.append(f"- `{b}`\n")
        md.append("\nAdd one by appending it to `MANIFEST` in "
                  "`scripts/finalize_manuscript_figures.py` and re-running.\n")

    (OUT / "FIGURES.md").write_text("".join(md))
    with (OUT / "MANIFEST.tsv").open("w") as fh:
        fh.write("label\tsource\tpng\tpdf\tcited_as\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    n = len(list(OUT.glob("*.png")))
    print(f"\nwrote {OUT.relative_to(ROOT)}/  "
          f"({n} figures, FIGURES.md, MANIFEST.tsv)")

    lint = []
    for label, _, _, _, cap, _ in resolved:
        lint += lint_caption(label, cap)
    if lint:
        print(f"\nCAPTION LINT -- {len(lint)} line(s) carry repo-internal references.")
        print("Captions ship verbatim, so these reach whoever you send FIGURES.md to.")
        print(f"Fix them in {CLEAN.name} and re-run; this script will not edit prose.")
        for msg in lint:
            print(f"  {msg}")
    else:
        print("\ncaption lint: clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
