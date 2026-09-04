#!/usr/bin/env python3
"""Recover the reported r-squared values from the original simulation figures.

    python3 scripts/extract_sim_reported_r2.py
      -> supplemental_data/deconvolution/simulation_reported_r2.tsv

WHY THIS EXISTS. The simulation compared NNLS-estimated strain frequencies
against the known input frequencies, and only the comparison's *summary*
survives: the expected input values themselves are not in the archive, and
neither is the simulation script. What does survive is one PDF per trait in
data/experiments/initial_sims/plots, each an eight-panel observed-against-
expected figure with "r^2= x" drawn into every panel.

This reads those numbers back out of the PDF content streams, so panel A of
SUPP_FIG_XX_simulation_depth.R is rebuildable from the archive rather than
transcribed by hand. It is the only route to those values.

PANEL ORDER. The generating script drew depths in the order 500, 100, 50, 30,
10, 5, 3, 1 -- the row order of nnls_estFreqs.RDS -- and PDF text operators
appear in drawing order, so the eight r-squared strings map onto DEPTHS below
in sequence. Verified against PC1.pdf read visually: 1, 1, 1, 1, 0.99, 0.98,
0.98, 0.91 left to right.

R cannot do this conveniently: the streams contain embedded NULs, which
rawToChar() refuses.
"""
import csv
import glob
import os
import re
import sys
import zlib

PLOTS = "data/experiments/initial_sims/plots"
OUT = "supplemental_data/deconvolution/simulation_reported_r2.tsv"
DEPTHS = [500, 100, 50, 30, 10, 5, 3, 1]
TJ = re.compile(rb"\((?:[^()\\]|\\.)*\)\s*Tj")


def shown_text(pdf_path):
    """Every string the PDF draws, in drawing order."""
    raw = open(pdf_path, "rb").read()
    out = []
    for m in re.finditer(rb"stream\r?\n(.*?)endstream", raw, re.S):
        blob = m.group(1)
        for wbits in (15, -15):          # zlib header, then raw deflate
            try:
                data = zlib.decompress(blob, wbits)
                break
            except zlib.error:
                data = None
        if data is None:
            continue
        for tok in TJ.findall(data):
            out.append(tok[1:tok.rfind(b")")].decode("latin-1", "replace"))
    return out


def main():
    pdfs = sorted(glob.glob(os.path.join(PLOTS, "*.pdf")))
    if not pdfs:
        sys.exit(f"no PDFs under {PLOTS} -- this needs the Dryad archive")
    rows = []
    for f in pdfs:
        vals = [float(s.split("=")[1])
                for s in (t.strip() for t in shown_text(f)) if "r^2" in s]
        if len(vals) != len(DEPTHS):
            sys.exit(f"{f}: found {len(vals)} r-squared strings, expected "
                     f"{len(DEPTHS)} -- panel layout has changed, check by eye")
        if vals[0] != 1.0:
            sys.exit(f"{f}: first panel r-squared is {vals[0]}, not 1.0 -- the "
                     "panel order assumption (500x first) may be wrong")
        trait = os.path.basename(f)[:-4]
        rows += [(trait, d, v) for d, v in zip(DEPTHS, vals)]

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "depth", "r2"])
        w.writerows(rows)
    traits = sorted({r[0] for r in rows})
    print(f"wrote {OUT}: {len(rows)} rows, {len(traits)} traits")
    for t in traits:
        vs = [f"{v:g}" for _, d, v in rows if _ == t]
        print(f"  {t:<24} {' '.join(vs)}")


if __name__ == "__main__":
    main()
