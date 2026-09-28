#!/usr/bin/env python3
"""sid-2 variants carried by any set of wild isolates.

    python3 scripts/sid2_strain_variants.py XZ1516 ECA191 MY2741
    python3 scripts/sid2_strain_variants.py JU1793 JU2466 --tsv out.tsv
    python3 scripts/sid2_strain_variants.py CB4856 --all

Reads only tracked files, so it runs from a clone:

    supplemental_data/genotypes/sid2_region.{bed,bim,fam}   540 isotypes,
        81 variants over III:13,679,014-13,681,764, from CaeNDR 20210121
    supplemental_data/structure/sid2_population_missense.tsv  BCSQ missense
    supplemental_data/structure/sid2_variants_cendr.tsv       curated labels

ORIENTATION IS AGAINST N2, NOT AGAINST PLINK
`plink2 --export A` names each column with the allele it happened to count,
which is not always the reference allele, so a dosage of 2 does not mean "alt".
N2 is the reference strain, so its own call at each site IS the reference
allele, and every genotype here is reported as agreeing with N2 or differing
from it. Pass --reference to compare against some other strain instead.

MISSING CALLS ARE REPORTED, NEVER SILENTLY COUNTED AS REFERENCE. A strain with
no data at a site is shown as `./.`, and the summary counts it separately. The
distinction matters: "same as N2" and "we do not know" look identical if you
only test for difference.

RESIDUE 151 IS RESOLVED PER STRAIN, ALWAYS AGAINST N2
bcftools csq is haplotype-aware, so III:13,680,412 carries two records: with
the partner SNP at 13,680,413 the codon reads 151I, without it 151T, and
without either it is the reference 151A. The population label keeps both
alternates ("A151I/T"); this script reports the residue each strain actually
carries. That resolution is done against N2 whatever --reference is set to,
because the residue a strain carries is a property of its own haplotype, not
of whoever it is being compared with. N2 is therefore always loaded, even when
it is not one of the strains asked for.

ANNOTATION COVERAGE, which is the real limit here
Genotypes are the 20210121 release; consequences come from the 20231213 BCSQ
annotation. Positions align, but a variant with no consequence record is shown
with a blank protein change rather than as non-coding -- absence of an
annotation is not evidence of a silent variant. The curated table covers only
variants segregating among the four cross parents, so isolates outside that
set routinely carry changes it does not name.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEM = os.path.join(ROOT, "supplemental_data", "genotypes", "sid2_region")
MISSENSE = os.path.join(ROOT, "supplemental_data", "structure",
                        "sid2_population_missense.tsv")
CURATED = os.path.join(ROOT, "supplemental_data", "structure",
                       "sid2_variants_cendr.tsv")

# What the tracked fileset is expected to hold. Pinned so that a rebuilt or
# swapped panel is caught here rather than producing a quietly different answer.
N_VARIANTS, N_SAMPLES = 81, 540
SPAN = (13679014, 13681764)

# Sites where the residue depends on a second SNP in the same codon.
# pos -> (partner position, residue when alone, residue when both)
HAPLO = {13680412: (13680413, "T", "I")}


def die(msg):
    sys.exit(f"ERROR: {msg}")


def load_fileset():
    for ext in ("bed", "bim", "fam"):
        if not os.path.exists(f"{STEM}.{ext}"):
            die(f"missing {STEM}.{ext} -- run scripts/make_supplemental_data.R")
    bim = {}
    for ln in open(f"{STEM}.bim"):
        _c, vid, _cm, pos, a1, a2 = ln.split()
        bim[vid] = (int(pos), a1, a2)
    fam = [ln.split()[0] for ln in open(f"{STEM}.fam")]
    if len(bim) != N_VARIANTS or len(fam) != N_SAMPLES:
        die(f"{os.path.basename(STEM)} holds {len(bim)} variants over {len(fam)} "
            f"samples; expected {N_VARIANTS} over {N_SAMPLES}. The panel changed, "
            "so the pinned numbers in this script and in METHODS need revisiting.")
    lo = min(p for p, _, _ in bim.values())
    hi = max(p for p, _, _ in bim.values())
    if (lo, hi) != SPAN:
        die(f"span is {lo}-{hi}, expected {SPAN[0]}-{SPAN[1]}")
    return bim, fam


def annotations():
    miss, cur = {}, {}
    if os.path.exists(MISSENSE):
        for r in csv.DictReader(open(MISSENSE), delimiter="\t"):
            miss[int(r["pos"])] = r
    if os.path.exists(CURATED):
        for r in csv.DictReader(open(CURATED), delimiter="\t"):
            cur[int(r["pos"])] = r
    return miss, cur


def export_dosages(strains, work):
    if not shutil.which("plink2"):
        die("plink2 not on PATH")
    keep = os.path.join(work, "keep.txt")
    with open(keep, "w") as fh:
        for s in strains:
            fh.write(f"{s} {s}\n")
    out = os.path.join(work, "gt")
    r = subprocess.run(["plink2", "--bfile", STEM, "--keep", keep, "--export", "A",
                        "--out", out, "--allow-extra-chr"],
                       capture_output=True, text=True)
    if not os.path.exists(out + ".raw"):
        die("plink2 produced no .raw\n" + r.stdout[-800:] + r.stderr[-800:])
    lines = open(out + ".raw").read().splitlines()
    hdr = lines[0].split("\t")
    cols = [(i, h) for i, h in enumerate(hdr) if h.startswith("III:")]
    dose = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        f = ln.split("\t")
        dose[f[1]] = {h: f[i] for i, h in cols}
    return dose, cols


def genotype(dose, strain, col, counted, other):
    """Homozygous call as a single letter, 'het', or None when not called."""
    d = dose.get(strain, {}).get(col, "NA")
    if d in ("NA", ""):
        return None
    d = int(float(d))
    return counted if d == 2 else other if d == 0 else "het"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("strains", nargs="+", help="isotype names, e.g. XZ1516 ECA191")
    ap.add_argument("--reference", default="N2",
                    help="strain to compare against (default N2, the reference strain)")
    ap.add_argument("--all", action="store_true",
                    help="show every site, not only those that differ")
    ap.add_argument("--tsv", help="also write the table here")
    args = ap.parse_args()

    bim, fam = load_fileset()
    miss_ann, curated = annotations()

    known = {s.upper(): s for s in fam}
    wanted = []
    for s in args.strains:
        if s in fam:
            wanted.append(s)
        elif s.upper() in known:
            wanted.append(known[s.upper()])
        else:
            near = [f for f in fam if f.upper().startswith(s.upper()[:3])][:8]
            die(f"{s!r} is not in the panel of {len(fam)} isotypes"
                + (f". Did you mean: {', '.join(near)}?" if near else ""))
    if args.reference not in fam:
        die(f"reference {args.reference!r} is not in the panel")
    # N2 is always loaded: it anchors the residue resolution below regardless
    # of what --reference is, and it is the reference strain of the assembly.
    order = list(dict.fromkeys([args.reference, "N2"] + wanted))
    if "N2" not in fam:
        die("N2 is not in the panel; residue resolution needs it")

    with tempfile.TemporaryDirectory() as work:
        dose, cols = export_dosages(order, work)

    show = [s for s in wanted if s != args.reference] or [args.reference]

    # Residue at a codon carrying two SNPs, per strain, anchored on N2.
    # by_pos maps a position to its column and allele pair, so both the focal
    # site and its partner can be read for any strain.
    by_pos = {}
    for i, h in cols:
        vid, counted = h.rsplit("_", 1)
        p, a1, a2 = bim[vid]
        by_pos[p] = (h, counted, a2 if counted == a1 else a1)

    def residue_at(pos, strain):
        """Which residue `strain` carries at a HAPLO site: ref, alone, or both."""
        partner, alone, both = HAPLO[pos]
        if pos not in by_pos or partner not in by_pos:
            return None
        h, counted, other = by_pos[pos]
        focal_s, focal_n2 = (genotype(dose, x, h, counted, other) for x in (strain, "N2"))
        ph, pcounted, pother = by_pos[partner]
        part_s, part_n2 = (genotype(dose, x, ph, pcounted, pother) for x in (strain, "N2"))
        if focal_s is None:
            return "?"
        if focal_s == focal_n2:
            # no change at the focal site: the reference residue, whose letter
            # is the one the label starts with, e.g. the A of A151I/T
            return "ref"
        if part_s is None:
            return "?"
        return both if part_s != part_n2 else alone

    rows = []
    for i, h in cols:
        vid, counted = h.rsplit("_", 1)
        pos, a1, a2 = bim[vid]
        other = a2 if counted == a1 else a1
        ref = genotype(dose, args.reference, h, counted, other)
        calls = {s: genotype(dose, s, h, counted, other) for s in show}
        pc = (curated[pos]["label"] if pos in curated
              else miss_ann[pos]["aa_change"] if pos in miss_ann else "")
        # per-strain residue where the codon holds a second SNP
        resolved = {}
        if pos in HAPLO and pc:
            ref_letter = pc[0] if pc else "?"
            for s in show:
                r = residue_at(pos, s)
                if r is not None:
                    resolved[s] = ref_letter if r == "ref" else r
        af = (curated[pos]["af"] if pos in curated
              else miss_ann[pos]["af"] if pos in miss_ann else "")
        topo = miss_ann[pos]["topology"] if pos in miss_ann else ""
        rows.append(dict(pos=pos, ref=ref, calls=calls, pc=pc, af=af,
                         topo=topo, resolved=resolved))
    rows.sort(key=lambda r: r["pos"])

    sel = [r for r in rows
           if args.all or any(v is not None and v != r["ref"] for v in r["calls"].values())]

    w = f"{'position':>10}  {'ref':3} {'protein change':15} {'AF':>6}  " + \
        "  ".join(f"{s:>9}" for s in show)
    print(f"sid-2, III:{SPAN[0]:,}-{SPAN[1]:,} -- {len(rows)} variants in the panel, "
          f"reference = {args.reference}")
    print(w)
    print("-" * len(w))
    per = {s: {"diff": [], "coding": [], "nocall": 0} for s in show}
    for r in sel:
        cells = []
        for s in show:
            v = r["calls"][s]
            if v is None:
                cells.append("./.")
                per[s]["nocall"] += 1
            elif v == r["ref"]:
                cells.append("=")
            else:
                cells.append(r["resolved"].get(s) and f"{v}" or v)
                per[s]["diff"].append(r["pos"])
                if r["pc"]:
                    lbl = r["pc"]
                    if s in r["resolved"]:
                        lbl = f"{lbl} ({r['resolved'][s]})"
                    per[s]["coding"].append(lbl)
        af = f"{float(r['af']):.3f}" if r["af"] else ""
        print(f"{r['pos']:>10}  {r['ref'] or '?':3} {r['pc']:15} {af:>6}  "
              + "  ".join(f"{c:>9}" for c in cells))

    print(f"\n{len(sel)} of {len(rows)} variants "
          + ("shown" if args.all else f"differ from {args.reference}"))
    for s in show:
        d = per[s]
        line = (f"  {s:<9} {len(d['diff']):>3} variants, "
                f"{len(d['coding'])} protein-altering")
        if d["coding"]:
            line += ": " + ", ".join(d["coding"])
        if d["nocall"]:
            line += f"  [{d['nocall']} sites not called]"
        print(line)

    if args.tsv:
        with open(args.tsv, "w", newline="") as fh:
            out = csv.writer(fh, delimiter="\t")
            out.writerow(["pos", f"ref_{args.reference}", "protein_change",
                          "residue_resolved", "topology", "af_cendr"] + show)
            for r in sel:
                out.writerow([r["pos"], r["ref"] or "", r["pc"],
                              ";".join(f"{k}={v}" for k, v in r["resolved"].items()),
                              r["topo"], r["af"]]
                             + [r["calls"][s] or "./." for s in show])
        print(f"\nwrote {args.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
