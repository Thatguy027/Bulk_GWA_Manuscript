#!/usr/bin/env python3
"""Assemble the data deposit, and archive what nothing reads.

    python3 scripts/make_deposit.py build     -> deposit/
    python3 scripts/make_deposit.py archive   -> moves unread files to archive_unused/
    python3 scripts/make_deposit.py verify    -> runs a figure inside deposit/

BUILD USES HARD LINKS FOR THE BULK, which is the whole design. data/ is ~12 GB
and this volume has ~15 GB free, so a copy is impossible; a hard link gives the
same bytes two names and costs nothing. Both paths stay valid, editing through
one is visible through the other, and deleting one does not touch the data.
Small text -- scripts, supplemental_data, the prose files -- is COPIED instead,
so nobody edits the deposit and silently changes the repository.

WHAT GOES IN
  data/               the files scripts actually read, per DATA_DEPENDENCY_AUDIT.tsv
  external/           the ~503 MB read from outside the repository
  scripts/            the analysis code, copied
  supplemental_data/  the staged reproduction set, copied
  *.txt, *.md         methods, captions, data availability

WHAT STAYS OUT, and how a reader gets it
The two CaeNDR releases are third-party and are cited, not redistributed;
scripts/fetch_cendr_genotypes.sh rebuilds both derived inputs from the public
files. The 207-strain expression set and the WormBase WS283 annotation are
likewise cited. See DATA_AVAILABILITY.md.

TARRING THE RESULT MATERIALISES THE BYTES. Hard links cost nothing on disk but
a tarball writes every file, so `tar czf deposit.tar.gz deposit/` needs ~9 GB
free on top of what is already used. Upload from the tree, or tar to an
external volume.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEPOSIT = os.path.join(ROOT, "deposit")
ARCHIVE = os.path.join(ROOT, "archive_unused")
AUDIT = os.path.join(ROOT, "DATA_DEPENDENCY_AUDIT.tsv")

# Read from outside the repository. Each entry is (source, path inside
# deposit/external, the environment variable that points at it, why).
EXTERNAL = [
    ("/Users/Stefan/UCLA/Projects/bulkGWAS/xqtl_analysis/NJX_rnai/plots/"
     "JU2466_XZ1516_F2-2_contrast_pos1-par1_10000_plot_DF.tsv",
     "jx_cross/JU2466_XZ1516_F2-2_contrast_pos1-par1_10000_plot_DF.tsv",
     "JX_PLOTS", "JU2466 x XZ1516 cross, pos-1 against par-1"),
    ("/Users/Stefan/UCLA/Projects/bulkGWAS/xqtl_analysis/NJX_rnai/plots/"
     "JU2466_XZ1516_F2-2_contrast_pos1-mig6_10000_plot_DF.tsv",
     "jx_cross/JU2466_XZ1516_F2-2_contrast_pos1-mig6_10000_plot_DF.tsv",
     "JX_PLOTS", "JU2466 x XZ1516 cross, pos-1 against mig-6"),
    ("/Users/Stefan/UCLA/Projects/bulkGWAS/xqtl_analysis/NJX_rnai/plots/"
     "JU2466_XZ1516_F2-2_contrast_mig6-par1_10000_plot_DF.tsv",
     "jx_cross/JU2466_XZ1516_F2-2_contrast_mig6-par1_10000_plot_DF.tsv",
     "JX_PLOTS", "JU2466 x XZ1516 cross, mig-6 against par-1"),
    ("/Users/Stefan/UCLA/Projects/bulkGWAS/baugh_wgs/cluster_data/"
     "20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData",
     "baugh_wgs/20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData",
     "BAUGH_BOOT", "the 103-strain L1 bootstrap input the platform comparison starts from"),
    ("/Users/Stefan/UCLA/Projects/bulkGWAS/traits_with_validated_qtl",
     "traits_with_validated_qtl",
     "SIM_TRAITS", "the trait tables the simulation draws its fitness vectors from"),
    ("/Users/Stefan/github_repos/xQTLSims/data/geneticMapXQTLsnplist.rds",
     "xqtl_sims/geneticMapXQTLsnplist.rds",
     "XQTL_GMAP", "the genetic map the cross null simulation expands"),
]

# Copied rather than linked: small, and a reader editing the deposit must not
# reach back into the repository.
COPY_TREES = ["scripts", "supplemental_data", "githooks"]
COPY_FILES = ["METHODS.txt", "METHODS_CONDENSED.txt", "FIGURE_CAPTIONS.txt",
              "MANUSCRIPT_CAPTIONS.txt", "DATA_AVAILABILITY.md", "README.md",
              "FIGURE_REPORT.Rmd", "CLAUDE.md"]

# Variables whose default already points at a CaeNDR file the deposit does not
# carry; env.sh leaves them for the reader to set after running the fetch script.
FETCHED = {
    "CENDR_BCSQ": "bcsq.vcf.gz, from scripts/fetch_cendr_genotypes.sh bcsq",
    "CENDR_DIR": "the directory holding bcsq.vcf.gz and the WBcel235 GFF3",
    "CENDR_PLINK": "filtered_geno10, from the CaeNDR 20231213 release",
    "CENDR_DIVERGENT": "20231213_c_elegans_divergent_regions_strain.bed",
    "CE207_DIR": "Ce207expression.csv and ce207_qtl.tsv (cited, not deposited)",
    "WS_GFF3": "c_elegans.PRJNA13758.WS283.csq.gff3.gz (WormBase, cited)",
}


def msg(*a):
    print("[deposit]", *a, flush=True)


def read_audit():
    if not os.path.exists(AUDIT):
        msg("no audit yet; running scripts/audit_data_dependencies.py")
        subprocess.run([sys.executable,
                        os.path.join(ROOT, "scripts", "audit_data_dependencies.py")],
                       check=True, cwd=ROOT)
    keep, drop = [], []
    with open(AUDIT) as fh:
        next(fh)
        for line in fh:
            cls, path, nbytes, _ = line.rstrip("\n").split("\t", 3)
            (keep if cls.startswith("required") else drop).append((path, int(nbytes)))
    return keep, drop


def link(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        if os.path.samefile(src, dst):
            return 0
        os.remove(dst)
    try:
        os.link(src, dst)
    except OSError:                      # different volume, or a permission wall
        shutil.copy2(src, dst)
    return 1


def build(args):
    keep, drop = read_audit()
    os.makedirs(DEPOSIT, exist_ok=True)

    msg(f"linking {len(keep)} files under data/")
    linked = total = 0
    for p, n in keep:
        linked += link(os.path.join(ROOT, p), os.path.join(DEPOSIT, p))
        total += n

    msg("linking external inputs")
    ext_bytes = 0
    for src, rel, _var, _why in EXTERNAL:
        dst = os.path.join(DEPOSIT, "external", rel)
        if os.path.isdir(src):
            for f in sorted(os.listdir(src)):
                s = os.path.join(src, f)
                if os.path.isfile(s):
                    link(s, os.path.join(dst, f))
                    ext_bytes += os.path.getsize(s)
        elif os.path.exists(src):
            link(src, dst)
            ext_bytes += os.path.getsize(src)
        else:
            msg(f"  WARNING: absent, not deposited -- {src}")

    msg("copying code and prose")
    for t in COPY_TREES:
        s = os.path.join(ROOT, t)
        if os.path.isdir(s):
            shutil.copytree(s, os.path.join(DEPOSIT, t), dirs_exist_ok=True)
    for f in COPY_FILES:
        s = os.path.join(ROOT, f)
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(DEPOSIT, f))
    for d in ("plots", "plots/assets", "plots/diagnostics"):
        os.makedirs(os.path.join(DEPOSIT, d), exist_ok=True)
    src_assets = os.path.join(ROOT, "plots", "assets")
    if os.path.isdir(src_assets):
        shutil.copytree(src_assets, os.path.join(DEPOSIT, "plots", "assets"),
                        dirs_exist_ok=True)

    write_env()
    write_readme(len(keep), total, ext_bytes, len(drop),
                 sum(n for _, n in drop))
    write_manifest(checksums=not args.no_checksums)
    msg(f"built {os.path.relpath(DEPOSIT, ROOT)}/ -- "
        f"{total / 2**30:.2f} GB linked, {ext_bytes / 2**20:.0f} MB external, "
        f"{linked} new links")
    msg("disk cost of the link step is ~0; `du -sh deposit` counts shared bytes twice")


def write_env():
    lines = [
        "# Point the analysis at the files inside this deposit.",
        "#",
        "#   cd <this directory> && source env.sh",
        "#",
        "# Every external input the scripts read is reachable through one of these",
        "# variables; each script falls back to the original machine's path, which",
        "# is why sourcing this file matters.",
        "",
        'DEPOSIT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"',
        "",
        "# --- carried in this deposit ---------------------------------------",
        'export JX_PLOTS="$DEPOSIT_ROOT/external/jx_cross"',
        'export BAUGH_BOOT="$DEPOSIT_ROOT/external/baugh_wgs/'
        '20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData"',
        'export SIM_TRAITS="$DEPOSIT_ROOT/external/traits_with_validated_qtl"',
        'export XQTL_GMAP="$DEPOSIT_ROOT/external/xqtl_sims/geneticMapXQTLsnplist.rds"',
        "",
        "# --- cited, NOT carried: set these yourself once you have them ------",
        "# See DATA_AVAILABILITY.md and scripts/fetch_cendr_genotypes.sh.",
    ]
    for var, why in FETCHED.items():
        lines.append(f"# {var:<16} {why}")
        lines.append(f'# export {var}="/path/to/..."')
    lines.append("")
    lines.append('echo "deposit environment set for $DEPOSIT_ROOT"')
    with open(os.path.join(DEPOSIT, "env.sh"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_manifest(checksums=True):
    rows = []
    for dirpath, _, files in os.walk(DEPOSIT):
        for f in sorted(files):
            p = os.path.join(dirpath, f)
            rel = os.path.relpath(p, DEPOSIT)
            if rel in ("MANIFEST.tsv",):
                continue
            n = os.path.getsize(p)
            if checksums:
                h = hashlib.sha256()
                with open(p, "rb") as fh:
                    for blk in iter(lambda: fh.read(1 << 22), b""):
                        h.update(blk)
                digest = h.hexdigest()
            else:
                digest = ""
            rows.append((rel, n, digest))
    rows.sort()
    with open(os.path.join(DEPOSIT, "MANIFEST.tsv"), "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for rel, n, d in rows:
            fh.write(f"{rel}\t{n}\t{d}\n")
    msg(f"manifest: {len(rows)} files" + ("" if checksums else " (no checksums)"))


def write_readme(n_keep, b_keep, b_ext, n_drop, b_drop):
    from textwrap import dedent
    txt = dedent(f"""\
    # Data and code deposit

    Everything needed to reproduce the analyses in this manuscript, with one
    deliberate exception: three public third-party datasets are **cited rather
    than redistributed**, and a script here rebuilds what the analysis derived
    from them. That is the only thing you may have to fetch.

    ## Layout

    ```
    README.md              this file
    env.sh                 source it first -- points the scripts at this tree
    MANIFEST.tsv           every file, its size and its SHA-256
    scripts/               all analysis code
    supplemental_data/     the staged inputs the figures read
    data/                  the provenance tree: the files the scripts actually read
    external/              inputs that lived outside the original repository
    plots/                 where figures are written (assets/ ships pre-rendered)
    METHODS.txt            methods, with the caveats
    FIGURE_CAPTIONS.txt    working captions, with the numbers and their sources
    MANUSCRIPT_CAPTIONS.txt  the publication captions
    DATA_AVAILABILITY.md   what is here, what is cited, and why
    ```

    ## Running it

    Scripts resolve their inputs **relative to the working directory**, so run
    them from this directory, not from `scripts/`:

    ```bash
    cd <this directory>
    source env.sh                      # sets JX_PLOTS, BAUGH_BOOT, SIM_TRAITS, XQTL_GMAP
    Rscript scripts/Figure3_quad.R     # writes plots/Figure3_quad.{{pdf,png}}
    ```

    `source env.sh` is not optional for the staging scripts. Each script falls
    back to an absolute path on the machine the analysis was done on, so without
    the environment they will look somewhere that does not exist here and stop
    with a clear error rather than guessing.

    ## Two tiers, and only the first needs anything fetched

    **The figures build from `supplemental_data/` alone.** All 27 of them. This
    is the tier most readers want, it needs no external data, and it is the one
    verified by deleting `data/` and rebuilding:

    ```bash
    Rscript scripts/Figure1_pos1.R
    Rscript scripts/Figure2.R
    Rscript scripts/Figure3_quad.R
    Rscript scripts/Figure4_sid2.R
    ```

    **Rebuilding `supplemental_data/` from the provenance tree** is the second
    tier, and it is what `data/` and `external/` are for. Those staging scripts
    are the `make_*.R` and `make_*.py` ones. Some also need the cited datasets
    below.

    ## What is cited rather than deposited, and how to get it

    | dataset | what the analysis derives from it | how to rebuild |
    |---|---|---|
    | CaeNDR release **20210121** | the PLINK genotype panel, 2,917,997 markers over 540 isotypes | `scripts/fetch_cendr_genotypes.sh plink` |
    | CaeNDR release **20231213** | `bcsq.vcf.gz`, the haplotype-aware consequence annotation | `scripts/fetch_cendr_genotypes.sh bcsq` |
    | 207-strain expression + eQTL set | the expression check on candidate genes | two files, named in `DATA_AVAILABILITY.md` |
    | WormBase **WS283** GFF3 | gene models for the interval tables | one file, named in `DATA_AVAILABILITY.md` |

    These are other people's data, and they are large: the two CaeNDR releases
    are 8.3 GB each. The fetch script does not merely download them, it
    reproduces the exact derived files, using the original commands recovered
    from the PLINK logs and the VCF header. Its `csq` step is verified against
    the shipped annotation; its PLINK step asserts the per-chromosome variant
    and sample counts and stops on a mismatch.

    After fetching, export the matching variables -- `env.sh` lists them,
    commented out, with what each one wants.

    ## Software

    R with tidyverse, data.table, patchwork, ggtext, ggrepel, scales; Python 3
    with numpy, pandas, matplotlib, biopython; PLINK 1.9 **and** PLINK 2;
    bcftools; GEMMA 0.98.5. The versions used were PLINK v1.90b6.21,
    PLINK v2.00a3, bcftools 1.11. The structure renders in `plots/assets/` were
    made with PyMOL and are shipped pre-rendered, so Figure 4 builds without it.

    ## What is in `data/`, and what was left out

    `data/` holds **{n_keep:,} files, {b_keep / 2**30:.2f} GB** -- the ones a script
    actually reads, identified by `scripts/audit_data_dependencies.py`, which
    classifies every file by whether a script names its path, globs its
    directory, or mentions its name. A further **{n_drop:,} files
    ({b_drop / 2**30:.2f} GB)** in the original tree are read by nothing and are not
    deposited; most of that is exploratory structure modelling that never
    reached the manuscript, and one PDB from it that did is included.

    `external/` adds **{b_ext / 2**20:.0f} MB** that lived outside the original
    repository: the three JU2466 x XZ1516 contrast tables, the L1 bootstrap
    input, the validated-QTL trait tables, and the xQTL genetic map.

    ## Caveats worth reading before reusing the data

    These are stated where they matter, in `METHODS.txt` and
    `FIGURE_CAPTIONS.txt`, but three affect anyone reusing the tables:

    - **The hatching assays are one plate per strain per condition.** The
      intervals are Wilson binomial intervals on that plate's embryo count, so
      they describe counting uncertainty, not between-plate variability, and the
      Fisher tests compare plates with plate confounded with genotype.
    - **RNAi dose is not constant across figures** -- 50% for Figure 3C and 4A,
      25% for Figure 4B -- so hatching percentages are not comparable between
      them.
    - **The JU2466 x XZ1516 cross is incomplete**: no HT115 control pool, so
      every contrast is one RNAi condition against another and a peak can be
      driven by either side. It is reported as a diagnostic, not a result.

    ## Checking what you received

    ```bash
    awk -F'\\t' 'NR>1 {{print $3"  "$1}}' MANIFEST.tsv | shasum -a 256 -c
    ```
    """)
    with open(os.path.join(DEPOSIT, "README.md"), "w") as fh:
        fh.write(txt)


def archive(args):
    _, drop = read_audit()
    if not drop:
        msg("nothing to archive")
        return
    moved = freed = 0
    for p, n in drop:
        src = os.path.join(ROOT, p)
        if not os.path.exists(src):
            continue
        dst = os.path.join(ARCHIVE, p)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if args.dry_run:
            print(f"  would move {p}")
        else:
            shutil.move(src, dst)
        moved += 1
        freed += n
    verb = "would move" if args.dry_run else "moved"
    msg(f"{verb} {moved} files ({freed / 2**30:.2f} GB) to "
        f"{os.path.relpath(ARCHIVE, ROOT)}/")
    if not args.dry_run:
        msg("nothing is deleted. Verify a full rebuild, then delete that directory.")


def verify(args):
    script = args.script
    msg(f"running {script} inside the deposit")
    env = dict(os.environ)
    env.update({
        "JX_PLOTS": os.path.join(DEPOSIT, "external", "jx_cross"),
        "BAUGH_BOOT": os.path.join(
            DEPOSIT, "external", "baugh_wgs",
            "20220908_Baugh_BulkL1_Bootstrap_Input_flippedCommon_NAfix.RData"),
        "SIM_TRAITS": os.path.join(DEPOSIT, "external", "traits_with_validated_qtl"),
        "XQTL_GMAP": os.path.join(DEPOSIT, "external", "xqtl_sims",
                                  "geneticMapXQTLsnplist.rds"),
    })
    r = subprocess.run(["Rscript", script], cwd=DEPOSIT, env=env)
    msg("exit", r.returncode)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build");   b.add_argument("--no-checksums", action="store_true")
    a = sub.add_parser("archive"); a.add_argument("--dry-run", action="store_true")
    v = sub.add_parser("verify");  v.add_argument("--script", default="scripts/Figure3_quad.R")
    args = ap.parse_args()
    return {"build": build, "archive": archive, "verify": verify}[args.cmd](args) or 0


if __name__ == "__main__":
    sys.exit(main())
