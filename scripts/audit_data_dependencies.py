#!/usr/bin/env python3
"""Which files under data/ does the analysis actually need?

    python3 scripts/audit_data_dependencies.py
      -> DATA_DEPENDENCY_AUDIT.tsv   one row per file under data/

Answers the question the Dryad deposit turns on: of the ~9,700 files in the
provenance tree, which are read by a script, which are read only by a
superseded script under scripts/legacy/, and which are read by nothing.

HOW A FILE IS CLASSIFIED, and why it errs toward keeping
--------------------------------------------------------
Three passes, each less precise than the one before, and a file is kept if ANY
of them claims it:

  1  resolved path   a script names the path, or names a variable that resolves
                     to it. Variable definitions are followed to a fixed point,
                     so PLOTS <- file.path(EXPORT, "plot_data") resolves.
  2  covering dir    the file sits under a directory a script names. A script
                     that globs a directory reads an unknown subset of it, so
                     the whole directory is kept.
  3  basename        the file's own name appears anywhere in any script. This
                     catches paths built in ways pass 1 cannot follow, at the
                     cost of false positives on common names.

Under-inclusion breaks reproduction and over-inclusion only makes the deposit
larger, so the passes are deliberately generous. Pass 3 in particular is a
NET, not a claim: a file it flags is reported as `required_basename` so it can
be looked at rather than silently trusted.

WHAT THIS DOES NOT SEE
Absolute paths outside the repository -- the CaeNDR releases under
~/UCLA/Genomics_Data and the earlier bulkGWAS project trees -- are reported
separately at the end, because no deposit built only from data/ contains them.
"""
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "DATA_DEPENDENCY_AUDIT.tsv")
SCRIPT_DIRS = ("scripts", "scripts/legacy")
EXT = (".R", ".r", ".py", ".sh", ".Rmd")

ASSIGN = re.compile(r'^\s*([A-Za-z_][A-Za-z0-9_.]*)\s*(?:<-|=)\s*(.+?)\s*$')
LIT = re.compile(r'^["\']([^"\'\n]+)["\']$')
JOIN = re.compile(r'^(?:file\.path|os\.path\.join)\(\s*(.+?)\s*\)$')
PASTE = re.compile(r'^paste0?\(\s*(.+?)\s*\)$')
TOKEN = re.compile(r'["\']([^"\'\n]*)["\']|([A-Za-z_][A-Za-z0-9_.]*)')
GLOB = re.compile(r'(list\.files|Sys\.glob|list\.dirs|glob\.glob)\s*\(\s*([A-Za-z_][A-Za-z0-9_.]*)')
EXTERNAL = re.compile(r'["\'](/Users/[^"\'\n]+)["\']')

# A bare "data" appears inside a console message in check_repo_invariants.sh and
# is not a path; taken as one it covers the whole tree and the audit says
# everything is required.
NOT_A_PATH = {"data", "data/"}


def scripts():
    out = []
    for d in SCRIPT_DIRS:
        p = os.path.join(ROOT, d)
        if not os.path.isdir(p):
            continue
        for f in sorted(os.listdir(p)):
            q = os.path.join(p, f)
            if f.endswith(EXT) and os.path.isfile(q):
                out.append(os.path.relpath(q, ROOT))
    return out


def resolve_tokens(expr, env, sep):
    parts = []
    for m in TOKEN.finditer(expr):
        if m.group(1) is not None:
            parts.append(m.group(1))
            continue
        name = m.group(2)
        if name in env:
            parts.append(env[name])
        elif name in ("file", "path", "os", "join", "paste0", "paste", "sep"):
            continue
        else:
            return None
    return sep.join(parts) if parts else None


def script_env(txt):
    """Variable -> string value, iterated so chained definitions resolve."""
    env = {}
    for _ in range(6):
        grew = False
        for ln in txt.splitlines():
            m = ASSIGN.match(ln)
            if not m or m.group(1) in env:
                continue
            rhs = m.group(2).rstrip(",")
            val = None
            if LIT.match(rhs):
                val = LIT.match(rhs).group(1)
            elif JOIN.match(rhs):
                val = resolve_tokens(JOIN.match(rhs).group(1), env, "/")
            elif PASTE.match(rhs):
                val = resolve_tokens(PASTE.match(rhs).group(1), env, "")
            if val:
                env[m.group(1)] = val
                grew = True
        if not grew:
            break
    return env


def main():
    os.chdir(ROOT)
    paths = defaultdict(set)      # data path -> scripts naming it
    globbed = defaultdict(set)
    external = defaultdict(set)
    corpus = []

    for s in scripts():
        txt = open(s, encoding="utf-8", errors="replace").read()
        corpus.append(txt)
        env = script_env(txt)

        for m in re.finditer(r'["\']([^"\'\n]*data/[^"\'\n]*)["\']', txt):
            p = m.group(1)
            if p.startswith("data/") and p not in NOT_A_PATH:
                paths[os.path.normpath(p)].add(s)
        for v in env.values():
            if v.startswith("data/") and v not in NOT_A_PATH:
                paths[os.path.normpath(v)].add(s)
        for rx, sep in ((r'(?:file\.path|os\.path\.join)\(([^()]*)\)', "/"),
                        (r'paste0\(([^()]*)\)', "")):
            for m in re.finditer(rx, txt):
                val = resolve_tokens(m.group(1), env, sep)
                if val and val.startswith("data/") and val not in NOT_A_PATH:
                    paths[os.path.normpath(val)].add(s)
        for m in GLOB.finditer(txt):
            v = env.get(m.group(2), "")
            if v.startswith("data/"):
                globbed[os.path.normpath(v)].add(s)
        for m in EXTERNAL.finditer(txt):
            external[m.group(1)].add(s)

    corpus = "\n".join(corpus)
    live = {p for p, sc in paths.items()
            if os.path.exists(p) and not all(x.startswith("scripts/legacy/") for x in sc)}
    legacy = {p for p, sc in paths.items()
              if os.path.exists(p) and all(x.startswith("scripts/legacy/") for x in sc)}

    def under(f, roots):
        return any(f == r or f.startswith(r + os.sep) for r in roots)

    rows, tally = [], defaultdict(lambda: [0, 0])
    for dirpath, _, filenames in os.walk("data"):
        for fn in filenames:
            f = os.path.join(dirpath, fn)
            try:
                sz = os.path.getsize(f)
            except OSError:
                continue
            if under(f, live):
                lab = "required"
            elif len(fn) >= 5 and fn != ".DS_Store" and fn in corpus:
                lab = "required_basename"
            elif under(f, legacy):
                lab = "legacy_only"
            else:
                lab = "unreferenced"
            readers = sorted({s for p, sc in paths.items() if under(f, {p}) for s in sc})
            rows.append((lab, f, sz, ";".join(readers)))
            tally[lab][0] += 1
            tally[lab][1] += sz

    rows.sort(key=lambda r: (r[0], -r[2]))
    with open(OUT, "w") as fh:
        fh.write("class\tpath\tbytes\tread_by\n")
        for lab, f, sz, who in rows:
            fh.write(f"{lab}\t{f}\t{sz}\t{who}\n")

    print(f"wrote {os.path.relpath(OUT, ROOT)}  ({len(rows)} files under data/)\n")
    print(f"{'class':20} {'files':>7} {'GB':>8}")
    for lab in ("required", "required_basename", "legacy_only", "unreferenced"):
        n, b = tally[lab]
        print(f"{lab:20} {n:7d} {b / 2**30:8.2f}")
    keep = sum(tally[k][1] for k in ("required", "required_basename"))
    drop = sum(tally[k][1] for k in ("legacy_only", "unreferenced"))
    print(f"\nkeep {keep / 2**30:.2f} GB, archivable {drop / 2**30:.2f} GB")

    if globbed:
        print("\ndirectories read by glob (whole directory kept):")
        for d in sorted(globbed):
            print(f"  {d}")
    print("\ninputs OUTSIDE the repository, which no data/ deposit contains:")
    seen = set()
    for p in sorted(external):
        base = p
        for known in ("/Users/Stefan/UCLA/Genomics_Data/CeNDR/20231213",
                      "/Users/Stefan/UCLA/Genomics_Data/CeNDR/expression",
                      "/Users/Stefan/UCLA/Genomics_Data/Annotations"):
            if p.startswith(known):
                base = known
        if base in seen:
            continue
        seen.add(base)
        if os.path.exists(base):
            sz = subprocess.run(["du", "-sh", base], capture_output=True,
                                text=True).stdout.split("\t")[0]
        else:
            sz = "ABSENT"
        print(f"  {sz:>8}  {base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
