#!/usr/bin/env python3
"""Every number in the manuscript draft must be one the code computed.

    python3 scripts/check_manuscript_numbers.py [draft.md ...]

Default target is MANUSCRIPT.md. Exits non-zero on any unattested number, so
the pre-push hook refuses the push.

WHY THIS EXISTS. Every number in the captions, the methods and the report is
either computed at knit time or asserted against the deposit. A manuscript
drafted anywhere else -- another Claude surface, a co-author's Word file, your
own memory -- carries those numbers as prose with no assertion, which is exactly
the failure mode the rest of this tooling exists to prevent: not a crash, but a
plausible sentence with a wrong number in it.

WHAT COUNTS AS ATTESTED. FIGURE_REPORT.md, and nothing else. It is regenerated
from supplemental_data/ on every knit, its tables are computed rather than
transcribed, and it is the document the manuscript is supposed to agree with.

  This checks DRIFT, not truth. A number wrong in FIGURE_REPORT.md will attest
  an identically wrong number in the draft. Nothing here substitutes for
  reviewing the report itself -- it only guarantees the two cannot diverge.

MATCHING. A draft number matches if some attested number equals it when rounded
to the draft's own precision, so "r = 0.997" matches a computed 0.99716 and
"RMSE 0.04" matches 0.0376. Percent forms are tried both ways, so "0.86%"
matches 0.0086. Exponent spellings are normalised: 1.4e-05, 1.4e-5, 1.4E-5 and
1.4 x 10^-5 are one number.

WHAT IS IGNORED, and why each is safe to ignore:
  strain names        JU1793, NIC256, wSZ203, CB4856 ... -- identifiers
  variant names       T96K, P153T, D34, N94A, V5L        -- identifiers
  gene names          sid-2, pos-1, mig-6                -- identifiers
  figure references   Figure 3, Figure S12, panel 4B     -- pointers
  years and citations 2026, (Li & Ji 2005)               -- bibliography
  chromosomes         I, II, III, IV, V, X               -- not numerals here
  bare 0 and 1        thresholds, ratios, "one of two"   -- uninformative
  small integers <=12 counts, panel numbers, replicates   -- too common to gate on

An unattested number that is legitimate goes in
manuscript_number_exceptions.txt, one per line as `number  # why`. That file is
the record of every number in the manuscript that the code does not produce,
which is a useful thing to be able to read on its own.
"""
import os
import re
import sys

ATTEST = "FIGURE_REPORT.md"
EXCEPTIONS = "manuscript_number_exceptions.txt"
# The draft lives in Google Docs. Export it (File > Download > Plain text, or
# Markdown) into manuscript/ and this picks it up; any path can also be passed
# explicitly. Nothing here needs the draft to be authored in the repository.
DEFAULT_GLOBS = ["manuscript/*.txt", "manuscript/*.md", "MANUSCRIPT.md"]

# identifiers and pointers, stripped before any number is read
STRIP = [
    r"\b(?:JU|NIC|ECA|CB|CX|ED|EG|MY|PB|PS|QG|QW|QX|RC|WN|XZ|DL|KR|LKC|GXW|BRC|AB|N2|wSZ)\d+[A-Za-z]?\b",
    r"\b[A-Z]\d+[A-Z]\b",                  # T96K, D34A
    r"\b[A-Z]\d+\b",                       # D34, H168, T96
    r"\b(?:sid|pos|mig|rde|sago|ppw)-\d+\b",
    r"\bFigures?\s+S?\d+[A-D]?(?:\s*(?:,|and)\s*S?\d+[A-D]?)*",
    r"\bpanels?\s+[A-E]\b",
    r"\b(?:19|20)\d{2}\b",                 # years
    r"\bBC[1-7]\b",
    r"\bchr(?:om)?\s*(?:I{1,3}V?|IV|V|X)\b",
    r"\bWS\d+\b", r"\bWBcel\d+\b", r"\bG5EEV9\b",
    r"\b[pP]\.?\s*=\s*NS\b",
]
NUM = re.compile(r"(?<![A-Za-z0-9_.])"
                 r"(-?\d+(?:,\d{3})*(?:\.\d+)?"
                 r"(?:\s*[eE]\s*[-+]?\d+|\s*[x×]\s*10\^?-?\d+)?)"
                 r"\s*(%?)")


def normalise(tok, pct):
    t = tok.replace(",", "").replace("−", "-").replace(" ", "")
    m = re.match(r"^(-?[\d.]+)[x×]10\^?(-?\d+)$", t)
    if m:
        t = f"{m.group(1)}e{m.group(2)}"
    try:
        v = float(t)
    except ValueError:
        return None
    return (v, bool(pct), t)


def extract(text):
    for pat in STRIP:
        text = re.sub(pat, " ", text)
    out = []
    for m in NUM.finditer(text):
        n = normalise(m.group(1), m.group(2))
        if n is None:
            continue
        v, pct, raw = n
        out.append((v, pct, raw, m.start()))
    return out


def decimals(raw):
    raw = raw.split("e")[0].split("E")[0]
    return len(raw.split(".")[1]) if "." in raw else 0


def attested_values(path):
    """Numbers the report states, plus the fraction a percent implies.

    Deliberately NOT v * 100 for every plain number: that made the pool so
    permissive that a fabricated 417 matched the eigen threshold 4.17 and the
    check passed a draft with invented figures in it.
    """
    vals = set()
    for v, pct, raw, _ in extract(open(path, encoding="utf-8").read()):
        vals.add(v)
        if pct:
            vals.add(v / 100.0)
    return vals


def matches(v, pct, dp, pool):
    """Match at the draft's own precision, so 0.997 accepts a computed 0.99716.

    A percent in the draft may also match the corresponding fraction, but a
    plain number is never rescaled -- rescaling everything is what let the first
    version of this check pass fabricated numbers.
    """
    cands = [v] + ([v / 100.0] if pct else [])
    for c in cands:
        for a in pool:
            try:
                if round(a, dp) == round(c, dp):
                    return True
                if pct and round(a, dp + 2) == round(c, dp + 2):
                    return True
            except (ValueError, OverflowError):
                continue
    return False


def load_exceptions(path):
    keep = set()
    if not os.path.exists(path):
        return keep
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if not line:
            continue
        try:
            keep.add(float(line.replace(",", "")))
        except ValueError:
            pass
    return keep


def main():
    import glob
    if sys.argv[1:]:
        present = [d for d in sys.argv[1:] if os.path.exists(d)]
    else:
        found = sorted(set(sum((glob.glob(g) for g in DEFAULT_GLOBS), [])))
        # manuscript/README.md documents the workflow and contains example
        # numbers; it is not a draft. Anything named README is skipped.
        present = [f for f in found
                   if not os.path.basename(f).upper().startswith("README")]
    if not present:
        print("  note  no exported draft under manuscript/; nothing to check")
        return 0
    if not os.path.exists(ATTEST):
        print(f"  FAIL  {ATTEST} is missing -- knit the report before checking")
        return 1

    pool = attested_values(ATTEST)
    allow = load_exceptions(EXCEPTIONS)
    fails = 0
    for d in present:
        text = open(d, encoding="utf-8").read()
        lines = text.splitlines()
        offsets, run = [], 0
        for ln in lines:
            offsets.append(run)
            run += len(ln) + 1
        bad = []
        seen = set()
        for v, pct, raw, pos in extract(text):
            # The first version skipped abs(v) <= 1, which discarded every
            # correlation, RMSE, p value and frequency in the manuscript -- i.e.
            # nearly everything worth checking. Only exact 0 and 1 and small
            # integer counts are uninformative enough to skip.
            if v in (0.0, 1.0) or (float(v).is_integer() and abs(v) <= 12):
                continue
            key = (v, pct)
            if key in seen:
                continue
            seen.add(key)
            if v in allow or (pct and v / 100.0 in allow):
                continue
            if matches(v, pct, decimals(raw), pool):
                continue
            lineno = max(i for i, o in enumerate(offsets) if o <= pos) + 1
            bad.append((lineno, raw + ("%" if pct else "")))
        if bad:
            fails += len(bad)
            print(f"  FAIL  {d}: {len(bad)} number(s) not produced by the code")
            for lineno, raw in bad[:25]:
                snippet = lines[lineno - 1].strip()
                if len(snippet) > 78:
                    snippet = snippet[:75] + "..."
                print(f"          line {lineno}: {raw}")
                print(f"            {snippet}")
            if len(bad) > 25:
                print(f"          ... and {len(bad) - 25} more")
        else:
            print(f"  ok    {d}: every number is one {ATTEST} attests")
    if fails:
        print(f"\n        Fix the draft, or record the number in {EXCEPTIONS}")
        print("        as `number  # why it is not computed here`.")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
