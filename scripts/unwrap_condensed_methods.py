#!/usr/bin/env python3
"""Unwrap METHODS_CONDENSED.txt so each paragraph is one line.

    python3 scripts/unwrap_condensed_methods.py

Hard-wrapped text pastes badly into a word processor -- every line break becomes
a break mid-sentence -- so the condensed Methods keeps one line per paragraph.
This is the tool that does it, and it is idempotent: run it again after editing
and it re-joins anything that got re-wrapped.

WHAT IS LEFT ALONE. Section rules, section titles, blank lines, and any line
starting with two spaces, which is how the tables are marked (the repair
templates, the threshold table, the RNAi doses). Those stay exactly as laid out.

It asserts the result is word-for-word identical to the input, so it can only
change line breaks and never content.
"""

import re, sys

src = "METHODS_CONDENSED.txt"
original = open(src).read()
lines = original.split("\n")

def is_rule(l):   return bool(re.fullmatch(r"-{60,}|={60,}", l.strip()))
def is_eq(l):     return bool(re.fullmatch(r"={60,}", l.strip()))
def is_indented(l): return l.startswith("  ") and l.strip()

# the "====" title block at the top is a header, not prose: keep it verbatim
n_eq = [i for i, l in enumerate(lines) if is_eq(l)]
title_end = n_eq[1] if len(n_eq) >= 2 else -1

out, buf = [], []
def flush():
    if buf:
        out.append(" ".join(x.strip() for x in buf))
        buf.clear()

for i, l in enumerate(lines):
    if i <= title_end:           # title block: verbatim
        flush(); out.append(l)
    elif not l.strip():          # blank: paragraph break
        flush(); out.append("")
    elif is_rule(l):             # section rule
        flush(); out.append(l)
    elif is_indented(l):         # table / list block: verbatim
        flush(); out.append(l)
    else:
        buf.append(l)
flush()

# collapse runs of >2 blank lines to 1
txt = "\n".join(out)
txt = re.sub(r"\n{3,}", "\n\n", txt)
open(src, "w").write(txt.rstrip("\n") + "\n")

# content must be preserved word for word
before = original.split()
after = txt.split()
if before != after:
    import difflib
    d = [x for x in difflib.unified_diff(before, after, lineterm="", n=0)][:20]
    sys.exit("WORDS CHANGED:\n" + "\n".join(d))
print(f"unwrapped: word-for-word identical ({len(after)} words)")
