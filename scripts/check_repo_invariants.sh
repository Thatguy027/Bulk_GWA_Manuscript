#!/usr/bin/env bash
# Mechanical invariant checks for this repository.
#
#   bash scripts/check_repo_invariants.sh          # fast checks only
#   bash scripts/check_repo_invariants.sh --full   # adds isolation and determinism
#
# The fast checks take seconds and are safe to wire into a hook or a pre-push.
# --full rebuilds every figure twice and once more with data/ hidden, which
# takes minutes and is what you want on a schedule rather than on every edit.
#
# Exit status is the number of failed checks, so a hook can gate on it.
set -uo pipefail
cd "$(dirname "$0")/.."
FAIL=0
ok()   { printf '  ok    %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; FAIL=$((FAIL+1)); }
note() { printf '        %s\n' "$1"; }

echo "== 1. no figure script reads data/ =="
# builders are allowed to: they read the archive and write the deposit table
# Only the FIGURE scripts have to run from a clone. Deposit builders, prep
# scripts and diagnostics legitimately read the Dryad archive; they are what
# turn it into the small tables the figure scripts read.
offenders=$(grep -ln '"data/' scripts/Figure*.R scripts/SUPP_FIG*.R 2>/dev/null || true)
if [ -z "$offenders" ]; then ok "no figure script references data/"
else bad "figure scripts reading data/ would break a fresh clone:"; echo "$offenders" | sed 's/^/        /'; fi
archive_readers=$(grep -ln '"data/' scripts/*.R scripts/*.py 2>/dev/null \
  | grep -v -E 'scripts/(Figure|SUPP_FIG)' | wc -l | tr -d ' ')
note "$archive_readers non-figure scripts read the archive, which is expected"

echo "== 2. randomness is seeded =="
unseeded=$(grep -ln 'geom_jitter\|position_jitter\|geom_text_repel\|geom_label_repel' scripts/*.R 2>/dev/null | while read -r f; do
  # a file is suspect if it uses a jitter/repel geom but never sets a seed
  grep -q 'seed *=' "$f" || echo "$f"
done)
if [ -z "$unseeded" ]; then ok "every jitter/repel geom sets a seed"
else bad "jitter or repel without a seed (figure will not reproduce):"; echo "$unseeded" | sed 's/^/        /'; fi

echo "== 3-4. figure lists agree with plots/ and with the captions =="
if Rscript scripts/check_figure_lists.R; then :; else
  FAIL=$((FAIL + $?)); fi

echo "== 5. scripts that self-assert their pinned numbers still pass =="
for s in scripts/SUPP_FIG_XX_dilution_validation.R; do
  if Rscript "$s" >/tmp/pin_$(basename "$s").log 2>&1; then ok "$(basename "$s") pins agree"
  else bad "$(basename "$s") failed — pins stale or a real error"
       tail -4 /tmp/pin_$(basename "$s").log | sed 's/^/        /'; fi
done

if [ "${1:-}" = "--full" ]; then
  echo "== 6. determinism: every figure byte-identical across two runs =="
  mkdir -p /tmp/det && rm -f /tmp/det/*
  for f in plots/*.png; do cp "$f" "/tmp/det/$(basename "$f")"; done
  for s in scripts/Figure*.R scripts/SUPP_FIG*.R; do Rscript "$s" >/dev/null 2>&1; done
  n=0
  for f in plots/*.png; do
    cmp -s "$f" "/tmp/det/$(basename "$f")" || { bad "not deterministic: $(basename "$f")"; n=$((n+1)); }
  done
  [ "$n" = 0 ] && ok "all figures byte-identical across runs"

  echo "== 7. isolation: every figure builds with data/ absent =="
  if [ -d data ]; then
    mv data /tmp/data_hidden_check
    n=0
    for s in scripts/Figure*.R scripts/SUPP_FIG*.R; do
      Rscript "$s" >/dev/null 2>&1 || { bad "needs data/: $(basename "$s")"; n=$((n+1)); }
    done
    mv /tmp/data_hidden_check data
    [ "$n" = 0 ] && ok "all figures build from supplemental_data/ alone"
  else note "data/ not present; isolation already implied"; fi
fi

echo
echo "$FAIL check(s) failed"
exit "$FAIL"
