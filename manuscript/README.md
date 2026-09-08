# The manuscript draft

The draft is written in **Google Docs**. Nothing here asks you to author it in
the repository — this directory holds periodic *exports*, so the prose can be
checked against the numbers the code produces and versioned alongside them.

## The loop

1. Write in Google Docs as normal.
2. When you want a check: **File → Download → Plain text (.txt)** (Markdown
   works too) and save it in here. Keep the same filename each time so the
   diffs are readable — `draft.txt` is fine.
3. Run the check:

   ```sh
   python3 scripts/check_manuscript_numbers.py
   ```

   or just push: the pre-push hook runs it along with the other invariants.

Every number in the export must be one `FIGURE_REPORT.md` states. Anything else
is reported with its line and the sentence it sits in.

## Why this direction

The check compares the draft against the **current** report. So when you tweak
a figure and a number moves, the stale sentence in the draft is what gets
flagged — which is the direction you want. You do not need to re-export after
changing a figure; you need to re-export after changing the *prose*.

## When the check flags something legitimate

Plenty of numbers in a real manuscript are not computed here — reagent
concentrations, incubation times, catalogue numbers, a co-author's institution
postcode. Record those in `manuscript_number_exceptions.txt` at the repository
root, one per line with a reason:

```
37       # incubation temperature, degrees C
1000     # dilution factor for the antibody stock
```

That file then reads as a complete list of every number in the manuscript the
code does not produce, which is worth being able to review on its own.

## What this does NOT do

It checks **drift**, not truth. A number that is wrong in `FIGURE_REPORT.md`
will happily attest an identically wrong number in the draft. Reviewing the
report itself is still the thing that catches a genuine error; this only
guarantees the two cannot silently diverge.

It also cannot see your Google Doc. If you have not exported recently, the
check passes on a stale file and tells you nothing. The export is the weak
link, and it is manual on purpose — an automated pull would need Docs
credentials in the repository, which is worse.
