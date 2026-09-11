## Consistency of the figure lists, for check_repo_invariants.sh -------------
## Parses FIGS and LABELS out of FIGURE_REPORT.Rmd and checks that they agree
## with each other, with plots/, and with FIGURE_CAPTIONS.txt. Exits non-zero
## on any mismatch so the shell caller can count failures.
s <- readLines("FIGURE_REPORT.Rmd", warn = FALSE)

vec <- function(tag) {
  i <- grep(tag, s)
  if (!length(i)) return(character(0))
  i <- i[1]
  j <- i
  depth <- 0L
  repeat {                                   # walk to the balanced close paren
    depth <- depth + lengths(regmatches(s[j], gregexpr("\\(", s[j]))) -
                     lengths(regmatches(s[j], gregexpr("\\)", s[j])))
    if (depth <= 0L || j >= length(s)) break
    j <- j + 1L
  }
  blk <- paste(s[i:j], collapse = " ")
  gsub('"', "", unlist(regmatches(blk, gregexpr('"[^"]+"', blk))))
}

figs <- vec("^FIGS *<- *c\\(")
labs_all <- vec("^LABELS *<- *c\\(")
## LABELS is name = "Figure N"; the quoted halves are the labels, the names are
## the stems, so pull the stems from the block's `stem =` tokens instead
i <- grep("^LABELS *<- *c\\(", s)[1]
j <- i; depth <- 0L
repeat {
  depth <- depth + lengths(regmatches(s[j], gregexpr("\\(", s[j]))) -
                   lengths(regmatches(s[j], gregexpr("\\)", s[j])))
  if (depth <= 0L || j >= length(s)) break
  j <- j + 1L
}
blk <- paste(s[i:j], collapse = " ")
lab_stems <- gsub("[ =]", "", unlist(regmatches(blk,
  gregexpr("[A-Za-z0-9_]+ *=", blk))))

fail <- 0L
say <- function(ok, msg) {
  cat(if (ok) "  ok    " else "  FAIL  ", msg, "\n", sep = "")
  if (!ok) fail <<- fail + 1L
}

say(length(figs) > 0, paste0("parsed FIGS: ", length(figs), " figures"))
miss_png <- figs[!file.exists(file.path("plots", paste0(figs, ".png")))]
say(!length(miss_png),
    if (length(miss_png)) paste0("listed but no PNG: ", paste(miss_png, collapse = ", "))
    else paste0("all ", length(figs), " listed figures have a PNG in plots/"))

say(setequal(figs, lab_stems),
    if (setequal(figs, lab_stems)) "FIGS and LABELS cover the same figures"
    else paste0("FIGS/LABELS disagree: ",
                paste(union(setdiff(figs, lab_stems), setdiff(lab_stems, figs)),
                      collapse = ", ")))

say(!anyDuplicated(lab_stems) &&
    length(unique(gsub("[^0-9S]", "", lab_stems))) == length(lab_stems) ||
    !anyDuplicated(lab_stems),
    "no duplicated figure stems in LABELS")

caps <- paste(readLines("FIGURE_CAPTIONS.txt", warn = FALSE), collapse = "\n")
no_cap <- figs[!vapply(figs, grepl, logical(1), x = caps, fixed = TRUE)]

## FIGURE_CAPTIONS.txt also opens with its own roster of the curated set, which
## a caption further down does not satisfy -- that block had 20 entries under a
## "twenty-three" heading, so check the block itself, not just the whole file
roster <- sub(".*CURATED SET[^\n]*\n", "",
              sub("Captions below cover those.*", "", caps))
in_roster <- vapply(figs, function(f)
  grepl(paste0("(^|\n)  ", f, " *(\n|$)"), roster), logical(1))
say(all(in_roster),
    if (all(in_roster))
      paste0("FIGURE_CAPTIONS.txt's curated-set roster lists all ",
             length(figs))
    else paste0("missing from the FIGURE_CAPTIONS.txt roster: ",
                paste(figs[!in_roster], collapse = ", ")))
say(!length(no_cap),
    if (length(no_cap)) paste0("no caption entry for: ", paste(no_cap, collapse = ", "))
    else "every curated figure appears in FIGURE_CAPTIONS.txt")

## ---------------------------------------------------------------------------
## The curated count is also stated in prose in three documents, and in a table
## in README.md. Those had drifted three ways -- README twenty-one, the deposit
## overview eighteen, plots/ twenty-three -- because nothing checked them, and
## the deposit's self-containment claim is stated AS a verified count, so a
## stale number there attests something narrower than what is in the tree.
WORD <- c("eighteen", "nineteen", "twenty", "twenty-one", "twenty-two",
          "twenty-three", "twenty-four", "twenty-five", "twenty-six",
          "twenty-seven", "twenty-eight", "twenty-nine", "thirty")
names(WORD) <- 18:30
want <- WORD[[as.character(length(figs))]]

## README's curated table: one row per figure, and no figure missing from it
readme <- readLines("README.md", warn = FALSE)
tbl <- gsub("^\\| `([A-Za-z0-9_]+)` \\| `[A-Za-z0-9_.]+` \\|.*$", "\\1",
            grep("^\\| `[A-Za-z0-9_]+` \\| `[A-Za-z0-9_.]+` \\|", readme,
                 value = TRUE))
say(setequal(tbl, figs),
    if (setequal(tbl, figs))
      paste0("README's curated table lists the same ", length(figs), " figures")
    else paste0("README table vs FIGS: ",
                paste(union(setdiff(figs, tbl), setdiff(tbl, figs)),
                      collapse = ", ")))

## the spelled-out count, in every document that states it. The (?!-) is
## load-bearing: "twenty" is a word-boundary match inside "twenty-three".
has_word <- function(txt, w)
  grepl(paste0("\\b", w, "\\b(?!-)"), txt, perl = TRUE,
        ignore.case = TRUE)

## The check is deliberately blunt: ANY count word in these files must be the
## current one. That is right today because every occurrence refers to the
## curated set; if one of these documents ever needs one of these words for a
## different quantity, write the digit or rephrase rather than loosening this.
## all five documents that state the count in words. It had drifted FIVE ways:
## README twenty-one, DATA_AVAILABILITY twenty-one, the deposit overview
## eighteen, FIGURE_CAPTIONS twenty, FIGURE_REPORT.Rmd twenty-one.
for (f in c("README.md", "DATA_AVAILABILITY.md",
            "supplemental_data/SUPPLEMENTAL_DATA_OVERVIEW.md",
            "FIGURE_CAPTIONS.txt", "FIGURE_REPORT.Rmd")) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  stale <- setdiff(WORD[vapply(WORD, has_word, logical(1), txt = txt)], want)
  say(!length(stale) && has_word(txt, want),
      if (length(stale))
        paste0(f, " states a stale figure count: ", paste(stale, collapse = ", "),
               " (should be ", want, ")")
      else if (!has_word(txt, want))
        paste0(f, " never states the figure count (expected \"", want, "\")")
      else paste0(f, " states the count as ", want))
}

## figures present in plots/ but not in the curated list -- not a failure, but
## worth naming, since a stray figure is how a set silently grows
stray <- setdiff(sub("\\.png$", "", basename(Sys.glob("plots/*.png"))), figs)
if (length(stray))
  cat("        note: in plots/ but not curated: ",
      paste(stray, collapse = ", "), "\n", sep = "")

quit(status = fail)
