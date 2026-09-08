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
say(!length(no_cap),
    if (length(no_cap)) paste0("no caption entry for: ", paste(no_cap, collapse = ", "))
    else "every curated figure appears in FIGURE_CAPTIONS.txt")

## figures present in plots/ but not in the curated list -- not a failure, but
## worth naming, since a stray figure is how a set silently grows
stray <- setdiff(sub("\\.png$", "", basename(Sys.glob("plots/*.png"))), figs)
if (length(stray))
  cat("        note: in plots/ but not curated: ",
      paste(stray, collapse = ", "), "\n", sep = "")

quit(status = fail)
