## Which design matrix did the archived Baugh deconvolution use? ------------
##
##   Rscript scripts/DIAG_baugh_design_check.R
##     -> plots/diagnostics/baugh_design_check.tsv
##
## THE QUESTION. NNLS can be set up two ways here. Plain G writes one equation
## per marker, about the alternate allele: the alternate frequency at that
## marker is the summed frequency of the strains carrying it. The
## complement-stacked design rbind(G, 1 - G) writes a second equation per
## marker about the reference allele. The two are not redundant -- they would
## be only if the weights were constrained to sum to the depth, and NNLS
## normalises after the fit, not during -- so they are different estimators. In
## simulation the difference is worth about 0.14 of r-squared at 1x.
##
## WHY IT MATTERED. scripts/legacy/haploReg_original.R fits a point estimate on
## the stacked design and its bootstrap on the plain one, so METHODS.txt
## carried an open item: if the Baugh run inherited that, the published
## intervals would belong to a different estimator from the published point
## estimates. This settles it from the archived input rather than by reading
## the legacy code, which is mostly simulation.
##
## WHAT IT SHOWS. Two things, either sufficient on its own.
##
##   1  Refitting the point estimate on PLAIN G reproduces the archived
##      frequencies exactly -- identical to double precision across all 2,346
##      strain-by-sample cells, not merely close.
##   2  The stacked design is not constructible from the archived input at all.
##      2024bootstrapINPUT.Rdata stores genotypes and ALTERNATE counts only;
##      rbind(G, 1 - G) needs a reference-count block of the same shape, and
##      none was stored.
##
## So the Baugh point estimate and its bootstrap use the same design, plain G,
## and the inconsistency in the legacy script belongs to its simulation arm.
##
## Reads the archive (data/baugh), which is not deposited, so this does not run
## from a clone -- it is a diagnostic, not a figure script. Writes only to
## plots/diagnostics/.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

IN   <- "data/baugh/2024bootstrapINPUT.Rdata"
ARCH <- "data/baugh/2024_processedBOOTs_with_MIP.RData"
DIAG <- "plots/diagnostics"
dir.create(DIAG, recursive = TRUE, showWarnings = FALSE)
for (f in c(IN, ARCH))
  if (!file.exists(f)) stop("needs the archive: ", f, call. = FALSE)

e <- new.env(); load(IN, e)
stopifnot(length(e$flipped_bootstrap_input) == 2)
gt0 <- e$flipped_bootstrap_input[[1]]
ct0 <- e$flipped_bootstrap_input[[2]]
cat(sprintf("archived input: genotypes %s, alternate counts %s\n",
            paste(dim(gt0), collapse = " x "), paste(dim(ct0), collapse = " x ")))
## plain, not stacked: one row per marker. Stacked would be twice as tall.
stopifnot(nrow(gt0) == nrow(ct0), ncol(gt0) == 102, ncol(ct0) == 23)

## markers with any missing genotype are dropped, as baugh_frequencies() does
keep <- which(apply(gt0, 1, function(x) sum(is.na(x))) == 0)
cat(sprintf("markers with no missing call: %s of %s\n",
            format(length(keep), big.mark = ","),
            format(nrow(gt0), big.mark = ",")))
gt <- gt0[keep, ]; a_ct <- ct0[keep, ]

GGp  <- crossprod(gt)
Gy   <- crossprod(gt, a_ct)
pred <- apply(Gy, 2, function(x)
  as.vector(RcppML::nnls(GGp, matrix(x), fast_nnls = TRUE)))
pred <- apply(pred, 2, function(x) x / sum(x))
rownames(pred) <- sub("_.*$", "", colnames(GGp))

a <- new.env(); load(ARCH, a)
arch <- as.data.frame(a$wgs_mip_results)
j <- arch %>% select(strain, sample, frq) %>% filter(is.finite(frq)) %>%
  inner_join(as.data.frame(pred) %>% tibble::rownames_to_column("strain") %>%
               pivot_longer(-strain, names_to = "sample", values_to = "refit"),
             by = c("strain", "sample"))

maxabs <- max(abs(j$frq - j$refit))
exact  <- identical(j$frq, j$refit)
cat(sprintf("\nrefit on PLAIN G against the archived frequencies:\n"))
cat(sprintf("  cells compared        %d\n", nrow(j)))
cat(sprintf("  max |refit - archive| %s\n", format(maxabs)))
cat(sprintf("  identical to double precision %s\n", exact))
cat("\nthe stacked design needs a reference-count block of the same shape;\n")
cat("the archived input stores genotypes and alternate counts only, so it\n")
cat("cannot be built from what the original run saved.\n")

## PINNED. These are what METHODS.txt states.
stopifnot(nrow(j) == 2346, maxabs == 0, exact,
          nrow(gt0) == 1237106L, length(keep) == 1237106L)
cat("\npins agree\n")

readr::write_tsv(tibble::tibble(
  design = "plain G", cells = nrow(j), max_abs_diff = maxabs,
  identical_double = exact, markers = nrow(gt0),
  markers_complete = length(keep), stacked_constructible = FALSE),
  file.path(DIAG, "baugh_design_check.tsv"))
cat("wrote ", file.path(DIAG, "baugh_design_check.tsv"), "\n", sep = "")
