## Supplement -- the N2 residue-96 swap across a pos-1 dose series -----------
##
##   Rscript scripts/SUPP_FIG_XX_n2_swap_dose.R
##     -> plots/SUPP_FIG_XX_n2_swap_dose.{pdf,png}
##
##   A  every genotype at every dose
##   B  the 25% dose alone, which is the only dose with dynamic range
##
## See scripts/n2_swap_panels.R for the experiment, the dose table and the
## caveats. The short version: 96K makes N2 more sensitive to pos-1 RNAi, the
## effect is confined to the 25% dose because every other dose sits at the
## ceiling or the floor, and two independently edited lines agree.
## ---------------------------------------------------------------------------

source("scripts/n2_swap_panels.R")
suppressPackageStartupMessages(library(patchwork))

OUT <- "plots"
panel_title <- function(letter, txt = NULL) {
  lt <- paste0("<span style='font-size:13pt;color:#111111'>**", letter,
               "**</span>")
  if (is.null(txt)) lt else paste0(lt, " ", txt)
}
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_markdown(size = base_size + 0.5),
          plot.subtitle = element_markdown(size = base_size - 2.5,
                                           colour = "grey30"),
          plot.title.position = "plot",
          legend.key.size = grid::unit(9, "pt"))
}
wrap_md <- function(txt, width = 78)
  paste(strwrap(txt, width = width), collapse = "<br>")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

d <- n2_swap_data()
tt <- n2_swap_tests(d)
msg("rows ", nrow(d), " | strains ", n_distinct(d$strain), " | doses ",
    paste(sort(unique(d$dose)), collapse = "/"))

## pooled comparison per dose, for the caption and for panel A's annotation
pool <- d %>% group_by(dose, allele) %>%
  summarise(h = sum(hatched), n = sum(n_plated), .groups = "drop") %>%
  pivot_wider(names_from = allele, values_from = c(h, n)) %>%
  rowwise() %>%
  mutate(p_96T = h_96T / n_96T, p_96K = h_96K / n_96K,
         fisher_p = fisher.test(matrix(c(h_96K, n_96K - h_96K,
                                         h_96T, n_96T - h_96T),
                                       nrow = 2, byrow = TRUE))$p.value) %>%
  ungroup()

cat("\n== pooled 96K vs N2[96T], per dose ==\n")
print(as.data.frame(pool %>% transmute(dose, n_96T, n_96K,
        `96T` = round(p_96T, 3), `96K` = round(p_96K, 3),
        difference = round(p_96T - p_96K, 3),
        fisher_p = signif(fisher_p, 3))), row.names = FALSE)
cat("\n== each edited line vs N2, per dose ==\n")
print(as.data.frame(tt %>% transmute(dose, strain, hatched = round(p, 3),
        n = n_plated, fisher_p = signif(fisher_p, 3)) %>%
        arrange(dose, strain)), row.names = FALSE)

## --- A: the whole dose series ---------------------------------------------
sig <- pool %>% filter(fisher_p < 0.001)
pA <- ggplot(d, aes(dose, p, colour = allele, group = line)) +
  geom_line(aes(linetype = strain), linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 2.4, linewidth = 0.4) +
  geom_point(size = 2.1) +
  ## mark the one dose that separates the genotypes
  geom_richtext(data = sig %>% mutate(y = 0.62),
                aes(dose, y, label = paste0("difference ",
                    sprintf("%+.2f", p_96T - p_96K), "<br>",
                    fmt_p(fisher_p))),
                inherit.aes = FALSE, size = 2.7, colour = "grey20",
                fill = "white", label.color = NA,
                label.padding = grid::unit(c(1, 2, 1, 2), "pt")) +
  scale_colour_manual(values = c(`96T` = COL_N2, `96K` = "#B0562A"),
                      name = "*sid-2* residue 96") +
  scale_linetype_manual(values = c(N2 = "solid", wSZ203 = "22",
                                   wSZ204 = "42"), name = NULL) +
  scale_x_continuous(breaks = sort(unique(d$dose)),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.02), expand = expansion(mult = c(0, 0))) +
  labs(x = "*pos-1* RNAi bacteria in the lawn", y = "Embryos hatched",
       title = panel_title("A", "**Only the 25% dose has any dynamic range**"),
       subtitle = wrap_md(paste0(
         "N2 carries 96T; wSZ203 and wSZ204 are independent lines edited to ",
         "96K. At 0% every genotype hatches and at 50% and above every ",
         "genotype is dead, so those doses cannot report a difference. One ",
         "plate per strain per dose, Wilson 95% intervals."))) +
  theme_pub() +
  theme(axis.title.x = element_markdown(),
        legend.title = element_markdown(size = 9),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(t = -2),
        legend.spacing.y = grid::unit(6, "pt"),
        legend.box.just = "left")

## --- B: the 25% dose alone ------------------------------------------------
d25 <- d %>% filter(dose == FOCAL_DOSE)
t25 <- tt %>% filter(dose == FOCAL_DOSE)
Y <- 0.44; TIP <- 0.018
seg <- t25 %>% mutate(x1 = 1, x2 = match(line, levels(d$line))) %>%
  rowwise() %>%
  reframe(x = c(x1, x1, x2), xend = c(x2, x1, x2),
          y = c(Y, Y, Y) + (x2 - 2) * 0.06,
          yend = c(Y, Y + TIP, Y + TIP) + (x2 - 2) * 0.06)
lb <- t25 %>% mutate(x2 = match(line, levels(d$line)),
                     x = (1 + x2) / 2, y = Y + TIP + 0.008 + (x2 - 2) * 0.06)

pB <- ggplot(d25, aes(line, p, fill = line)) +
  geom_col(colour = "black", width = 0.66, linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.22, colour = "black",
                linewidth = 0.45) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_richtext(data = lb, aes(x = x, y = y, label = fmt_p(fisher_p)),
                inherit.aes = FALSE, size = 2.9, vjust = 0, fill = NA,
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_fill_manual(values = setNames(c(COL_N2, COL_EDIT, COL_EDIT),
                                      levels(d$line)), guide = "none") +
  scale_x_discrete(labels = function(l)
    sprintf("%s\nn = %d", sub("  ", "\n", l),
            d25$n_plated[match(l, as.character(d25$line))])) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.60), expand = expansion(0)) +
  labs(x = NULL, y = "Embryos hatched",
       title = panel_title("B", "**96K makes N2 more sensitive**"),
       subtitle = wrap_md(paste0(
         "The 25% dose alone. Both edited lines drop from 32% to about 4% ",
         "hatching, an eightfold loss of survival, and they agree with each ",
         "other. Same direction as JU2466 (96K, sensitive) against JU1793 ",
         "(96T, resistant), now in a third background."), 62)) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8.6, lineheight = 0.95))

fig <- pA + pB + plot_layout(widths = c(1.25, 1))
ggsave(file.path(OUT, "SUPP_FIG_XX_n2_swap_dose.pdf"), fig,
       width = 11.0, height = 5.2, device = cairo_pdf)
ggsave(file.path(OUT, "SUPP_FIG_XX_n2_swap_dose.png"), fig,
       width = 11.0, height = 5.2, dpi = 300, bg = "white")
msg("wrote SUPP_FIG_XX_n2_swap_dose.{pdf,png}")
