From a 231-strain panel to a single residue
================
Assembled 2026-09-09

-   [Results, as a narrative](#results-as-a-narrative)
-   [Conventions that cross every
    figure](#conventions-that-cross-every-figure)
    -   [RNAi dose is not constant](#rnai-dose-is-not-constant)
    -   [Two thresholds, always both](#two-thresholds-always-both)
-   [I. The assay](#i-the-assay)
    -   [Figure S1 — the simulation the project started
        from](#figure-s1--the-simulation-the-project-started-from)
    -   [Figure S2 — the designed DNA
        mixture](#figure-s2--the-designed-dna-mixture)
    -   [Figure 1 — the validated assay, the phenotype it produces, and
        the
        map](#figure-1--the-validated-assay-the-phenotype-it-produces-and-the-map)
    -   [Figure S3 — the sample-level measurements behind the Figure 1
        slopes](#figure-s3--the-sample-level-measurements-behind-the-figure-1-slopes)
    -   [Figure S4 — the bootstrap propagation
        checks](#figure-s4--the-bootstrap-propagation-checks)
    -   [Figure S5 — the sequencing-depth
        requirement](#figure-s5--the-sequencing-depth-requirement)
    -   [Figure S6 — replicate reproducibility of the
        phenotype](#figure-s6--replicate-reproducibility-of-the-phenotype)
    -   [Figure S7 — the plate assay, validated
        externally](#figure-s7--the-plate-assay-validated-externally)
-   [II. The map](#ii-the-map)
    -   [Figure 2 — pooled GWAS and cross mapping identify overlapping
        loci](#figure-2--pooled-gwas-and-cross-mapping-identify-overlapping-loci)
    -   [Figure S8 — why these cross
        parents](#figure-s8--why-these-cross-parents)
    -   [Figure S9 — the expanded view behind Figure 2’s
        tracks](#figure-s9--the-expanded-view-behind-figure-2s-tracks)
-   [III. The interval](#iii-the-interval)
    -   [Figure 3 — NILs fine-map the chromosome III QTL to 37
        kb](#figure-3--nils-fine-map-the-chromosome-iii-qtl-to-37-kb)
    -   [Figure S10 — the complete NIL hatching
        experiment](#figure-s10--the-complete-nil-hatching-experiment)
-   [IV. The residue](#iv-the-residue)
    -   [Figure 4 — editing *sid-2* residue 96 moves sensitivity in
        three
        backgrounds](#figure-4--editing-sid-2-residue-96-moves-sensitivity-in-three-backgrounds)
    -   [Effect sizes across the three chromosome III
        experiments](#effect-sizes-across-the-three-chromosome-iii-experiments)
    -   [Figure S11 — the full N2 dose
        series](#figure-s11--the-full-n2-dose-series)
    -   [Figure S12 — everything held back from Figure
        4A](#figure-s12--everything-held-back-from-figure-4a)
    -   [Figure S13 — the honest negative
        check](#figure-s13--the-honest-negative-check)
    -   [Figure S14 — surface charge of the
        ectodomain](#figure-s14--surface-charge-of-the-ectodomain)
-   [Diagnostics](#diagnostics)
    -   [Leakage in the MIP-seq
        validation](#leakage-in-the-mip-seq-validation)
    -   [GWAS interval admission](#gwas-interval-admission)
    -   [Coverage against reference
        size](#coverage-against-reference-size)
    -   [Off-pool leakage against
        coverage](#off-pool-leakage-against-coverage)
    -   [Gene content at the four strongest mig-6
        QTL](#gene-content-at-the-four-strongest-mig-6-qtl)
    -   [Figure S19 — SID-2 across two species, with elegans variation
        on
        top](#figure-s19--sid-2-across-two-species-with-elegans-variation-on-top)
    -   [eQTL and parental expression at the censused
        loci](#eqtl-and-parental-expression-at-the-censused-loci)
    -   [Figure S15 — where T96’s pocket sits in the charge
        distribution](#figure-s15--where-t96s-pocket-sits-in-the-charge-distribution)
    -   [Figure S16 — model confidence, and the proximity
        null](#figure-s16--model-confidence-and-the-proximity-null)
    -   [Figure S17 — Figure 2 without the cross
        QTL](#figure-s17--figure-2-without-the-cross-qtl)
    -   [Figure S18 — what the 37 kb interval
        contains](#figure-s18--what-the-37-kb-interval-contains)
    -   [The panel split at two loci](#the-panel-split-at-two-loci)
-   [Open before submission](#open-before-submission)
-   [Figure manifest](#figure-manifest)

<!--
FIGURE_REPORT.Rmd -- the twenty-one manuscript figures with their captions, ordered
by the argument rather than by build order.

  Rscript -e 'rmarkdown::render("FIGURE_REPORT.Rmd", "all")'

builds both targets:

  FIGURE_REPORT.md     the GitHub-rendered version. Figures are referenced from
                       plots/ rather than embedded, so GitHub displays them and
                       the file stays a few tens of kB. This is the tracked one.
  FIGURE_REPORT.html   the styled, self-contained version (~13 MB, figures
                       base64-embedded, click-to-zoom). Git-ignored; rebuild it
                       when you want the designed reading copy.

The CSS, the webfont link and the zoom overlay are guarded on the output
format, so none of them leak into the Markdown.

Run from the repository root; the Rmd lives there so that "plots/..." and
"supplemental_data/..." resolve without any setwd().

WHAT IS DERIVED AND WHAT IS PROSE. Every table below marked "derived" is
computed from supplemental_data/ when this file knits, so it cannot drift from
the deposit. Four such tables were checked against the numbers in
FIGURE_CAPTIONS.txt and reproduce them exactly. Numbers inside the caption prose
came from the figure scripts' console output and are literal text; the figure
scripts remain the authority for those.

The figures themselves are read from plots/ and base64-embedded by pandoc, so
the knitted HTML is a single self-contained file.
-->

The manuscript figure set for RNAi sensitivity in wild *Caenorhabditis
elegans*, ordered by the argument rather than by build order. Each
figure carries its caption, with the caveats kept attached to the panel
they qualify. Click any figure to enlarge it.

The argument narrows by roughly an order of magnitude at each step:
**231 wild isolates** given a quantitative phenotype, **six
chromosomes** scanned in two crosses, **37 kb** resolved by a NIL
series, **one residue** edited in three backgrounds.

# Results, as a narrative

A draft Results section in manuscript voice, written to be a reference
point rather than final text. Every number in it is one this repository
computes, and each is traceable to the figure entry below that produces
it. Where the data will not carry a sentence, the sentence is not made —
the closing box lists the claims to avoid.

Figure numbers follow the curated set: four main figures and fourteen
supplements, `S1`–`S17`, in the order the argument uses them.

<div class="ms">

### A pooled sequencing assay measures RNAi response across wild isolates

To measure the RNAi response of many wild isolates at once, we pooled
strains, exposed the pool to RNAi, and sequenced it, inferring each
strain’s frequency from bulk allele counts by non-negative least squares
(NNLS) regression. Because every downstream result depends on that
inference, we validated it three times, each time conceding something
the previous test controlled.

We first asked whether the inference is possible in principle. We
assigned each wild isolate a fitness value drawn from an inverse χ²
distribution, computed the pooled allele frequencies such a population
would produce, simulated observed alt-allele counts by binomial sampling
across 1–500× coverage, and deconvolved the simulated counts back to
per-strain frequencies <span class="cite">(Figure S1)</span>. Using
seven published traits with validated QTL to seed seven independent
populations of 327 strains, recovery of the known input was essentially
exact at high coverage (r² = 1.00 at 500× for all seven traits) and
degraded gradually as coverage fell. At 1× coverage recovery was
trait-dependent, ranging from r² = 0.91 to 0.52 (median 0.79); 10× was
sufficient for five of the seven traits, and 30× was the lowest coverage
at which all seven reached r² ≥ 0.95.

We then asked whether real libraries behave as the simulation predicts.
We divided 174 wild isolates into four sets, pooled genomic DNA from
each, and sequenced the pools both pure and as a seven-step titration of
one set against another <span class="cite">(Figure S2)</span>. Pure
pools returned 0.775–0.853 of their own set. Across the titration the
two titrated sets traded off monotonically — set B rising from 0.12 to
0.72 and set C falling from 0.70 to 0.10 (Spearman ρ = +1 and −1 against
titration step) — while the two untitrated sets remained flat. Because
only two sets were titrated against each other, their combined share of
the pool must remain constant however the DNA was mixed; it did, to a
standard deviation of 0.86% and a maximum departure of 1.38%. Against
the designed proportions themselves the recovered fractions were
accurate to a root mean squared error of 0.038 (Pearson r = 0.997), with
the largest single deviation at the step whose 0.1 µL of set B was the
smallest volume pipetted; each dilution was prepared once, so pipetting
error is unreplicated and enters that figure in full.

Two limits of the inference emerged from the same experiment. Roughly a
fifth of each pure pool was assigned to strains absent from it, and that
misassignment was predicted by relatedness: a strain’s frequency in
pools it did not belong to rose with its identity-by-state to the
nearest other strain in the reference (Spearman ρ = 0.326, p = 1.4 ×
10⁻⁵, n = 170). The fraction of its own pool a strain recovered, by
contrast, was uncorrelated with relatedness (ρ = 0.046, p = 0.55).
Genetic similarity therefore causes strains to absorb signal from one
another without systematically depleting their own estimates, and the
shortfall in pure-pool recovery has some other cause.

Finally we compared the inference against an independent measurement of
the same material. For the one experiment in which both pooled
whole-genome sequence and published targeted MIP-seq exist — an L1
starvation time course — per-strain rates of change in pool frequency
agreed closely across platforms (Spearman ρ = 0.97, p \< 10⁻⁴, n = 98
strains) <span class="cite">(Figure 1A)</span>. Repeating the entire
slope calculation inside each of 100 bootstrap replicates of the
deconvolution gave intervals narrower than the platform disagreement for
50 of 98 strains and 0.57× it on average <span class="cite">(Figure
S4)</span>, so the residual scatter reflects disagreement between
platforms rather than noise in the inference. Agreement was lower and
more variable at the level of individual samples than of fitted slopes
(median ρ = 0.84, minimum 0.517) <span class="cite">(Figures S3,
S5)</span>, and saturated by 3× coverage, consistent with the
simulation.

Applying the assay to 231 wild isotypes exposed to *pos-1* RNAi produced
a continuously distributed response phenotype on a variance-stabilised
scale <span class="cite">(Figure 1B)</span>, reproducible across
replicate pools <span class="cite">(Figure S6)</span> and correlated in
the expected direction with manual plate scoring of the same strains
(Spearman ρ = 0.41, n = 111, p = 7.8 × 10⁻⁶) and with published
embryonic-lethality measurements (ρ = −0.55, n = 19, p = 0.014) <span
class="cite">(Figure S7)</span>.

### Pooled association and cross mapping converge on overlapping loci

A genome-wide association scan of the *pos-1* response across 231
isotypes and 464,045 markers identified ten markers exceeding a
Bonferroni threshold (maximum −log₁₀*p* = 8.84) and 465 exceeding a
threshold corrected for the 1,972 effective independent tests <span
class="cite">(Figure 1C)</span>. Eight of the ten lay on chromosome IV,
in two clusters near 13.41 and 15.32 Mb, with single markers on
chromosome III at 5.97 Mb (−log₁₀*p* = 8.68) and chromosome X at 4.88
Mb.

To separate loci affecting the RNAi response generally from those
specific to one target, we measured the pooled response to several RNAi
targets across 84 isotypes and mapped the same responses in F2
bulk-segregant crosses <span class="cite">(Figure 2)</span>. Cross
parents were drawn from opposite extremes of the pooled assay rather
than for convenience <span class="cite">(Figure S8)</span>. Nine cross
QTL exceeded LOD 100 across the two crosses, on chromosomes I, II, III,
IV, V and X. Contrasting the response to *pos-1* knockdown against the
response to knockdown of an unrelated target distinguished the two
classes <span class="cite">(Figure S9)</span>: a locus on the right arm
of chromosome III behaved as a general RNAi-response locus, while loci
on chromosomes I, V and X were target-specific.

Association and cross mapping agreed at the level of locus rather than
of marker. No cross interval contained its matching association peak,
the closest correspondences being 0.84 Mb and 0.96 Mb apart, which is
the resolution a panel of 84 phenotyped strains supports. We therefore
pursued the chromosome III locus by introgression rather than by
association.

### Near-isogenic lines resolve the chromosome III locus to 37 kb

To fine-map the general RNAi-response locus, we constructed
near-isogenic lines carrying overlapping introgressions across the right
arm of chromosome III in the two cross-parent backgrounds and scored
embryonic hatching on *pos-1* RNAi <span class="cite">(Figure 3)</span>.
Hatching formed a graded series across the introgression series — 99.4%
in the resistant parent, 97.3%, 79.4% and 58.5% in successive lines, and
35.7% in the sensitive parent — and the smallest interval distinguishing
a resistant from a sensitive line spanned **37 kb**, from 13.658 to
13.695 Mb. Control hatching was 97–100% for every line <span
class="cite">(Figure S10)</span>, so the differences are attributable to
the RNAi exposure rather than to the introgressions themselves.

### A missense variant in *sid-2* accounts for part of the response

The resolved interval contains *sid-2*, which encodes an intestinal
transmembrane protein required for the uptake of ingested
double-stranded RNA. Of eight annotated protein-altering variants
segregating in *sid-2* in the wild population, only two differ between
the cross parents, one of which is a threonine-to-lysine substitution at
residue 96 <span class="cite">(Figure 4D)</span>. A second common
variant, P153T, is in near-complete linkage disequilibrium with T96K
across the wild population (r² = 0.935) but is carried by both cross
parents and therefore cannot contribute to the mapped difference.

Reciprocal editing of residue 96 in both parental backgrounds moved
hatching in both directions <span class="cite">(Figure 4A)</span>.
Introducing 96K into the resistant parent reduced hatching from 94.8% to
53.1% (Fisher’s exact test, p = 1.7 × 10⁻²⁴), and restoring 96T in the
sensitive parent raised hatching from 5.4% to 18.4% (p = 3.9 × 10⁻⁵).
Residue 96 therefore contributes in both directions without accounting
for the whole difference between the parents. The same substitution
introduced into the laboratory reference background reduced hatching
from 32.3% to 4.4% and 3.7% in two independently derived lines (p = 7.8
× 10⁻¹⁷ and 6.2 × 10⁻¹⁹) <span class="cite">(Figure 4B)</span> — a third
genetic background, and the only one with independent edits, with the
direction unchanged throughout: 96K sensitive, 96T resistant. That
comparison required a sub-maximal RNAi dose: in the reference background
25% *pos-1* bacteria is the only dilution with dynamic range, since
every genotype hatches without RNAi and every genotype is inviable from
50% upward <span class="cite">(Figure S11)</span>. Hatching percentages
are therefore not comparable between the two backgrounds.

Two further observations bound the interpretation. First, T96K has no
marginal effect across the mapping panel: the variant reaches p = 0.24
in the association scan, ranking 18,662 of 64,423 markers on chromosome
III, and splitting the pooled phenotype by residue 96 gives no shift
(Wilcoxon p = 0.99, r² = 0.007) <span class="cite">(Figure S13)</span>.
*sid-2* was identified by the crosses and the introgression series, not
by association, and a variant at 36% frequency whose effect is this
context-dependent would not be expected to surface in a marginal test.
Second, residue 96 lies on the lumenal face of the predicted ectodomain,
on the same face as three residues with published effects on dsRNA
uptake <span class="cite">(Figure 4C)</span>, and the substitution
raises the domain’s net charge at gut-lumen pH from −0.2 e to +0.8 e
<span class="cite">(Figure S14)</span>. Both observations are spatial
and electrostatic context; neither is evidence of a shared binding site,
and 41% of the ectodomain lies within 20 Å of residue 96 (binomial p =
0.19, permutation p = 0.30). A predicted N-glycosylation sequon is
removed by T96K, but an edit that removes the same sequon while leaving
residue 96 intact does not phenocopy it <span class="cite">(Figure
S12)</span>, so no glycosylation mechanism is proposed.

</div>

<div class="msnote">

#### Claims to avoid, each because a figure contradicts or fails to support it

-   **“Accurate at 1× coverage.”** True for the best traits, false for
    the worst (r² 0.91 to 0.52). Quote the median, or quote 10× and 30×.
-   **“The dilution series recovers the intended proportions.”** The
    intended proportions are not recorded. Say monotonic and
    complementary, and quote the conservation bound (0.86% sd) as the
    accuracy statement.
-   **“The association scan found nothing.”** It found ten
    Bonferroni-significant markers, eight of them on chromosome IV. What
    it did not find is the chromosome III interval the crosses resolved.
-   **“The association scan identified *sid-2*.”** It did not — the
    variant ranks 18,662 of 64,423, and the chromosome III association
    peak at 5.97 Mb lies 7.7 Mb away. The crosses and the NILs
    identified it.
-   **“The chromosome III peak coincides with the NIL interval.”** The
    peak marker at 13.784 Mb is the chromosome’s terminal marker, so its
    position is a boundary artefact. The 37 kb NIL interval is the
    claim.
-   **Putting Figure 4A and 4B hatching percentages side by side.** 50%
    versus 25% *pos-1* RNAi; the numbers are not comparable.
-   **“T96 sits at the dsRNA binding site.”** Proximity only, and not
    significant (p = 0.19–0.30).
-   **“T96K acts by removing a glycosylation site.”** The N94A control
    fails.
-   **Any statement about what the bootstrap resampled.** Only its
    outputs were archived.

</div>

# Conventions that cross every figure

## RNAi dose is not constant

Every embryo-hatching experiment ran on a lawn of *pos-1* RNAi bacteria
diluted with HT115, and the dilution differs.

| Experiment                    | Dose            |
|:------------------------------|:----------------|
| Figure 3C, and its supplement | 50% *pos-1*     |
| Figure 4A, and its supplement | 50% *pos-1*     |
| Figure 4B                     | **25% *pos-1*** |

Hatching percentages are therefore **not comparable between Figure 4A
and Figure 4B**, and the text should not put those numbers side by side.
The 25% dose is not an oversight: in the N2 background it is the only
dilution with any dynamic range — at 0% every genotype hatches, and from
50% upward every genotype is dead. The dose series below the Figure 4
entry shows this.

## Two thresholds, always both

Bonferroni over every tested marker assumes markers are independent,
which they are not in this panel. The eigen threshold divides α by the
effective number of independent tests from the eigenvalues of the marker
correlation matrix (Li & Ji 2005).

| Panel            |   n | Markers | M<sub>eff</sub> | Bonferroni | Eigen |
|:-----------------|----:|--------:|----------------:|-----------:|------:|
| pooled_RNAi_expt |  84 | 322,010 |             732 |       6.81 |  4.17 |
| pos1_2023        | 231 | 464,045 |           1,972 |       6.97 |  4.60 |

Thresholds as −log<sub>10</sub>*p*. Cross scans use LOD instead, with a
genome-wide threshold of LOD 3.57 (α = 0.05 over 2,000 effective tests).
**vst** throughout is the variance-stabilised pooled phenotype — the
exact trait the GEMMA associations were run on, so phenotype and
Manhattan panels always share a scale.

# I. The assay

<div class="scale">

231 isolates

</div>

<div class="claim">

The method everything else rests on is validated three times, each
conceding something the one before it controlled: in simulation, where
both the input frequencies and the counts are synthetic; on a designed
DNA mixture, where the input is known and the counts are real; and
against published MIP-seq on the same samples, where neither is
controlled. It then yields a quantitative *pos-1* RNAi phenotype across
the wild panel that a genome-wide scan can be run on.

</div>

<div class="aside">

<span class="ch">Read these three in order</span>

Quote the **simulation** for the depth requirement, the **DNA dilution**
for the fact that real libraries behave, and **Figure 1A** for agreement
with an independent published measurement. No one of them carries the
claim alone.

</div>

## Figure S1 — the simulation the project started from

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_simulation_depth.R`<br> **Validates**
NNLS against a known input, with synthetic counts<br> **Scope** 7 traits
× 8 depths (1–500×) × 327 strains<br> **Depth for r² ≥ 0.95** six of
seven by 30×, all seven by 50×

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_simulation_depth.png" alt="Two panels: r-squared against the known input frequency by simulated sequencing depth for seven traits, and estimated frequencies against the known input faceted by depth." width="100%" />
<p class="filecap">
SUPP_FIG_XX_simulation_depth
</p>

</div>

Each wild isolate was assigned a fitness value taken from one of seven
published traits with validated QTL — the trait value itself, not a
draw, shifted by the trait’s own minimum so fitness is non-negative —
giving seven independent populations of 327 strains. A strain with no
published value for a trait was set to fitness zero, so it was absent
from that trait’s pool. The expected pooled allele frequencies such a
population would produce were computed; alt-allele counts were simulated
by binomial sampling at 1, 3, 5, 10, 30, 50, 100 and 500×; and those
counts were deconvolved back to per-strain frequencies by NNLS and
compared with the known input.

<div class="panel">

<span class="pl">A</span> r² of estimated against known input frequency,
against simulated depth, one line per trait — the **mean of ten seeded
replicates, with the band spanning them**, so the band’s width is the
sampling variability of one experiment at that depth. At 1× the spread
across traits is wide — PC1 `0.91`, value `0.88`, amsacrine_f.L1 `0.81`,
assay_norm `0.76`, etoposide_median.TOF `0.75`, Albendazole_q75.TOF
`0.56`, mtDNA_ratio `0.51`; median `0.76` — and the bands are widest
there too, up to `0.149` for assay_norm.

**Six of seven traits reach 0.95 by 30×, and 50× clears all seven.**
Counting a trait only where it clears the bar in *every* replicate: 0
traits at 1×, 1 at 3×, 2 at 5×, **5 at 10×**, **6 at 30×**, **7 at 50×**
and above. The exception at 30× is mtDNA_ratio, which sits *on* the line
— mean `0.9504`, span `0.9443`–`0.9610`, clearing 0.95 in four of the
ten. At 50× the lowest r² anywhere in the ten replicates is `0.9639`.
**No trait reaches r² of 1.000 at 500×** — the range there is `0.9966`
(mtDNA_ratio) to `0.9998` (PC1).

</div>

<div class="panel">

<span class="pl">B</span> The estimates against the **known input**,
faceted by depth, all seven traits pooled — eight facets, because 500×
is an ordinary depth rather than the reference it used to be. **One
replicate is drawn, not an average of the ten**: averaging point clouds
would narrow the scatter and misrepresent a single experiment. Axes
share limits so the dashed <span class="m">y = x</span> line means the
same thing in every facet. Pooled r²: `0.80` at 1×, `0.91` at 3×, `0.95`
at 5×, `0.97` at 10×, `0.99` at 30× and 50×, `1.00` at 100× and 500×.
The band along <span class="m">y = 0</span> is the strains with no
published value for that trait, absent from the pool; NNLS assigns them
frequency anyway, and since pool mass is conserved the strains that were
present are underestimated by the same amount — which is why the cloud
sits below <span class="m">y = x</span> at low depth. Leaked mass runs
`7.3%` at 1× down to `0.3%` at 500×, and the slope of input on estimate
rises from `0.883` to `1.004`.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Trait
</th>
<th style="text-align:right;">
1×
</th>
<th style="text-align:right;">
3×
</th>
<th style="text-align:right;">
5×
</th>
<th style="text-align:right;">
10×
</th>
<th style="text-align:right;">
30×
</th>
<th style="text-align:right;">
50×
</th>
<th style="text-align:right;">
100×
</th>
<th style="text-align:right;">
500×
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PC1
</td>
<td style="text-align:right;">
0.912<br><span style="font-size:80%;color:#5A6664">0.880–0.932</span>
</td>
<td style="text-align:right;">
0.968<br><span style="font-size:80%;color:#5A6664">0.954–0.979</span>
</td>
<td style="text-align:right;">
0.981<br><span style="font-size:80%;color:#5A6664">0.973–0.986</span>
</td>
<td style="text-align:right;">
0.991<br><span style="font-size:80%;color:#5A6664">0.988–0.995</span>
</td>
<td style="text-align:right;">
0.997<br><span style="font-size:80%;color:#5A6664">0.996–0.998</span>
</td>
<td style="text-align:right;">
0.998<br><span style="font-size:80%;color:#5A6664">0.997–0.999</span>
</td>
<td style="text-align:right;">
0.999<br><span style="font-size:80%;color:#5A6664">0.999–0.999</span>
</td>
<td style="text-align:right;">
1.000<br><span style="font-size:80%;color:#5A6664">1.000–1.000</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
value
</td>
<td style="text-align:right;">
0.885<br><span style="font-size:80%;color:#5A6664">0.852–0.910</span>
</td>
<td style="text-align:right;">
0.956<br><span style="font-size:80%;color:#5A6664">0.949–0.964</span>
</td>
<td style="text-align:right;">
0.971<br><span style="font-size:80%;color:#5A6664">0.962–0.981</span>
</td>
<td style="text-align:right;">
0.986<br><span style="font-size:80%;color:#5A6664">0.983–0.990</span>
</td>
<td style="text-align:right;">
0.996<br><span style="font-size:80%;color:#5A6664">0.995–0.997</span>
</td>
<td style="text-align:right;">
0.997<br><span style="font-size:80%;color:#5A6664">0.997–0.998</span>
</td>
<td style="text-align:right;">
0.999<br><span style="font-size:80%;color:#5A6664">0.998–0.999</span>
</td>
<td style="text-align:right;">
1.000<br><span style="font-size:80%;color:#5A6664">1.000–1.000</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
amsacrine_f.L1
</td>
<td style="text-align:right;">
0.811<br><span style="font-size:80%;color:#5A6664">0.778–0.853</span>
</td>
<td style="text-align:right;">
0.930<br><span style="font-size:80%;color:#5A6664">0.913–0.949</span>
</td>
<td style="text-align:right;">
0.952<br><span style="font-size:80%;color:#5A6664">0.948–0.956</span>
</td>
<td style="text-align:right;">
0.975<br><span style="font-size:80%;color:#5A6664">0.971–0.981</span>
</td>
<td style="text-align:right;">
0.992<br><span style="font-size:80%;color:#5A6664">0.991–0.993</span>
</td>
<td style="text-align:right;">
0.995<br><span style="font-size:80%;color:#5A6664">0.995–0.996</span>
</td>
<td style="text-align:right;">
0.998<br><span style="font-size:80%;color:#5A6664">0.997–0.998</span>
</td>
<td style="text-align:right;">
1.000<br><span style="font-size:80%;color:#5A6664">0.999–1.000</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
assay_norm
</td>
<td style="text-align:right;">
0.758<br><span style="font-size:80%;color:#5A6664">0.664–0.812</span>
</td>
<td style="text-align:right;">
0.894<br><span style="font-size:80%;color:#5A6664">0.879–0.917</span>
</td>
<td style="text-align:right;">
0.937<br><span style="font-size:80%;color:#5A6664">0.922–0.943</span>
</td>
<td style="text-align:right;">
0.964<br><span style="font-size:80%;color:#5A6664">0.952–0.978</span>
</td>
<td style="text-align:right;">
0.990<br><span style="font-size:80%;color:#5A6664">0.987–0.992</span>
</td>
<td style="text-align:right;">
0.994<br><span style="font-size:80%;color:#5A6664">0.993–0.994</span>
</td>
<td style="text-align:right;">
0.997<br><span style="font-size:80%;color:#5A6664">0.996–0.998</span>
</td>
<td style="text-align:right;">
0.999<br><span style="font-size:80%;color:#5A6664">0.999–0.999</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
etoposide_median.TOF
</td>
<td style="text-align:right;">
0.753<br><span style="font-size:80%;color:#5A6664">0.722–0.788</span>
</td>
<td style="text-align:right;">
0.902<br><span style="font-size:80%;color:#5A6664">0.867–0.928</span>
</td>
<td style="text-align:right;">
0.937<br><span style="font-size:80%;color:#5A6664">0.918–0.951</span>
</td>
<td style="text-align:right;">
0.968<br><span style="font-size:80%;color:#5A6664">0.958–0.975</span>
</td>
<td style="text-align:right;">
0.991<br><span style="font-size:80%;color:#5A6664">0.988–0.993</span>
</td>
<td style="text-align:right;">
0.994<br><span style="font-size:80%;color:#5A6664">0.993–0.996</span>
</td>
<td style="text-align:right;">
0.997<br><span style="font-size:80%;color:#5A6664">0.995–0.998</span>
</td>
<td style="text-align:right;">
0.999<br><span style="font-size:80%;color:#5A6664">0.999–1.000</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
Albendazole_q75.TOF
</td>
<td style="text-align:right;">
0.556<br><span style="font-size:80%;color:#5A6664">0.504–0.595</span>
</td>
<td style="text-align:right;">
0.773<br><span style="font-size:80%;color:#5A6664">0.717–0.827</span>
</td>
<td style="text-align:right;">
0.847<br><span style="font-size:80%;color:#5A6664">0.829–0.875</span>
</td>
<td style="text-align:right;">
0.919<br><span style="font-size:80%;color:#5A6664">0.905–0.932</span>
</td>
<td style="text-align:right;">
0.971<br><span style="font-size:80%;color:#5A6664">0.963–0.977</span>
</td>
<td style="text-align:right;">
0.983<br><span style="font-size:80%;color:#5A6664">0.979–0.987</span>
</td>
<td style="text-align:right;">
0.991<br><span style="font-size:80%;color:#5A6664">0.990–0.993</span>
</td>
<td style="text-align:right;">
0.998<br><span style="font-size:80%;color:#5A6664">0.998–0.999</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
mtDNA_ratio
</td>
<td style="text-align:right;">
0.510<br><span style="font-size:80%;color:#5A6664">0.458–0.546</span>
</td>
<td style="text-align:right;">
0.719<br><span style="font-size:80%;color:#5A6664">0.675–0.761</span>
</td>
<td style="text-align:right;">
0.770<br><span style="font-size:80%;color:#5A6664">0.732–0.815</span>
</td>
<td style="text-align:right;">
0.867<br><span style="font-size:80%;color:#5A6664">0.837–0.903</span>
</td>
<td style="text-align:right;">
0.950<br><span style="font-size:80%;color:#5A6664">0.944–0.961</span>
</td>
<td style="text-align:right;">
0.967<br><span style="font-size:80%;color:#5A6664">0.964–0.973</span>
</td>
<td style="text-align:right;">
0.982<br><span style="font-size:80%;color:#5A6664">0.980–0.985</span>
</td>
<td style="text-align:right;">
0.997<br><span style="font-size:80%;color:#5A6664">0.996–0.997</span>
</td>
</tr>
</tbody>
</table>

<div class="derived">

Mean over ten seeded replicates, with the replicate range beneath each
value, from supplemental_data/deconvolution/simulation_seeded_r2.tsv
(scripts/make_simulation_seeded.R, seed 20260909).

The 2021 run is kept as a comparison, not discarded: its estimates are
still deposited, all 56 of its r² recompute exactly from them
(scripts/simulation_recompute_r2.R, run on every push), and its values
fall inside the ten-replicate band in **43 of 56** trait-depth cells —
consistent with the replicated run rather than anomalous.

</div>

<div class="derived">

<span class="ch">Provenance — this used to be the figure’s governing
caveat</span>

**The figure is generated from source rather than recovered.** The
simulation’s fitness input was thought lost, so panel A was transcribed
— read out of text embedded in the original per-trait PDFs by
`scripts/extract_sim_reported_r2.py` — and panel B substituted the 500×
estimate for the truth. Neither is necessary: the seven-trait arm of
`scripts/legacy/haploReg_original.R` uses published trait values as
fitness, those files are deposited as
`supplemental_data/deconvolution/simulation_fitness_traits.tsv`, and the
whole simulation now runs from the genotype panel **with a seed** — ten
replicates, `scripts/make_simulation_seeded.R`. That closes the
reproducibility gap the missing seed left, and the seeded run carries
**no negative coefficients** against the archive’s 139 of 18,312.

One caveat survives the change: the r² are computed over all 327 strains
with the unmeasured ones held at exactly zero, which flatters the
correlation. The strains actually carrying a published value run from
327 (mtDNA_ratio) down to 84 (PC1), and restricted to those the pooled
r² falls from `0.79` to `0.70` at 1× and from `0.97` to `0.95` at 10×.

</div>

<div class="caveat">

<span class="ch">On the wording of the claim</span>

“NNLS can accurately infer strain frequencies with as little as 1×
sequencing depth” holds for the best-behaved traits and not for the
worst: at 1×, r² runs from **0.91 down to 0.52**. What the whole figure
supports is that 1× recovers most of the signal for most traits (median
r² 0.79), that 10× suffices for five traits of seven, and that **50×**
is the lowest depth at which all seven sit at or above 0.95. Earlier
drafts said 30×, which came from applying the 0.95 bar to the
transcribed two-decimal values: mtDNA_ratio at 30× is `0.9477`, which
displays as 0.95 but is below the bar.

Separately, 139 of 18,312 archived coefficients (0.8%) are negative, at
depths 5 through 100. A strict non-negative solver cannot return those,
so the archive carries solver or post-processing noise of order 1e-2 on
the count scale. It changes no conclusion, but the archived values are
not exactly a clean NNLS output.

</div>

## Figure S2 — the designed DNA mixture

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_dilution_validation.R`<br> **Validates**
NNLS against a known input, with real counts<br> **Design** 174 isolates
in 4 sets; pure pools in triplicate + a 7-step B-into-C titration<br>
**Counts** GATK ASEReadCounter

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_dilution_validation.png" alt="Three panels: pure-pool recovery per set, the B-into-C titration on the pool-wide reference with values labelled, and the same titration renormalised within sets B and C." width="100%" />
<p class="filecap">
SUPP_FIG_XX_dilution_validation
</p>

</div>

<div class="panel">

<span class="pl">A</span> Pure pools, three libraries each. The fraction
assigned to the pool’s own set is `0.785` (A), `0.775` (B), `0.796` (C),
`0.853` (D) — filled points. The fraction misassigned to each other set
runs `0.020`–`0.127`, mean `0.066` (open points); total misassigned per
pool `0.135`–`0.242`.

</div>

<div class="panel">

<span class="pl">B</span> The titration on the 170-strain pool
reference, each B and C value labelled. Set B rises monotonically 0.12 →
0.72 and set C falls monotonically 0.70 → 0.10 (Spearman ρ `+1` and `−1`
against step). The two untitrated sets stay flat — A `0.061`–`0.086`, D
`0.093`–`0.116`, drawn dotted.

</div>

<div class="panel">

<span class="pl">C</span> The same seven samples with the reference
restricted to the 84 strains of sets B and C and renormalised — the
analysis the original figure showed: B 0.17 → 0.84 against C 0.83 →
0.16, crossing between BC4 and BC5.

</div>

<div class="panel">

<span class="pl">D</span> Whether the strains that resolve badly are the
genetically similar ones. For each of the 170 strains, x is
identity-by-state to the closest *other* strain in the reference and y
is its **leakage** — the mean frequency it picks up across the nine
libraries of the three sets it is **not** in. That must be zero whatever
the DNA input was, so unlike own-pool recovery it assumes nothing about
the design. Grey squares are bin medians; points are coloured by the
strain’s own set.

</div>

**Leakage rises with genetic similarity: Spearman ρ `+0.326`,
`p = 1.4e-05`, n = 170.** Mean identity to all strains gives the same
answer more weakly (ρ `+0.258`, `p = 6.8e-04`).

<table>
<thead>
<tr>
<th style="text-align:left;">
Bin
</th>
<th style="text-align:right;">
Strains
</th>
<th style="text-align:right;">
Median leakage
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
\[0.66,0.9\]
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.00000
</td>
</tr>
<tr>
<td style="text-align:left;">
(0.9,0.95\]
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
0.00010
</td>
</tr>
<tr>
<td style="text-align:left;">
(0.95,0.97\]
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
0.00023
</td>
</tr>
<tr>
<td style="text-align:left;">
(0.97,0.99\]
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
0.00143
</td>
</tr>
<tr>
<td style="text-align:left;">
(0.99,1\]
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.00106
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from dilution_predictions_poolref.tsv.gz +
dilution_strain_similarity.tsv

</div>

<div class="aside">

<span class="ch">The ceiling and the leakage are one quantity, not
two</span>

For a pure pool the fractions sum to 1, so **own-set share = 1 − (mass
assigned to strains not in the pool)**. Checked: `own + out = 1` to
4.4e-16 across all twelve pure pools, mean own-set share `0.802` against
out-of-pool `0.198`. The \~0.8 ceiling in panel A *is* the leakage in
panel D, seen from the other side.

**And the leakage is relatedness-directed.** Across the 1,530
strain-by-pool combinations where a strain is absent from the pool, mass
assigned to it averages `1.33` per mille when its nearest neighbour is
elsewhere and `2.25` per mille when its nearest neighbour is *in* that
pool — a 1.70× enrichment, Mann-Whitney `p = 3.6e-05`. Those pairs are
24.3% of the combinations but carry **35.3%** of all out-of-pool mass.
Rescaling, relatedness accounts for about **14%** of the out-of-pool
mass, so it is a real contributor to the ceiling.

That 14% is a **lower bound**, because “nearest neighbour in the pool”
is the crudest possible measure: a strain has many relatives among a
pool’s \~43 members, not one. The proper test is mass against
relatedness to *all* pool members, which needs the full genotype matrix.

What the per-strain null actually showed.
`ρ(nn_ibs, own-pool recovery) = +0.046, p = 0.55` is about *per-strain*
recovery, which is **bidirectional** — a strain both loses mass to
confusable partners and gains it from them, and a trade with a partner
in the same set does not move the set total at all. So per-strain
recovery can show no correlation while relatedness drives every
transfer. Leakage is unidirectional, which is why the signal survives
there.

</div>

<div class="aside">

<span class="ch">Restricting the reference does recover the input ratios
better</span>

Against the designed B fraction, in RMSE:

| estimate                                      |     RMSE |      bias |
|:----------------------------------------------|---------:|----------:|
| 170-strain reference, raw                     | `0.0788` | `−0.0593` |
| 170-strain reference, renormalised within B+C | `0.0415` | `+0.0297` |
| 84-strain B+C reference                       | `0.0376` | `+0.0208` |

Most of the damage is undone by renormalising — that removes the \~18%
of mass sitting on sets A and D — and restricting the reference to the
strains actually present recovers a further **9%**. So the panel B
against panel C comparison reads the way you would expect: the fewer
absent candidates the solver is offered, the closer the set frequencies
land to the input.

One thing these data still **cannot** support: that the damage
concentrates in near-identical pairs *split across sets*. Only 30
reciprocal nearest-neighbour pairs exist, two above IBS 0.99, so that
comparison is underpowered and answers differently at different
thresholds.

</div>

Both panels’ plotted fractions are pinned as literals in the script and
checked against the deposit on every run, so these numbers cannot drift
from the figure without the script failing.

### The renormalisation, step by step

There are two normalisations here and they do different things. The
second one is why panel C looks the way it does, and it is worth being
explicit about.

<div class="panel">

<span class="pl">1</span> **Within each sample, after the solve.** NNLS
fits the observed alt-allele counts of one sample as a non-negative
combination of the reference strains’ genotypes, returning one
coefficient per strain. Those coefficients are on an arbitrary scale —
they grow with sequencing depth — so each sample’s vector is divided by
its own sum. After this every sample’s strain frequencies sum to exactly
1 and read as “the fraction of this pool contributed by that strain”.
The division is *within* a sample, so it cannot move signal between
samples.

</div>

<div class="panel">

<span class="pl">2</span> **Which strains the sum runs over.** Step 1’s
denominator covers whatever is in the reference, so the reference
decides what “the whole pool” means. On the **170-strain pool
reference** (panels A, B) the denominator includes all four sets, so B +
C comes to about `0.82` and the remaining `0.18` is assigned to the
untitrated sets A and D. On the **84-strain B+C reference** (panel C)
the reference contains only sets B and C, so the denominator runs over
those strains alone and **B + C = 1 by construction**.

</div>

<div class="panel">

<span class="pl">3</span> **Set-level sums.** Per-strain frequencies are
added within a set. Nothing is renormalised at this stage.

</div>

<div class="aside">

<span class="ch">What this means for reading panel C</span>

**Panel C’s two series are exact complements because the normalisation
makes them so, not because the inference discovered it.** B and C
summing to 1 there carries no information. The informative quantity is
in panel B, where the denominator is the whole panel and B + C is free
to move — and it barely does, which is metric 1 below.

Before either normalisation, markers were restricted to mean depth above
5 and below 100 across samples, high-frequency variants were recoded so
the minor allele is counted throughout, and sites with any missing
genotype were dropped.

</div>

### How far off was the inference?

<table>
<thead>
<tr>
<th style="text-align:left;">
Sample
</th>
<th style="text-align:right;">
Set B
</th>
<th style="text-align:right;">
Set C
</th>
<th style="text-align:right;">
B + C
</th>
<th style="text-align:right;">
A + D
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BC1
</td>
<td style="text-align:right;">
0.121
</td>
<td style="text-align:right;">
0.701
</td>
<td style="text-align:right;">
0.822
</td>
<td style="text-align:right;">
0.178
</td>
</tr>
<tr>
<td style="text-align:left;">
BC2
</td>
<td style="text-align:right;">
0.157
</td>
<td style="text-align:right;">
0.669
</td>
<td style="text-align:right;">
0.826
</td>
<td style="text-align:right;">
0.174
</td>
</tr>
<tr>
<td style="text-align:left;">
BC3
</td>
<td style="text-align:right;">
0.211
</td>
<td style="text-align:right;">
0.614
</td>
<td style="text-align:right;">
0.825
</td>
<td style="text-align:right;">
0.175
</td>
</tr>
<tr>
<td style="text-align:left;">
BC4
</td>
<td style="text-align:right;">
0.395
</td>
<td style="text-align:right;">
0.413
</td>
<td style="text-align:right;">
0.808
</td>
<td style="text-align:right;">
0.192
</td>
</tr>
<tr>
<td style="text-align:left;">
BC5
</td>
<td style="text-align:right;">
0.530
</td>
<td style="text-align:right;">
0.281
</td>
<td style="text-align:right;">
0.811
</td>
<td style="text-align:right;">
0.189
</td>
</tr>
<tr>
<td style="text-align:left;">
BC6
</td>
<td style="text-align:right;">
0.678
</td>
<td style="text-align:right;">
0.144
</td>
<td style="text-align:right;">
0.822
</td>
<td style="text-align:right;">
0.178
</td>
</tr>
<tr>
<td style="text-align:left;">
BC7
</td>
<td style="text-align:right;">
0.724
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
0.824
</td>
<td style="text-align:right;">
0.176
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from
supplemental_data/deconvolution/dilution_predictions_poolref.tsv.gz

</div>

**Metric 1 — conservation of the titrated pair. Assumes nothing about
the design, and is the one to quote.** Only B and C are titrated against
each other, so whatever the intended proportions were, the pair’s
combined share of the pool has to stay constant across BC1–BC7: moving
DNA from C to B cannot change how much of the pool is B-or-C. Any wobble
in `B + C` is therefore inference error, measurable without knowing a
single nominal ratio. Observed: mean 0.8198, sd 0.0070 (0.86% of the
mean), largest deviation 0.0113 (1.38%). **The pair is conserved to
better than 1.5% across a titration that moves each set sixfold**, which
bounds the inference error without invoking the design at all.

**Metric 2 — deviation from the designed series.** The design was
recovered from the lab record on 2026-09-08 and is in
`supplemental_data/deconvolution/dilution_design.tsv`: a two-fold
doubling series of set B against a fixed 1 µL of set C, made up to 10 µL
with water. Stocks measured 100 ng/µL (B1) against 99.9 ng/µL (C1), so
the DNA mass fraction equals the volume fraction to within `2.5e-4` —
the equal-concentration assumption is verified, not assumed, and the
correction is four orders of magnitude below the error being measured.

Recovery tracks the design at **Pearson r = 0.997**, **RMSE 0.038** in
fraction units, Spearman ρ `+1`:

<table>
<thead>
<tr>
<th style="text-align:left;">
Sample
</th>
<th style="text-align:right;">
B (µL)
</th>
<th style="text-align:right;">
Total DNA (ng)
</th>
<th style="text-align:right;">
Designed
</th>
<th style="text-align:right;">
Recovered
</th>
<th style="text-align:right;">
Deviation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BC1
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
+0.080
</td>
</tr>
<tr>
<td style="text-align:left;">
BC2
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.167
</td>
<td style="text-align:right;">
0.209
</td>
<td style="text-align:right;">
+0.043
</td>
</tr>
<tr>
<td style="text-align:left;">
BC3
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.286
</td>
<td style="text-align:right;">
0.285
</td>
<td style="text-align:right;">
-0.001
</td>
</tr>
<tr>
<td style="text-align:left;">
BC4
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.445
</td>
<td style="text-align:right;">
0.478
</td>
<td style="text-align:right;">
+0.033
</td>
</tr>
<tr>
<td style="text-align:left;">
BC5
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
0.616
</td>
<td style="text-align:right;">
0.620
</td>
<td style="text-align:right;">
+0.004
</td>
</tr>
<tr>
<td style="text-align:left;">
BC6
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
0.762
</td>
<td style="text-align:right;">
0.771
</td>
<td style="text-align:right;">
+0.009
</td>
</tr>
<tr>
<td style="text-align:left;">
BC7
</td>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
0.865
</td>
<td style="text-align:right;">
0.843
</td>
<td style="text-align:right;">
-0.022
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from dilution_design.tsv + dilution_predictions_bcref.tsv.gz

</div>

**Metric 3 — shape.** The recovered B share is linear in titration step
with R² `0.974`, and Spearman ρ `+1` against both step and designed
fraction.

<div class="panel">

<span class="pl">E</span> Recovered against designed fraction of set B,
point area giving the volume of B pipetted. This is the accuracy panel,
and it exists because the design was recovered.

</div>

<div class="caveat">

<span class="ch">The error is unreplicated, and one step dominates
it</span>

**There are no replicate dilutions**, so pipetting error enters the
comparison in full. The largest deviation is at **BC1** (`+0.080`),
whose B volume is 0.1 µL — the smallest and hardest-to-pipette volume in
the series. Dropping BC1 halves the error to RMSE `0.024`, max `0.043`.
Deviation does trend with 1/volume across the series, but not
significantly at n = 7 (Spearman `+0.50`, `p = 0.27`), so the pipetting
explanation is mechanistically plausible and statistically unproven.

**Total DNA is not constant** across the series either — 11 ng at BC1
rising to 74 ng at BC7, because only the B volume varied. Input mass and
B fraction are therefore perfectly confounded, and a deviation that
scaled with total DNA would be indistinguishable from one that scaled
with B fraction.

**The bias runs toward B**: mean deviation `+0.021` on the B+C
reference, `+0.030` on the pool reference. Stock concentration is
excluded as the cause (100 vs 99.9 ng/µL). Set B holds 46 of the 84
reference columns against set C’s 38, which is a candidate explanation
and is not tested here.

</div>

<div class="caveat">

<span class="ch">One isotype is assigned to two sets, and it is not
cosmetic</span>

Strains JU1580 and JU1793 are the same isotype, and the experiment put
them in sets B and D. The deconvolution returns one frequency for the
pair, so the source tables carry that isotype **twice** — once labelled
B, once D. Summed naively the single estimate is added to both sets,
where it is **5–19% of set D’s apparent frequency** in BC1–BC7 and
*rises across the titration*, making the untitrated set D appear to
drift upward. The figure resolves it to set B — what the original
analysis did, and what the estimate’s own behaviour says — and asserts
no other isotype is ambiguous. This is the same strain whose duplicated
genotype-reference entry is documented for the *pos-1* phenotype, and it
is a cross parent in Figures 2 and 3.

</div>

<div class="caveat">

<span class="ch">Two further limits</span>

Pure-pool recovery tops out near 0.8, not 1.0, so roughly a fifth of
each pool is assigned to strains that are not in it — and that fifth
*is* the leakage of panel D, since the two sum to 1 by construction.
Relatedness demonstrably directs it (1.70× enrichment when the absent
strain’s nearest neighbour is in the pool, `p = 3.6e-05`) and accounts
for at least \~14% of it. The rest is not attributed here:
cross-contamination between pools and unequal DNA input across strains
within a pool remain candidates these data cannot separate.

Panels A and B use the 170-strain pool reference and panel C the
84-strain B+C reference. After renormalising the pool reference within B
and C the two disagree by at most `0.054`, mean `0.029`, on the same
samples — the size of the reference-choice effect. The deposit also
carries the 540-strain full-panel and regenotyped references.

</div>

## Figure 1 — the validated assay, the phenotype it produces, and the map

<div class="meta">

**Script** `scripts/Figure1_pos1.R`<br> **Test bed** Baugh L1
starvation, the one dataset with pooled WGS and published MIP-seq on the
same samples<br> **n** 98 strains (A) · 231 strains (B, C) · 464,045
markers

</div>

<div class="plate">

<img src="plots/Figure1_pos1.png" alt="Three panels: MIP-seq against NNLS slopes, the pooled pos-1 vst phenotype distribution, and a genome-wide GEMMA association scan." width="100%" />
<p class="filecap">
plots/Figure1_pos1.pdf · .png
</p>

</div>

<div class="panel">

<span class="pl">A</span> Per-strain rate of change in pool frequency
during L1 starvation, MIP-seq against NNLS, one point per wild isolate
averaged over five replicate arms — `n = 98`, Spearman `ρ = 0.97`,
`p < 1e-4`. Slopes are the regression of frequency change on day, taking
day 1 as baseline and excluding day 17. The dashed line is `y = x`, not
a fit.

</div>

<div class="panel">

<span class="pl">B</span> Distribution of the 2023 pooled *pos-1* RNAi
response across wild isotypes on the vst scale (`n = 231` strains with a
vst value, of 366 in the trait file). The dashed line marks zero;
positive values are strains that gained pool frequency under *pos-1*
RNAi, i.e. resistant.

</div>

<div class="panel">

<span class="pl">C</span> GEMMA LOCO association scan for that
phenotype: 464,045 markers, `n = 231`. The grey dashed line is
Bonferroni over every marker (`6.97`); the green dotted line divides α
by the 1,972 effective independent tests (`4.60`). Markers clearing
Bonferroni are red, markers clearing only the eigen threshold are green.
**Ten markers clear Bonferroni and 465 clear the eigen threshold.** The
maximum is `8.84` on chromosome IV at 15.323 Mb; eight of the ten lie on
chromosome IV, in clusters near 13.41 and 15.32 Mb, with single markers
on chromosome III at 5.966 Mb (`8.68`) and chromosome X at 4.876 Mb
(`7.83`). **None is near the chromosome III interval the NILs resolve**
— the chromosome III association peak is 7.7 Mb from it.

</div>

<div class="caveat">

<span class="ch">Two things to state in the text</span>

N2 is excluded from panel A throughout: as the reference strain its
genotype shares sites with every other strain in the pool, so its NNLS
frequency is not identified on the same footing as the wild isolates’.

The shipped vst value for JU1793 still carries the duplicate-entry
averaging documented in `scripts/2023_pos1_analysis.R` — which matters,
because JU1793 is a cross parent in Figures 2 and 3.

</div>

## Figure S3 — the sample-level measurements behind the Figure 1 slopes

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_baugh_per_sample_frequencies.R`<br>
**Supports** Figure 1A · 15 samples used by the slope fits

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_baugh_per_sample_frequencies.png" alt="Fifteen faceted scatter plots, one per sample, of MIP-seq against NNLS strain frequency with bootstrap intervals." width="100%" />
<p class="filecap">
SUPP_FIG_XX_baugh_per_sample_frequencies
</p>

</div>

One facet per sample for the 15 samples the slope fits use
(non-baseline, day 17 excluded), one point per wild isolate, published
MIP-seq frequency against NNLS frequency, with 95% percentile intervals
over 100 bootstrap replicates of the deconvolution taken at the sample
level. Rows are one replicate arm’s time course; Spearman ρ is given per
facet. Dashed line is `y = x`, not a fit; N2 excluded.

Vertical bars only: the published MIP-seq frequencies are point values
with no distributed uncertainty, so an interval on that axis would be
invented.

Two things Figure 1 cannot show. First, `rep5_d13` agrees far worse than
any other sample (`ρ = 0.50` against 0.77–0.88 elsewhere) and is the
single low bar in the per-sample distribution. Second, 221 of 1,470
points (15%) are exact NNLS zeros — the non-negativity constraint
clamping a strain out of a pool — while the MIP-seq table reports no
zeros at all. A block of ties at zero cannot be ranked, which is the
likely reason per-sample ρ (median 0.84) is so much lower than the
slope-level 0.97.

## Figure S4 — the bootstrap propagation checks

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_bootstrap_propagation_checks.R`<br>
**Supports** the bootstrap intervals · 100 replicates

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_bootstrap_propagation_checks.png" alt="Six-panel audit of the bootstrap intervals: frequency-level and slope-level agreement, standardised bias, coverage, the correlation's bootstrap distribution, and interval width against deviation." width="100%" />
<p class="filecap">
SUPP_FIG_XX_bootstrap_propagation_checks
</p>

</div>

The intervals are obtained by redoing the whole slope calculation inside
each of the 100 bootstrap replicates of the deconvolution, not by
propagating a stored standard error through that chain. Each panel is a
check that could have failed.

<div class="panel">

<span class="pl">A</span> Frequency level: bootstrap mean against point
estimate, *r* = `0.999984`, `n = 2,346`. Confirms the array resamples
the estimator Figure 1 plots.

</div>

<div class="panel">

<span class="pl">B</span> The same after the transformation, *r* =
`0.999994`, `n = 101`. This is the panel that would break if the day-1
baseline were taken from the point estimate rather than from within each
replicate, or if the replicate arms were averaged in the wrong order.

</div>

<div class="panel">

<span class="pl">C</span> Standardised bias, (point estimate − bootstrap
mean) / bootstrap SD: median −0.029, largest 0.61 SD. Dotted lines at
±1.96.

</div>

<div class="panel">

<span class="pl">D</span> Coverage: each strain’s interval with its
point estimate marked. All 101 of 101 intervals contain their own point
estimate.

</div>

<div class="panel">

<span class="pl">E</span> The correlation’s own bootstrap distribution:
median 0.972, 95% percentile interval `[0.969, 0.975]`, full-data value
0.974.

</div>

<div class="panel">

<span class="pl">F</span> Bootstrap interval width against distance from
`y = x`, per strain. The uncertainty is smaller than the platform
disagreement for 50 of 98 strains, and 0.57× it on average — which is
why the bars are invisible at full scale.

</div>

<div class="caveat">

<span class="ch">Caveat, and it is not resolvable from the archived
data</span>

Only the bootstrap outputs were saved. The archive holds the 102 × 23 ×
100 array and a reshaped copy of it, not the resampling scheme, and the
generating code is not in this repository. These panels establish that
the array resamples the plotted estimator and that the transformation is
carried through it consistently. They cannot establish *what* was
resampled — reads, markers or strains. No statement stronger than “100
bootstrap replicates of the deconvolution” is supported by what is on
disk.

</div>

## Figure S5 — the sequencing-depth requirement

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_downsample_per_sample.R`<br>
**Supports** Figure 1C · depths 0.25× to 10×

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_downsample_per_sample.png" alt="Two panels of Spearman correlation against sequencing depth, per sample and at the slope level, on a shared y axis." width="100%" />
<p class="filecap">
SUPP_FIG_XX_downsample_per_sample
</p>

</div>

<div class="panel">

<span class="pl">A</span> Per-sample agreement against depth: one grey
line per sample (23 samples), dark line the median across samples,
dashed line the median full-depth agreement (0.835). Median ρ rises
0.724 at 0.25× to 0.827 at 5×, saturating by 3–5×.

</div>

<div class="panel">

<span class="pl">B</span> The aggregate, slope-level curve — Figure 1C —
with the 0.5× depth restored: 0.833, 0.870, 0.851, 0.903, 0.906, 0.904
for 0.25× through 10×.

</div>

Panel B correlates slopes, which average fifteen measurements, so it
sits above A at every depth and saturates sooner: averaging removes most
of what depth costs an individual sample. **Quoting the depth
requirement from B alone overstates how well any single sample is
measured.**

Two further points. The spread between samples exceeds the effect of
depth — the median moves about 0.10 across a fortyfold depth range while
samples span roughly 0.43 to 0.88 at any fixed depth. And `rep5_d13` is
worst at every depth (0.428–0.478) without improving, so its
disagreement is not a depth problem and no amount of sequencing would
have fixed it.

The 0.5× inversion visible in B — 0.870, above the 0.851 of the 1× that
follows — is why that depth is omitted from Figure 1C. No such inversion
appears per sample, so it looks like noise in one downsampling draw
rather than a property of that depth. No error bars here: the bootstrap
array holds replicates of the full-depth deconvolution only, so there is
nothing to resample at a reduced depth.

## Figure S6 — replicate reproducibility of the phenotype

<div class="meta">

**Script**
`scripts/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.R`<br>
**Supports** the Figure 1B phenotype · all 6 pairwise comparisons of 4
replicates

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_original_pos1_dfreq_rep_correlation.png" alt="Six faceted scatter plots giving all pairwise comparisons of four pos-1 replicate delta-frequency measurements." width="100%" />
<p class="filecap">
SUPP_FIG_XX_original_pos1_dfreq_rep_correlation
</p>

</div>

All six pairwise comparisons of the four *pos-1* replicates, at
read-depth cutoff 5, with Spearman ρ and n on each facet. Dashed line is
`y = x`; red is the fitted slope.

Depth cutoff 5 is used because it is the cutoff the shipped association
traits were built from: taking the mean *pos-1* delta per strain at
cutoff 5 reproduces `delta_ctrl_pos-1_T2` to `2.5e-16`, against *r* =
`0.9998` at cutoff 3 and `0.9988` at cutoff 10.

<div class="caveat">

<span class="ch">One strain needs a sentence</span>

Rows sharing a (strain, sample) key are summed before the delta is
taken. This affects JU1793 only, which appears twice per sample at every
depth cutoff and in both arms — one row carrying the real frequency and
the other `~1e-21`. That is NNLS splitting one strain’s abundance across
two identical entries in the genotype reference, not two measurements.
Summing roughly doubles JU1793’s *pos-1* delta, from 0.00857 to about
0.0171. For every other strain the collapse is a no-op.

</div>

## Figure S7 — the plate assay, validated externally

<div class="meta">

**Script** `scripts/SUPP_FIG_plate_vs_paaby_vs_pos1original.R`<br>
**Supports** the pooled phenotype, against two independent
measurements<br> **n** 111 strains (A) · 19 strains (B)

</div>

<div class="plate">

<img src="plots/SUPP_FIG_plate_vs_paaby_vs_pos1original.png" alt="Two scatter panels: manual plate score against the pooled VST phenotype, and against Paaby 2015 embryonic lethality." width="100%" />
<p class="filecap">
SUPP_FIG_plate_vs_paaby_vs_pos1original
</p>

</div>

The plate score is a six-level ordinal scale (0 = complete RNAi
response/sensitive, 5 = no response/resistant), so Spearman rank
correlation is used throughout: it requires only a monotonic
relationship, not shared units or linearity.

<div class="panel">

<span class="pl">A</span> Plate score against the pooled 2023 *pos-1*
response on the VST scale (`vst_ctrl_pos-1_T2`), the same trait the
association mapping was run on: `ρ = 0.41`, `n = 111`, `p = 7.8e-06`.
Positive is the expected direction — VST is positive for strains that
gained pool frequency under *pos-1* RNAi, and resistant strains score
high on the plate. Dashed line marks zero.

</div>

<div class="panel">

<span class="pl">B</span> Plate score against Paaby et al. 2015 mean
embryonic lethality for the *pos-1* clone, computed per well as
unhatched eggs over eggs plus larvae and averaged within a strain:
`ρ = −0.55`, `n = 19`, `p = 0.014`. Negative is the expected direction —
a low plate score means sensitive, and sensitive means high lethality.

</div>

Both agree with the plate assay in the predicted direction. **The Paaby
comparison rests on 19 shared strains and should be described as
consistent rather than as independent confirmation.**

<div class="caveat">

<span class="ch">Rebuilt on the VST trait</span>

The earlier version of panel A used a bootstrap frequency change from a
different export (`ρ = 0.42`, `n = 106`), which was on no scale used
elsewhere in the manuscript. The VST version puts it on the same scale
as every other phenotype panel.

</div>

# II. The map

<div class="scale">

genome-wide · two crosses

</div>

<div class="claim">

Pooled GWAS and F2 bulk-segregant cross mapping identify overlapping
loci, and the two cross contrasts separate general RNAi-response loci
from knockdown-specific ones.

</div>

## Figure 2 — pooled GWAS and cross mapping identify overlapping loci

<div class="meta">

**Script** `scripts/Figure2.R`<br> **GWAS** n = 84 strains · 322,010
markers<br> **Crosses** N2 × XZ1516 · JU1793 × JU2466<br> **Drawn** QTL
with peak LOD \> 100 (9 intervals)

</div>

<div class="plate">

<img src="plots/Figure2.png" alt="Mirrored Manhattan plot of the pooled GWAS across six chromosomes, with F2 cross QTL intervals overlaid as directional tracks." width="100%" />
<p class="filecap">
plots/Figure2.pdf · .png
</p>

</div>

Mirrored Manhattan of the pooled GWAS on the vst traits: the *mig-6*
response above the axis, the *pos-1* response below (`n = 84`, 322,010
markers). Grey dashed lines are Bonferroni over every marker (`6.81`);
green dotted lines divide α by the 732 effective independent tests
(`4.17`). Threshold labels sit in the chromosome I panel.

Arrowheads mark the F2 cross QTL, one per QTL at its peak marker,
coloured by cross and directional: *mig-6*-specific QTL above the axis,
the *pos-1* response below. The *mig-6* arrows come from the
*mig-6*-vs-*pos-1* contrast and the *pos-1* arrows from the
HT115-vs-*pos-1* contrast, so each side shows the comparison that
isolates the effect it is labelled with. Shaded vertical bands mark ±1
Mb around the pooled GWAS peaks.

<div class="aside">

<span class="ch">Why the intervals became arrows</span>

The arrows carry **position only** — no interval width, no peak LOD. The
previous version drew each QTL as an interval bar with opacity scaled to
LOD, which gave three visual channels to a mark that reliably carries
one.

The width in particular was not what it appeared to be. A cross interval
is often a few tens of kb against a 15–20 Mb axis, so every bar had to
be padded to a `0.09` Mb minimum just to be visible — at which point the
drawn width *was the padding*, not the interval, and the figure implied
a precision the data do not have. Widths and peak LODs are in the table
below, where a number needing three significant figures belongs.

Figure S17 is the same figure with the cross QTL removed altogether —
the mirrored Manhattan alone, tracks, labels and cross legend dropped —
written by the same script, so the arrows can be judged against the
panel without them.

</div>

| Track              | Cross         | Chr | Peak Mb | Interval    |  kb | LOD |
|:-------------------|:--------------|:----|--------:|:------------|----:|----:|
| *mig-6* vs *pos-1* | N2×XZ1516     | I   |    1.88 | 1.80–1.96   | 166 | 751 |
| *mig-6* vs *pos-1* | N2×XZ1516     | III |    0.43 | 0.28–0.59   | 304 | 165 |
| *mig-6* vs *pos-1* | N2×XZ1516     | IV  |   17.49 | 17.44–17.49 |  52 | 167 |
| *mig-6* vs *pos-1* | N2×XZ1516     | V   |   13.81 | 13.63–13.91 | 278 | 730 |
| *mig-6* vs *pos-1* | JU1793×JU2466 | X   |    5.96 | 5.66–6.24   | 582 | 361 |
| HT115 vs *pos-1*   | N2×XZ1516     | II  |    3.71 | 3.59–3.80   | 212 | 148 |
| HT115 vs *pos-1*   | N2×XZ1516     | III |   13.31 | 13.11–13.50 | 394 | 710 |
| HT115 vs *pos-1*   | JU1793×JU2466 | III |   13.78 | 13.73–13.78 |  51 | 140 |
| HT115 vs *pos-1*   | N2×XZ1516     | X   |   11.38 | 11.13–11.62 | 492 | 183 |

<div class="tnote">

Transcribed from `Figure2.R` console output; there is no interval table
in the deposit to derive this from. Only these nine are drawn — the
complete significant set is the legacy supplement
`SUPP_FIG_XX_cross_qtl_all`.

</div>

<div class="caveat">

<span class="ch">State this in the text rather than leave a reader to
notice it</span>

**No cross interval contains its matching pooled GWAS peak marker.** The
nearest correspondences are the *mig-6* peak at V:14.65 Mb against the
N2×XZ1516 chromosome V interval at 13.81 Mb (0.84 Mb apart), and the
*pos-1* peak at III:12.35 Mb against the N2×XZ1516 chromosome III
interval at 13.31 Mb (0.96 Mb) and the JU1793×JU2466 interval at 13.78
Mb (1.43 Mb). A GWAS on 84 strains is not expected to localise to a
cross interval; the claim is concordance of *locus*, not of marker.

</div>

## Figure S8 — why these cross parents

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_pooled_phenotype_ranks.R`<br>
**Supports** the choice of cross parents · n = 84 of 93 panel strains

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_pooled_phenotype_ranks.png" alt="Ranked per-strain pooled RNAi response for mig-6 and pos-1 on the vst scale, cross parents highlighted." width="100%" />
<p class="filecap">
SUPP_FIG_XX_pooled_phenotype_ranks
</p>

</div>

Per-strain pooled RNAi response for *mig-6* and *pos-1* on the vst
scale, ranked, with the cross parents highlighted. Establishes that the
crosses were built from the extremes of the pooled assay rather than
from convenience. Nine of the 93 panel strains have no vst value, so
`n = 84` and ranks are out of 84.

## Figure S9 — the expanded view behind Figure 2’s tracks

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_cross_contrast_panels.R`<br>
**Supports** Figure 2’s compressed tracks<br> **Contrasts**
HT115-vs-*pos-1* · HT115-vs-*mig-6* · *mig-6*-vs-*pos-1*

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_cross_contrast_panels.png" alt="Mirrored pooled GWAS above, then one panel per cross with three LOD contrasts overlaid: HT115-vs-pos-1, HT115-vs-mig-6, and mig-6-vs-pos-1." width="100%" />
<p class="filecap">
SUPP_FIG_XX_cross_contrast_panels
</p>

</div>

Mirrored pooled GWAS on top, then one panel per cross with all three
contrasts overlaid: each knockdown against the HT115 control — *pos-1*
in purple, *mig-6* in orange — and the difference between the two
knockdowns in blue.

Both HT115 traces high with the blue difference flat marks a **general**
RNAi-response locus — the machinery is affected whatever the target —
while one HT115 trace high with blue tracking it marks a
**knockdown-specific** one. This is the panel that justifies calling the
chromosome III locus general and the chromosome V, I and X loci
specific: at chromosome III both knockdowns respond, while chromosome V
and X move under *mig-6* with *pos-1* flat. Drawing the *mig-6* response
directly, rather than inferring it from the difference trace, is what
separates those two readings.

Figure 2 draws only the peaks above LOD 100, and only the tallest peak
per chromosome per contrast. `scripts/cross_qtl_full_summary.R` drops
both cutoffs and tabulates every peak above the genome-wide threshold of
LOD 3.57 into `plots/TABLE_cross_qtl_full.tsv`, with the model effect at
the peak, the parental allele frequency of each pool in a ±50 kb window,
and — for the three contrasts run in both crosses — whether the QTL is
shared. That table also flags which secondary peaks are separated from a
taller peak by a sub-threshold LOD trough and which are merely shoulders
of one sweep.

<table>
<thead>
<tr>
<th style="text-align:left;">
Quantity
</th>
<th style="text-align:right;">
Value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Contrasts scanned
</td>
<td style="text-align:right;">
18
</td>
</tr>
<tr>
<td style="text-align:left;">
QTL above the genome-wide threshold (LOD 3.57)
</td>
<td style="text-align:right;">
495
</td>
</tr>
<tr>
<td style="text-align:left;">
Of those, drawn in Figure 2 (LOD \> 100, tallest per chromosome)
</td>
<td style="text-align:right;">
46
</td>
</tr>
<tr>
<td style="text-align:left;">
Secondary peaks Figure 2’s MAX_PEAKS = 1 cannot call
</td>
<td style="text-align:right;">
387
</td>
</tr>
<tr>
<td style="text-align:left;">
— separated from the taller peak by a sub-threshold trough
</td>
<td style="text-align:right;">
79
</td>
</tr>
<tr>
<td style="text-align:left;">
— shoulders of one sweep, not independent QTL
</td>
<td style="text-align:right;">
308
</td>
</tr>
<tr>
<td style="text-align:left;">
Defensible QTL (rank 1, or a secondary with a clean trough)
</td>
<td style="text-align:right;">
187
</td>
</tr>
</tbody>
</table>

<div class="derived">

The 495 against 46 is the whole point of the table: Figure 2’s two
display constants, `LOD_MIN <- 100` and `MAX_PEAKS <- 1`, are doing more
of the filtering than the significance threshold is. Selection in these
pools is strong enough that one sweep can span most of a chromosome, so
the raw count overstates the number of independent loci — `trough.LOD`
is the column that separates a real secondary QTL from a shoulder, and
`separated` is that test applied.

</div>

Of the three contrasts run in both crosses, these are the QTL the other
cross also calls. `Match` says how the two calls were matched.
**interval** is the stricter claim: the crosses’ LOD-drop intervals
overlap. **position** means they do not overlap, but the other cross’s
scan still clears the genome-wide threshold at this peak’s own position
— the `Other LOD here` column.

<div class="derived">

Both kinds have to be shown, because interval overlap alone is too
strict here. These intervals are LOD-drop intervals, and pooled depths
put LOD in the hundreds, which collapses them: JU1793 × JU2466’s *pos-1*
peak on chromosome III has an interval 50 kb wide. A peak offset of a
few hundred kb then defeats the overlap test even when both crosses are
unambiguous at each other’s peak. On interval overlap alone the table
holds 19 QTL; adding position matches brings it to 45.

Chromosome III is the case that matters. The *sid-2* region is matched
by interval for HT115 vs *mig-6* in both directions, at LOD 940 in N2 ×
XZ1516 and 335 in JU1793 × JU2466. For HT115 vs *pos-1* it is matched
only by position: the two peaks sit at 13.31 Mb and 13.78 Mb, 469 kb
apart, so their 0.39 Mb and 0.05 Mb intervals miss each other by 234 kb
— but N2 × XZ1516 reaches LOD 710 at its peak with JU1793 × JU2466 at
LOD 88 there, and JU1793 × JU2466 reaches LOD 140 at its peak with N2 ×
XZ1516 at LOD 522 there. Requiring overlap would report the *sid-2*
locus as shared for *mig-6* and not for *pos-1*, which the scans do not
support.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Cross
</th>
<th style="text-align:left;">
Contrast
</th>
<th style="text-align:left;">
Chr
</th>
<th style="text-align:right;">
Peak (Mb)
</th>
<th style="text-align:right;">
LOD
</th>
<th style="text-align:right;">
Interval
</th>
<th style="text-align:right;">
Δfreq
</th>
<th style="text-align:left;">
Match
</th>
<th style="text-align:right;">
Other cross peak (Mb)
</th>
<th style="text-align:right;">
Other LOD
</th>
<th style="text-align:right;">
Other LOD here
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.88
</td>
<td style="text-align:right;">
728.2
</td>
<td style="text-align:right;">
1.80–1.98
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
5.0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:right;">
10.8
</td>
<td style="text-align:right;">
1.10–1.26
</td>
<td style="text-align:right;">
-0.154
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
318.9
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:right;">
100.8
</td>
<td style="text-align:right;">
3.57–3.80
</td>
<td style="text-align:right;">
0.148
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
16.1
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.47
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:right;">
13.37–13.53
</td>
<td style="text-align:right;">
0.140
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
13.72
</td>
<td style="text-align:right;">
27.0
</td>
<td style="text-align:right;">
26.3
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.24
</td>
<td style="text-align:right;">
31.2
</td>
<td style="text-align:right;">
14.77–15.24
</td>
<td style="text-align:right;">
-0.472
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
18.6
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
939.5
</td>
<td style="text-align:right;">
13.09–13.53
</td>
<td style="text-align:right;">
-0.320
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
13.45
</td>
<td style="text-align:right;">
334.7
</td>
<td style="text-align:right;">
328.3
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.45
</td>
<td style="text-align:right;">
334.7
</td>
<td style="text-align:right;">
13.23–13.66
</td>
<td style="text-align:right;">
-0.616
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
939.5
</td>
<td style="text-align:right;">
921.3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
109.2
</td>
<td style="text-align:right;">
0.32–0.67
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:right;">
69.5
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:right;">
0.45–0.90
</td>
<td style="text-align:right;">
-0.449
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
109.2
</td>
<td style="text-align:right;">
102.9
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.75
</td>
<td style="text-align:right;">
168.5
</td>
<td style="text-align:right;">
5.41–6.16
</td>
<td style="text-align:right;">
-0.538
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
47.3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.98
</td>
<td style="text-align:right;">
155.8
</td>
<td style="text-align:right;">
7.77–9.36
</td>
<td style="text-align:right;">
-0.153
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
8.40
</td>
<td style="text-align:right;">
87.8
</td>
<td style="text-align:right;">
72.8
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:right;">
17.44–17.49
</td>
<td style="text-align:right;">
0.131
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
27.7
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
30.7
</td>
<td style="text-align:right;">
0.16–0.37
</td>
<td style="text-align:right;">
-0.249
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
11.9
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.98
</td>
<td style="text-align:right;">
19.6
</td>
<td style="text-align:right;">
3.98–4.01
</td>
<td style="text-align:right;">
-0.221
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
26.7
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
6.36
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
6.31–6.42
</td>
<td style="text-align:right;">
-0.047
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
91.7
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.15
</td>
<td style="text-align:right;">
798.0
</td>
<td style="text-align:right;">
6.86–7.45
</td>
<td style="text-align:right;">
-0.836
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
7.05
</td>
<td style="text-align:right;">
39.5
</td>
<td style="text-align:right;">
51.3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.71
</td>
<td style="text-align:right;">
154.3
</td>
<td style="text-align:right;">
9.14–9.94
</td>
<td style="text-align:right;">
0.139
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
9.70
</td>
<td style="text-align:right;">
469.9
</td>
<td style="text-align:right;">
469.7
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
9.1
</td>
<td style="text-align:right;">
0.01–0.04
</td>
<td style="text-align:right;">
0.048
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
8.6
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
15.06
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
15.06–15.06
</td>
<td style="text-align:right;">
-0.158
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
15.07
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
15.4
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.71
</td>
<td style="text-align:right;">
148.4
</td>
<td style="text-align:right;">
3.59–3.80
</td>
<td style="text-align:right;">
0.239
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
3.52
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
6.0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.24
</td>
<td style="text-align:right;">
7.2
</td>
<td style="text-align:right;">
15.24–15.24
</td>
<td style="text-align:right;">
-0.252
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
6.1
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.52
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
3.17–3.87
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
3.71
</td>
<td style="text-align:right;">
148.4
</td>
<td style="text-align:right;">
125.3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
709.9
</td>
<td style="text-align:right;">
13.11–13.50
</td>
<td style="text-align:right;">
-0.287
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
88.4
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.78
</td>
<td style="text-align:right;">
139.7
</td>
<td style="text-align:right;">
13.73–13.78
</td>
<td style="text-align:right;">
-0.416
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
521.6
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.56
</td>
<td style="text-align:right;">
75.4
</td>
<td style="text-align:right;">
12.39–12.72
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
17.5
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.88
</td>
<td style="text-align:right;">
74.7
</td>
<td style="text-align:right;">
10.74–11.08
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
10.23
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
4.5
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
6.90
</td>
<td style="text-align:right;">
9.4
</td>
<td style="text-align:right;">
6.67–7.14
</td>
<td style="text-align:right;">
-0.112
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
12.0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.55
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
2.52–2.59
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
319.9
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
10.7
</td>
<td style="text-align:right;">
0.02–0.09
</td>
<td style="text-align:right;">
0.195
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
31.7
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
70.2
</td>
<td style="text-align:right;">
0.43–0.88
</td>
<td style="text-align:right;">
0.416
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
165.3
</td>
<td style="text-align:right;">
149.4
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.16
</td>
<td style="text-align:right;">
35.5
</td>
<td style="text-align:right;">
13.00–13.32
</td>
<td style="text-align:right;">
0.173
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
13.40
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
48.3
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
0.37–0.82
</td>
<td style="text-align:right;">
0.178
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
16.3
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:right;">
39.6
</td>
<td style="text-align:right;">
1.51–1.75
</td>
<td style="text-align:right;">
0.288
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
33.4
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.76
</td>
<td style="text-align:right;">
32.7
</td>
<td style="text-align:right;">
3.73–3.81
</td>
<td style="text-align:right;">
0.249
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
112.3
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.02
</td>
<td style="text-align:right;">
11.9
</td>
<td style="text-align:right;">
9.93–10.09
</td>
<td style="text-align:right;">
0.212
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
234.0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.96
</td>
<td style="text-align:right;">
361.4
</td>
<td style="text-align:right;">
5.66–6.24
</td>
<td style="text-align:right;">
0.736
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
10.5
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
43.7
</td>
<td style="text-align:right;">
0.79–0.98
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
6.5
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.62
</td>
<td style="text-align:right;">
36.6
</td>
<td style="text-align:right;">
13.56–13.69
</td>
<td style="text-align:right;">
0.121
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
11.6
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.25
</td>
<td style="text-align:right;">
28.0
</td>
<td style="text-align:right;">
9.17–9.34
</td>
<td style="text-align:right;">
-0.055
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
9.83
</td>
<td style="text-align:right;">
11.1
</td>
<td style="text-align:right;">
10.7
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
165.3
</td>
<td style="text-align:right;">
0.28–0.59
</td>
<td style="text-align:right;">
-0.162
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
70.2
</td>
<td style="text-align:right;">
67.3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.40
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
13.18–13.58
</td>
<td style="text-align:right;">
-0.016
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
13.16
</td>
<td style="text-align:right;">
35.5
</td>
<td style="text-align:right;">
31.0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
9.35
</td>
<td style="text-align:right;">
37.6
</td>
<td style="text-align:right;">
9.26–9.53
</td>
<td style="text-align:right;">
0.057
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
9.35
</td>
<td style="text-align:right;">
31.7
</td>
<td style="text-align:right;">
31.7
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.89
</td>
<td style="text-align:right;">
91.0
</td>
<td style="text-align:right;">
6.73–7.17
</td>
<td style="text-align:right;">
-0.117
</td>
<td style="text-align:left;">
interval
</td>
<td style="text-align:right;">
6.94
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
11.2
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.58
</td>
<td style="text-align:right;">
84.6
</td>
<td style="text-align:right;">
13.36–13.69
</td>
<td style="text-align:right;">
-0.093
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
166.1
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
63.8
</td>
<td style="text-align:right;">
2.76–3.13
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:left;">
position
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
30.1
</td>
</tr>
</tbody>
</table>

The complete table follows: every peak above the genome-wide threshold,
in every contrast, with both effect-size scales and the parental
frequencies they come from. `f1.a` and `f1.b` are parent 1’s frequency
in the first and second pool the contrast names, each pooled over ±50 kb
and weighted by read depth; `Δfreq` is their difference. Parent 1 is N2
in one cross and JU1793 in the other, so the sign of `Δfreq` is not
comparable between crosses without saying which parent is which. One QTL
has no frequency: its peak falls in a single-marker gap with no coverage
in either pool.

### General RNAi-response loci against target-specific ones

`scripts/cross_qtl_condition_sharing.R` uses a feature of the cross
design that Figure 2 does not: N2 × XZ1516 ran five knockdown pools
against one HT115 control — *pos-1*, *mig-6*, *par-1*, *rpn-12*, *vha-5*
— so “general” can be counted over targets rather than inferred from a
flat difference trace. The JU cross ran two, which is its ceiling.

<div class="derived">

The classification is made on **effect size**, not significance. Pooled
depths here put LOD in the hundreds, so at the genome-wide threshold
nearly every target is significant at nearly every locus: doing it that
way called 55 of 61 N2 × XZ1516 loci general, 19 in all five targets,
which is not credible. Two further facts fix the approach. Every
contrast in a cross shares one HT115 pool, so drift in that pool
imitates a general locus. And the timepoint-1 replicate — its own HT115,
*mig-6*, *par-1* and *rpn-12* pools — shows that replicate agreement of
the frequency shift depends sharply on the target: *mig-6* r = 0.975,
*rpn-12* r = 0.683, *par-1* r = 0.125. A *par-1* “response” largely does
not reproduce.

A target counts as responding when its pool differs from HT115 by
\|Δfreq\| ≥ 0.10 with LOD above threshold. A locus is general when ≥ 3
targets respond **in the same direction**. Requiring equal magnitude
instead is the wrong test — *mig-6* shifts about three times as far as
*par-1* at the chromosome V loci — and it called nothing general at all,
chromosome III included.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Locus
</th>
<th style="text-align:right;">
Interval
</th>
<th style="text-align:left;">
Targets responding
</th>
<th style="text-align:left;">
Direction
</th>
<th style="text-align:right;">
mig-6 fold
</th>
<th style="text-align:left;">
mig-6 dominant
</th>
<th style="text-align:right;">
Top LOD
</th>
<th style="text-align:left;">
Shared across crosses
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:right;">
13.03–13.54
</td>
<td style="text-align:left;">
mig6, pos1, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
1.11
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
939.5
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:9.28
</td>
<td style="text-align:right;">
9.23–9.28
</td>
<td style="text-align:left;">
mig6, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
1.43
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
327.7
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:right;">
10.28–11.08
</td>
<td style="text-align:left;">
mig6, par1, pos1, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
806.1
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:12.21
</td>
<td style="text-align:right;">
11.65–12.60
</td>
<td style="text-align:left;">
mig6, par1, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
1.68
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
725.1
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:13.81
</td>
<td style="text-align:right;">
12.98–15.26
</td>
<td style="text-align:left;">
mig6, par1, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
611.1
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:16.26
</td>
<td style="text-align:right;">
15.60–16.44
</td>
<td style="text-align:left;">
mig6, par1, rpn12, vha5
</td>
<td style="text-align:left;">
parent1
</td>
<td style="text-align:right;">
2.46
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
398.9
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
</tbody>
</table>

Six loci are general, and they split into two kinds. **Chromosome III at
13.31 Mb** — the *sid-2* interval — responds in four of five targets
with *mig-6* only 1.11× the next target, so every knockdown moves the
locus by a comparable amount. That is what a general RNAi-response locus
should look like, and it is matched in the JU cross, where both targets
tested respond (Δfreq 0.616 for *mig-6*, 0.452 for *pos-1*). The five
chromosome V loci are also concordant across three to five targets, but
*mig-6* runs 1.4–2.5× the next target and four of them are flagged
`mig6.dominant`: general in direction, carried mostly by *mig-6*.

<table>
<thead>
<tr>
<th style="text-align:left;">
Locus
</th>
<th style="text-align:right;">
Interval
</th>
<th style="text-align:right;">
Δfreq mig-6
</th>
<th style="text-align:right;">
Δfreq replicate 1
</th>
<th style="text-align:left;">
Distinct from every other target
</th>
<th style="text-align:right;">
Top LOD
</th>
<th style="text-align:left;">
Shared across crosses
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JU1793xJU2466 I:11.30
</td>
<td style="text-align:right;">
11.30–11.37
</td>
<td style="text-align:right;">
-0.301
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
24.1
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 II:0.02
</td>
<td style="text-align:right;">
0.02–0.51
</td>
<td style="text-align:right;">
0.289
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
20.4
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 II:2.02
</td>
<td style="text-align:right;">
2.02–2.02
</td>
<td style="text-align:right;">
0.143
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
17.4
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 II:10.68
</td>
<td style="text-align:right;">
10.12–10.68
</td>
<td style="text-align:right;">
0.177
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
20.5
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 II:12.25
</td>
<td style="text-align:right;">
11.69–12.25
</td>
<td style="text-align:right;">
0.140
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
23.7
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 III:0.68
</td>
<td style="text-align:right;">
0.45–0.90
</td>
<td style="text-align:right;">
0.449
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 III:9.36
</td>
<td style="text-align:right;">
9.26–9.48
</td>
<td style="text-align:right;">
0.311
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
54.1
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 V:0.28
</td>
<td style="text-align:right;">
0.16–0.37
</td>
<td style="text-align:right;">
0.249
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
30.7
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 V:1.64
</td>
<td style="text-align:right;">
1.47–1.73
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
28.7
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 V:2.94
</td>
<td style="text-align:right;">
2.90–2.98
</td>
<td style="text-align:right;">
0.185
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
29.8
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 V:3.98
</td>
<td style="text-align:right;">
3.98–4.01
</td>
<td style="text-align:right;">
0.221
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
19.6
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:5.86
</td>
<td style="text-align:right;">
5.64–5.86
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
663.2
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:8.45
</td>
<td style="text-align:right;">
8.15–8.68
</td>
<td style="text-align:right;">
0.742
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
569.1
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:0.80
</td>
<td style="text-align:right;">
0.01–0.80
</td>
<td style="text-align:right;">
-0.232
</td>
<td style="text-align:right;">
-0.215
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
249.2
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 II:13.47
</td>
<td style="text-align:right;">
12.66–13.53
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
-0.170
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 IV:12.54
</td>
<td style="text-align:right;">
12.39–12.72
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:right;">
0.015
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
80.5
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 IV:17.49
</td>
<td style="text-align:right;">
16.85–17.49
</td>
<td style="text-align:right;">
-0.131
</td>
<td style="text-align:right;">
-0.097
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:left;">
FALSE
</td>
</tr>
</tbody>
</table>

Seventeen loci are *mig-6*-specific — the only target that responds.
Four are in N2 × XZ1516, and three of those four reproduce in the
timepoint-1 replicate by sign and rough magnitude (I:0.80, −0.232
against −0.215; II:13.47, −0.140 against −0.170; IV:17.49, −0.131
against −0.097); IV:12.54 does not (+0.104 against +0.015) and should be
treated as unsupported. The thirteen in the JU cross carry weaker
evidence by construction: with only *pos-1* and *mig-6* run, “specific”
there means “not *pos-1*”, not “not any other target”.

<div class="derived">

The two classifications are not equally robust. Sweeping the response
threshold, *mig-6*-specific calls are stable — 13, 17, 18, 16 loci at
\|Δfreq\| ≥ 0.05, 0.10, 0.15, 0.20 — while general calls fall away as
the threshold rises: 29, 7, 2, 1. The general set therefore depends on
where the response cut is placed, and only chromosome III survives to
\|Δfreq\| ≥ 0.20.

</div>

<div style="max-height:560px;overflow:auto">

<table>
<thead>
<tr>
<th style="text-align:left;">
Cross
</th>
<th style="text-align:left;">
Contrast
</th>
<th style="text-align:left;">
Chr
</th>
<th style="text-align:right;">
Peak (Mb)
</th>
<th style="text-align:right;">
LOD
</th>
<th style="text-align:right;">
Interval
</th>
<th style="text-align:right;">
Width (kb)
</th>
<th style="text-align:right;">
Rank
</th>
<th style="text-align:right;">
β
</th>
<th style="text-align:right;">
f1.a
</th>
<th style="text-align:right;">
f1.b
</th>
<th style="text-align:right;">
Δfreq
</th>
<th style="text-align:right;">
Trough LOD
</th>
<th style="text-align:right;">
Separated
</th>
<th style="text-align:right;">
Other LOD at peak
</th>
<th style="text-align:right;">
Shared
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.68
</td>
<td style="text-align:right;">
39.1
</td>
<td style="text-align:right;">
7.31–10.30
</td>
<td style="text-align:right;">
2990
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.297
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.23
</td>
<td style="text-align:right;">
28.6
</td>
<td style="text-align:right;">
6.14–6.23
</td>
<td style="text-align:right;">
90
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.259
</td>
<td style="text-align:right;">
0.200
</td>
<td style="text-align:right;">
0.125
</td>
<td style="text-align:right;">
0.075
</td>
<td style="text-align:right;">
28.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.30
</td>
<td style="text-align:right;">
24.1
</td>
<td style="text-align:right;">
11.30–11.37
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.237
</td>
<td style="text-align:right;">
0.537
</td>
<td style="text-align:right;">
0.236
</td>
<td style="text-align:right;">
0.301
</td>
<td style="text-align:right;">
24.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
13.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
5.08
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
5.08–5.08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.173
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:right;">
0.429
</td>
<td style="text-align:right;">
0.071
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:right;">
10.8
</td>
<td style="text-align:right;">
1.10–1.26
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
0.663
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
-0.154
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
318.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.24
</td>
<td style="text-align:right;">
31.2
</td>
<td style="text-align:right;">
14.77–15.24
</td>
<td style="text-align:right;">
476
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.267
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.632
</td>
<td style="text-align:right;">
-0.472
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
18.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.72
</td>
<td style="text-align:right;">
27.0
</td>
<td style="text-align:right;">
13.26–13.72
</td>
<td style="text-align:right;">
459
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.249
</td>
<td style="text-align:right;">
0.200
</td>
<td style="text-align:right;">
0.111
</td>
<td style="text-align:right;">
0.089
</td>
<td style="text-align:right;">
27.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
50.7
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.25
</td>
<td style="text-align:right;">
23.7
</td>
<td style="text-align:right;">
11.69–12.25
</td>
<td style="text-align:right;">
561
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.233
</td>
<td style="text-align:right;">
0.333
</td>
<td style="text-align:right;">
0.474
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
23.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
10.8
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.68
</td>
<td style="text-align:right;">
20.5
</td>
<td style="text-align:right;">
10.12–10.68
</td>
<td style="text-align:right;">
561
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.217
</td>
<td style="text-align:right;">
0.258
</td>
<td style="text-align:right;">
0.435
</td>
<td style="text-align:right;">
-0.177
</td>
<td style="text-align:right;">
20.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
35.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
20.4
</td>
<td style="text-align:right;">
0.02–0.51
</td>
<td style="text-align:right;">
489
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.217
</td>
<td style="text-align:right;">
0.225
</td>
<td style="text-align:right;">
0.514
</td>
<td style="text-align:right;">
-0.289
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.45
</td>
<td style="text-align:right;">
334.7
</td>
<td style="text-align:right;">
13.23–13.66
</td>
<td style="text-align:right;">
430
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.655
</td>
<td style="text-align:right;">
0.349
</td>
<td style="text-align:right;">
0.964
</td>
<td style="text-align:right;">
-0.616
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
921.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
11.62
</td>
<td style="text-align:right;">
117.1
</td>
<td style="text-align:right;">
11.43–11.79
</td>
<td style="text-align:right;">
357
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.464
</td>
<td style="text-align:right;">
0.393
</td>
<td style="text-align:right;">
0.841
</td>
<td style="text-align:right;">
-0.448
</td>
<td style="text-align:right;">
89.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
90.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:right;">
0.45–0.90
</td>
<td style="text-align:right;">
454
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.398
</td>
<td style="text-align:right;">
0.282
</td>
<td style="text-align:right;">
0.731
</td>
<td style="text-align:right;">
-0.449
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
102.9
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
9.36
</td>
<td style="text-align:right;">
54.1
</td>
<td style="text-align:right;">
9.26–9.48
</td>
<td style="text-align:right;">
214
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.346
</td>
<td style="text-align:right;">
0.420
</td>
<td style="text-align:right;">
0.731
</td>
<td style="text-align:right;">
-0.311
</td>
<td style="text-align:right;">
37.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.75
</td>
<td style="text-align:right;">
168.5
</td>
<td style="text-align:right;">
5.41–6.16
</td>
<td style="text-align:right;">
749
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.555
</td>
<td style="text-align:right;">
0.287
</td>
<td style="text-align:right;">
0.825
</td>
<td style="text-align:right;">
-0.538
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
47.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.17
</td>
<td style="text-align:right;">
144.5
</td>
<td style="text-align:right;">
7.17–7.36
</td>
<td style="text-align:right;">
185
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.526
</td>
<td style="text-align:right;">
0.286
</td>
<td style="text-align:right;">
0.769
</td>
<td style="text-align:right;">
-0.483
</td>
<td style="text-align:right;">
144.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
154.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
3.66
</td>
<td style="text-align:right;">
141.9
</td>
<td style="text-align:right;">
3.56–3.78
</td>
<td style="text-align:right;">
217
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.522
</td>
<td style="text-align:right;">
0.215
</td>
<td style="text-align:right;">
0.710
</td>
<td style="text-align:right;">
-0.495
</td>
<td style="text-align:right;">
121.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
13.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.40
</td>
<td style="text-align:right;">
87.8
</td>
<td style="text-align:right;">
8.40–8.50
</td>
<td style="text-align:right;">
96
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.434
</td>
<td style="text-align:right;">
0.227
</td>
<td style="text-align:right;">
0.641
</td>
<td style="text-align:right;">
-0.414
</td>
<td style="text-align:right;">
87.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
153.2
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.02
</td>
<td style="text-align:right;">
72.7
</td>
<td style="text-align:right;">
9.69–10.38
</td>
<td style="text-align:right;">
696
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.400
</td>
<td style="text-align:right;">
0.219
</td>
<td style="text-align:right;">
0.768
</td>
<td style="text-align:right;">
-0.549
</td>
<td style="text-align:right;">
67.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
117.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
30.7
</td>
<td style="text-align:right;">
0.16–0.37
</td>
<td style="text-align:right;">
213
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.274
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:right;">
0.660
</td>
<td style="text-align:right;">
-0.249
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
11.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
2.94
</td>
<td style="text-align:right;">
29.8
</td>
<td style="text-align:right;">
2.90–2.98
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.270
</td>
<td style="text-align:right;">
0.308
</td>
<td style="text-align:right;">
0.493
</td>
<td style="text-align:right;">
-0.185
</td>
<td style="text-align:right;">
18.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.64
</td>
<td style="text-align:right;">
28.7
</td>
<td style="text-align:right;">
1.47–1.73
</td>
<td style="text-align:right;">
256
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.266
</td>
<td style="text-align:right;">
0.411
</td>
<td style="text-align:right;">
0.666
</td>
<td style="text-align:right;">
-0.255
</td>
<td style="text-align:right;">
18.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.98
</td>
<td style="text-align:right;">
19.6
</td>
<td style="text-align:right;">
3.98–4.01
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.221
</td>
<td style="text-align:right;">
0.336
</td>
<td style="text-align:right;">
0.557
</td>
<td style="text-align:right;">
-0.221
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
26.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
6.36
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
6.31–6.42
</td>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
0.284
</td>
<td style="text-align:right;">
0.332
</td>
<td style="text-align:right;">
-0.047
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
91.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.15
</td>
<td style="text-align:right;">
798.0
</td>
<td style="text-align:right;">
6.86–7.45
</td>
<td style="text-align:right;">
591
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.823
</td>
<td style="text-align:right;">
0.085
</td>
<td style="text-align:right;">
0.921
</td>
<td style="text-align:right;">
-0.836
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.86
</td>
<td style="text-align:right;">
663.2
</td>
<td style="text-align:right;">
5.75–5.86
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.798
</td>
<td style="text-align:right;">
0.079
</td>
<td style="text-align:right;">
0.896
</td>
<td style="text-align:right;">
-0.817
</td>
<td style="text-align:right;">
663.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
31.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.45
</td>
<td style="text-align:right;">
569.1
</td>
<td style="text-align:right;">
8.45–8.68
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.775
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
0.841
</td>
<td style="text-align:right;">
-0.742
</td>
<td style="text-align:right;">
569.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
132.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.70
</td>
<td style="text-align:right;">
469.9
</td>
<td style="text-align:right;">
9.70–9.85
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.744
</td>
<td style="text-align:right;">
0.073
</td>
<td style="text-align:right;">
0.826
</td>
<td style="text-align:right;">
-0.753
</td>
<td style="text-align:right;">
469.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
154.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
10.85
</td>
<td style="text-align:right;">
313.8
</td>
<td style="text-align:right;">
10.85–11.09
</td>
<td style="text-align:right;">
243
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.670
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
0.747
</td>
<td style="text-align:right;">
-0.646
</td>
<td style="text-align:right;">
313.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
71.5
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
14.4
</td>
<td style="text-align:right;">
0.73–1.08
</td>
<td style="text-align:right;">
352
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.180
</td>
<td style="text-align:right;">
0.610
</td>
<td style="text-align:right;">
0.795
</td>
<td style="text-align:right;">
-0.185
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.11
</td>
<td style="text-align:right;">
7.8
</td>
<td style="text-align:right;">
2.08–2.16
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
0.837
</td>
<td style="text-align:right;">
0.951
</td>
<td style="text-align:right;">
-0.115
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
15.06
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
15.06–15.06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.101
</td>
<td style="text-align:right;">
0.270
</td>
<td style="text-align:right;">
0.429
</td>
<td style="text-align:right;">
-0.158
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
15.4
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.15
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
8.78–9.53
</td>
<td style="text-align:right;">
744
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:right;">
0.471
</td>
<td style="text-align:right;">
0.278
</td>
<td style="text-align:right;">
0.193
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.24
</td>
<td style="text-align:right;">
7.2
</td>
<td style="text-align:right;">
15.24–15.24
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.132
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.412
</td>
<td style="text-align:right;">
-0.252
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.52
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
3.17–3.87
</td>
<td style="text-align:right;">
690
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.125
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
0.375
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
125.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
2.02
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
2.02–2.02
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
0.357
</td>
<td style="text-align:right;">
0.273
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
16.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.78
</td>
<td style="text-align:right;">
139.7
</td>
<td style="text-align:right;">
13.73–13.78
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.525
</td>
<td style="text-align:right;">
0.367
</td>
<td style="text-align:right;">
0.782
</td>
<td style="text-align:right;">
-0.416
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
521.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.73
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
12.68–12.73
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.355
</td>
<td style="text-align:right;">
0.358
</td>
<td style="text-align:right;">
0.755
</td>
<td style="text-align:right;">
-0.397
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
460.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
11.68
</td>
<td style="text-align:right;">
18.7
</td>
<td style="text-align:right;">
11.61–11.68
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.225
</td>
<td style="text-align:right;">
0.396
</td>
<td style="text-align:right;">
0.657
</td>
<td style="text-align:right;">
-0.261
</td>
<td style="text-align:right;">
18.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
60.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.42
</td>
<td style="text-align:right;">
56.7
</td>
<td style="text-align:right;">
5.10–7.29
</td>
<td style="text-align:right;">
2191
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.370
</td>
<td style="text-align:right;">
0.180
</td>
<td style="text-align:right;">
0.609
</td>
<td style="text-align:right;">
-0.428
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.10
</td>
<td style="text-align:right;">
44.6
</td>
<td style="text-align:right;">
3.92–4.10
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.327
</td>
<td style="text-align:right;">
0.143
</td>
<td style="text-align:right;">
0.357
</td>
<td style="text-align:right;">
-0.214
</td>
<td style="text-align:right;">
44.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.30
</td>
<td style="text-align:right;">
43.5
</td>
<td style="text-align:right;">
8.30–8.58
</td>
<td style="text-align:right;">
282
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.328
</td>
<td style="text-align:right;">
0.189
</td>
<td style="text-align:right;">
0.592
</td>
<td style="text-align:right;">
-0.403
</td>
<td style="text-align:right;">
43.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
2.92
</td>
<td style="text-align:right;">
32.8
</td>
<td style="text-align:right;">
2.82–2.92
</td>
<td style="text-align:right;">
95
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.277
</td>
<td style="text-align:right;">
0.198
</td>
<td style="text-align:right;">
0.440
</td>
<td style="text-align:right;">
-0.242
</td>
<td style="text-align:right;">
32.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
7.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.61
</td>
<td style="text-align:right;">
31.5
</td>
<td style="text-align:right;">
9.61–9.82
</td>
<td style="text-align:right;">
213
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.281
</td>
<td style="text-align:right;">
0.300
</td>
<td style="text-align:right;">
0.596
</td>
<td style="text-align:right;">
-0.296
</td>
<td style="text-align:right;">
31.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
23.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.53
</td>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
14.26–14.79
</td>
<td style="text-align:right;">
533
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.23
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
9.82–10.94
</td>
<td style="text-align:right;">
1120
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.093
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
0.129
</td>
<td style="text-align:right;">
0.126
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
55.5
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.91
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
12.50–13.25
</td>
<td style="text-align:right;">
749
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.082
</td>
<td style="text-align:right;">
0.172
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.074
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
8.53
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
7.97–8.80
</td>
<td style="text-align:right;">
835
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.089
</td>
<td style="text-align:right;">
0.231
</td>
<td style="text-align:right;">
0.141
</td>
<td style="text-align:right;">
0.089
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
6.90
</td>
<td style="text-align:right;">
9.4
</td>
<td style="text-align:right;">
6.67–7.14
</td>
<td style="text-align:right;">
468
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.118
</td>
<td style="text-align:right;">
0.081
</td>
<td style="text-align:right;">
0.193
</td>
<td style="text-align:right;">
-0.112
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
12.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
10.00
</td>
<td style="text-align:right;">
6.8
</td>
<td style="text-align:right;">
9.74–10.26
</td>
<td style="text-align:right;">
521
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.103
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
0.183
</td>
<td style="text-align:right;">
-0.115
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
150.4
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.15
</td>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
8.15–8.31
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.095
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
0.177
</td>
<td style="text-align:right;">
-0.074
</td>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
111.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.67
</td>
<td style="text-align:right;">
4.3
</td>
<td style="text-align:right;">
5.64–5.67
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.077
</td>
<td style="text-align:right;">
0.088
</td>
<td style="text-align:right;">
0.152
</td>
<td style="text-align:right;">
-0.065
</td>
<td style="text-align:right;">
4.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
6.8
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.65
</td>
<td style="text-align:right;">
16.6
</td>
<td style="text-align:right;">
7.18–8.25
</td>
<td style="text-align:right;">
1070
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.202
</td>
<td style="text-align:right;">
0.182
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:right;">
-0.318
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
10.17
</td>
<td style="text-align:right;">
16.1
</td>
<td style="text-align:right;">
9.61–10.80
</td>
<td style="text-align:right;">
1190
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.200
</td>
<td style="text-align:right;">
0.172
</td>
<td style="text-align:right;">
0.227
</td>
<td style="text-align:right;">
-0.055
</td>
<td style="text-align:right;">
14.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.14
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
6.09–6.14
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.174
</td>
<td style="text-align:right;">
0.345
</td>
<td style="text-align:right;">
0.400
</td>
<td style="text-align:right;">
-0.055
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.80
</td>
<td style="text-align:right;">
7.6
</td>
<td style="text-align:right;">
11.80–11.84
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.141
</td>
<td style="text-align:right;">
0.234
</td>
<td style="text-align:right;">
0.409
</td>
<td style="text-align:right;">
-0.175
</td>
<td style="text-align:right;">
7.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
79.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.55
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
2.52–2.59
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.085
</td>
<td style="text-align:right;">
0.805
</td>
<td style="text-align:right;">
0.899
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
319.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.36
</td>
<td style="text-align:right;">
15.6
</td>
<td style="text-align:right;">
10.84–12.07
</td>
<td style="text-align:right;">
1230
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.199
</td>
<td style="text-align:right;">
0.518
</td>
<td style="text-align:right;">
0.226
</td>
<td style="text-align:right;">
0.292
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.16
</td>
<td style="text-align:right;">
12.6
</td>
<td style="text-align:right;">
13.16–13.38
</td>
<td style="text-align:right;">
224
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.183
</td>
<td style="text-align:right;">
0.583
</td>
<td style="text-align:right;">
0.350
</td>
<td style="text-align:right;">
0.233
</td>
<td style="text-align:right;">
12.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.83
</td>
<td style="text-align:right;">
11.1
</td>
<td style="text-align:right;">
7.73–9.83
</td>
<td style="text-align:right;">
2099
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.170
</td>
<td style="text-align:right;">
0.463
</td>
<td style="text-align:right;">
0.288
</td>
<td style="text-align:right;">
0.176
</td>
<td style="text-align:right;">
11.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
10.7
</td>
<td style="text-align:right;">
0.02–0.09
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.168
</td>
<td style="text-align:right;">
0.514
</td>
<td style="text-align:right;">
0.318
</td>
<td style="text-align:right;">
0.195
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
31.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
6.73
</td>
<td style="text-align:right;">
9.2
</td>
<td style="text-align:right;">
6.18–6.73
</td>
<td style="text-align:right;">
555
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.156
</td>
<td style="text-align:right;">
0.399
</td>
<td style="text-align:right;">
0.321
</td>
<td style="text-align:right;">
0.078
</td>
<td style="text-align:right;">
9.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
70.2
</td>
<td style="text-align:right;">
0.43–0.88
</td>
<td style="text-align:right;">
447
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
0.729
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:right;">
0.416
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
149.4
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.16
</td>
<td style="text-align:right;">
35.5
</td>
<td style="text-align:right;">
13.00–13.32
</td>
<td style="text-align:right;">
311
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.218
</td>
<td style="text-align:right;">
0.949
</td>
<td style="text-align:right;">
0.776
</td>
<td style="text-align:right;">
0.173
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
48.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
9.35
</td>
<td style="text-align:right;">
31.7
</td>
<td style="text-align:right;">
9.25–9.44
</td>
<td style="text-align:right;">
187
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.279
</td>
<td style="text-align:right;">
0.723
</td>
<td style="text-align:right;">
0.460
</td>
<td style="text-align:right;">
0.263
</td>
<td style="text-align:right;">
8.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
37.6
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
11.49
</td>
<td style="text-align:right;">
31.3
</td>
<td style="text-align:right;">
11.38–11.62
</td>
<td style="text-align:right;">
237
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.250
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.614
</td>
<td style="text-align:right;">
0.247
</td>
<td style="text-align:right;">
8.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
1.89
</td>
<td style="text-align:right;">
16.7
</td>
<td style="text-align:right;">
1.89–1.94
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.209
</td>
<td style="text-align:right;">
0.450
</td>
<td style="text-align:right;">
0.280
</td>
<td style="text-align:right;">
0.170
</td>
<td style="text-align:right;">
16.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
74.7
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
3.64
</td>
<td style="text-align:right;">
17.5
</td>
<td style="text-align:right;">
3.57–3.69
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.214
</td>
<td style="text-align:right;">
0.707
</td>
<td style="text-align:right;">
0.529
</td>
<td style="text-align:right;">
0.178
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
0.37–0.82
</td>
<td style="text-align:right;">
446
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.204
</td>
<td style="text-align:right;">
0.486
</td>
<td style="text-align:right;">
0.308
</td>
<td style="text-align:right;">
0.178
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
16.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.72
</td>
<td style="text-align:right;">
15.0
</td>
<td style="text-align:right;">
5.48–5.91
</td>
<td style="text-align:right;">
429
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.189
</td>
<td style="text-align:right;">
0.779
</td>
<td style="text-align:right;">
0.609
</td>
<td style="text-align:right;">
0.170
</td>
<td style="text-align:right;">
9.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
31.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.94
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
6.94–7.17
</td>
<td style="text-align:right;">
225
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.164
</td>
<td style="text-align:right;">
0.907
</td>
<td style="text-align:right;">
0.679
</td>
<td style="text-align:right;">
0.227
</td>
<td style="text-align:right;">
11.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
90.5
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.39
</td>
<td style="text-align:right;">
11.1
</td>
<td style="text-align:right;">
12.34–12.53
</td>
<td style="text-align:right;">
197
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.177
</td>
<td style="text-align:right;">
0.613
</td>
<td style="text-align:right;">
0.481
</td>
<td style="text-align:right;">
0.132
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:right;">
39.6
</td>
<td style="text-align:right;">
1.51–1.75
</td>
<td style="text-align:right;">
231
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.320
</td>
<td style="text-align:right;">
0.668
</td>
<td style="text-align:right;">
0.380
</td>
<td style="text-align:right;">
0.288
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
33.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.76
</td>
<td style="text-align:right;">
32.7
</td>
<td style="text-align:right;">
3.73–3.81
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:right;">
0.538
</td>
<td style="text-align:right;">
0.289
</td>
<td style="text-align:right;">
0.249
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
112.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
28.1
</td>
<td style="text-align:right;">
0.22–0.41
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.274
</td>
<td style="text-align:right;">
0.674
</td>
<td style="text-align:right;">
0.388
</td>
<td style="text-align:right;">
0.286
</td>
<td style="text-align:right;">
20.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
17.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
7.04
</td>
<td style="text-align:right;">
13.3
</td>
<td style="text-align:right;">
6.92–7.14
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.175
</td>
<td style="text-align:right;">
0.360
</td>
<td style="text-align:right;">
0.192
</td>
<td style="text-align:right;">
0.168
</td>
<td style="text-align:right;">
7.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
114.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.02
</td>
<td style="text-align:right;">
11.9
</td>
<td style="text-align:right;">
9.93–10.09
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.153
</td>
<td style="text-align:right;">
0.345
</td>
<td style="text-align:right;">
0.133
</td>
<td style="text-align:right;">
0.212
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
234.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.96
</td>
<td style="text-align:right;">
361.4
</td>
<td style="text-align:right;">
5.66–6.24
</td>
<td style="text-align:right;">
582
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.713
</td>
<td style="text-align:right;">
0.906
</td>
<td style="text-align:right;">
0.169
</td>
<td style="text-align:right;">
0.736
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
10.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.24
</td>
<td style="text-align:right;">
358.5
</td>
<td style="text-align:right;">
7.24–7.70
</td>
<td style="text-align:right;">
450
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.710
</td>
<td style="text-align:right;">
0.905
</td>
<td style="text-align:right;">
0.214
</td>
<td style="text-align:right;">
0.691
</td>
<td style="text-align:right;">
329.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.70
</td>
<td style="text-align:right;">
297.3
</td>
<td style="text-align:right;">
8.70–9.24
</td>
<td style="text-align:right;">
535
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.678
</td>
<td style="text-align:right;">
0.892
</td>
<td style="text-align:right;">
0.208
</td>
<td style="text-align:right;">
0.684
</td>
<td style="text-align:right;">
297.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
10.24
</td>
<td style="text-align:right;">
207.5
</td>
<td style="text-align:right;">
10.24–10.38
</td>
<td style="text-align:right;">
142
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.611
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.168
</td>
<td style="text-align:right;">
0.693
</td>
<td style="text-align:right;">
207.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466
</td>
<td style="text-align:left;">
mig6 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
12.72
</td>
<td style="text-align:right;">
198.1
</td>
<td style="text-align:right;">
11.88–12.93
</td>
<td style="text-align:right;">
1053
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.601
</td>
<td style="text-align:right;">
0.787
</td>
<td style="text-align:right;">
0.222
</td>
<td style="text-align:right;">
0.564
</td>
<td style="text-align:right;">
183.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
22.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.88
</td>
<td style="text-align:right;">
728.2
</td>
<td style="text-align:right;">
1.80–1.98
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.386
</td>
<td style="text-align:right;">
0.725
</td>
<td style="text-align:right;">
0.315
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
249.2
</td>
<td style="text-align:right;">
0.01–0.80
</td>
<td style="text-align:right;">
795
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.238
</td>
<td style="text-align:right;">
0.635
</td>
<td style="text-align:right;">
0.403
</td>
<td style="text-align:right;">
0.232
</td>
<td style="text-align:right;">
249.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.99
</td>
<td style="text-align:right;">
153.8
</td>
<td style="text-align:right;">
2.99–3.00
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.833
</td>
<td style="text-align:right;">
0.691
</td>
<td style="text-align:right;">
0.142
</td>
<td style="text-align:right;">
153.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:right;">
100.8
</td>
<td style="text-align:right;">
3.57–3.80
</td>
<td style="text-align:right;">
220
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.154
</td>
<td style="text-align:right;">
0.651
</td>
<td style="text-align:right;">
0.503
</td>
<td style="text-align:right;">
0.148
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
16.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.47
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:right;">
13.37–13.53
</td>
<td style="text-align:right;">
151
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.142
</td>
<td style="text-align:right;">
0.415
</td>
<td style="text-align:right;">
0.275
</td>
<td style="text-align:right;">
0.140
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
26.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.40
</td>
<td style="text-align:right;">
41.8
</td>
<td style="text-align:right;">
11.27–11.52
</td>
<td style="text-align:right;">
241
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.094
</td>
<td style="text-align:right;">
0.349
</td>
<td style="text-align:right;">
0.253
</td>
<td style="text-align:right;">
0.096
</td>
<td style="text-align:right;">
8.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
21.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.80
</td>
<td style="text-align:right;">
26.6
</td>
<td style="text-align:right;">
4.80–4.89
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.079
</td>
<td style="text-align:right;">
0.435
</td>
<td style="text-align:right;">
0.371
</td>
<td style="text-align:right;">
0.064
</td>
<td style="text-align:right;">
23.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
15.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.95
</td>
<td style="text-align:right;">
26.4
</td>
<td style="text-align:right;">
9.89–10.07
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.076
</td>
<td style="text-align:right;">
0.366
</td>
<td style="text-align:right;">
0.282
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
23.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
19.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
939.5
</td>
<td style="text-align:right;">
13.09–13.53
</td>
<td style="text-align:right;">
436
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.314
</td>
<td style="text-align:right;">
0.677
</td>
<td style="text-align:right;">
0.996
</td>
<td style="text-align:right;">
-0.320
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
328.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.08
</td>
<td style="text-align:right;">
259.4
</td>
<td style="text-align:right;">
12.04–12.08
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.153
</td>
<td style="text-align:right;">
0.818
</td>
<td style="text-align:right;">
0.958
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
259.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
92.8
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
109.2
</td>
<td style="text-align:right;">
0.32–0.67
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.156
</td>
<td style="text-align:right;">
0.516
</td>
<td style="text-align:right;">
0.678
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
69.5
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.98
</td>
<td style="text-align:right;">
155.8
</td>
<td style="text-align:right;">
7.77–9.36
</td>
<td style="text-align:right;">
1587
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.148
</td>
<td style="text-align:right;">
0.108
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
-0.153
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
72.8
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.77
</td>
<td style="text-align:right;">
132.8
</td>
<td style="text-align:right;">
6.72–6.77
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.144
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.250
</td>
<td style="text-align:right;">
-0.135
</td>
<td style="text-align:right;">
132.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
150.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
11.17
</td>
<td style="text-align:right;">
111.5
</td>
<td style="text-align:right;">
10.78–11.42
</td>
<td style="text-align:right;">
631
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.127
</td>
<td style="text-align:right;">
0.124
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
-0.137
</td>
<td style="text-align:right;">
101.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
55.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
93.6
</td>
<td style="text-align:right;">
17.44–17.49
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.137
</td>
<td style="text-align:right;">
0.396
</td>
<td style="text-align:right;">
0.266
</td>
<td style="text-align:right;">
0.131
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
27.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.42
</td>
<td style="text-align:right;">
68.6
</td>
<td style="text-align:right;">
12.42–12.50
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
0.166
</td>
<td style="text-align:right;">
0.280
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
68.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
70.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.76
</td>
<td style="text-align:right;">
806.1
</td>
<td style="text-align:right;">
10.52–11.02
</td>
<td style="text-align:right;">
491
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.397
</td>
<td style="text-align:right;">
0.389
</td>
<td style="text-align:right;">
0.773
</td>
<td style="text-align:right;">
-0.384
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.21
</td>
<td style="text-align:right;">
725.1
</td>
<td style="text-align:right;">
12.02–12.60
</td>
<td style="text-align:right;">
583
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.385
</td>
<td style="text-align:right;">
0.341
</td>
<td style="text-align:right;">
0.728
</td>
<td style="text-align:right;">
-0.387
</td>
<td style="text-align:right;">
682.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.81
</td>
<td style="text-align:right;">
611.1
</td>
<td style="text-align:right;">
13.61–13.94
</td>
<td style="text-align:right;">
330
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.359
</td>
<td style="text-align:right;">
0.357
</td>
<td style="text-align:right;">
0.732
</td>
<td style="text-align:right;">
-0.374
</td>
<td style="text-align:right;">
525.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.94
</td>
<td style="text-align:right;">
480.4
</td>
<td style="text-align:right;">
14.94–15.26
</td>
<td style="text-align:right;">
320
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.320
</td>
<td style="text-align:right;">
0.262
</td>
<td style="text-align:right;">
0.587
</td>
<td style="text-align:right;">
-0.325
</td>
<td style="text-align:right;">
454.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
16.26
</td>
<td style="text-align:right;">
398.9
</td>
<td style="text-align:right;">
16.26–16.44
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.296
</td>
<td style="text-align:right;">
0.302
</td>
<td style="text-align:right;">
0.626
</td>
<td style="text-align:right;">
-0.323
</td>
<td style="text-align:right;">
398.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.71
</td>
<td style="text-align:right;">
154.3
</td>
<td style="text-align:right;">
9.14–9.94
</td>
<td style="text-align:right;">
796
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.124
</td>
<td style="text-align:right;">
0.941
</td>
<td style="text-align:right;">
0.802
</td>
<td style="text-align:right;">
0.139
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
469.7
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.14
</td>
<td style="text-align:right;">
105.4
</td>
<td style="text-align:right;">
8.05–8.14
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.095
</td>
<td style="text-align:right;">
0.951
</td>
<td style="text-align:right;">
0.866
</td>
<td style="text-align:right;">
0.085
</td>
<td style="text-align:right;">
105.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
610.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
10.95
</td>
<td style="text-align:right;">
63.0
</td>
<td style="text-align:right;">
10.95–10.98
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.077
</td>
<td style="text-align:right;">
0.932
</td>
<td style="text-align:right;">
0.856
</td>
<td style="text-align:right;">
0.077
</td>
<td style="text-align:right;">
63.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
307.9
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.99
</td>
<td style="text-align:right;">
51.9
</td>
<td style="text-align:right;">
11.99–12.04
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.066
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.881
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
50.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
290.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.05
</td>
<td style="text-align:right;">
39.5
</td>
<td style="text-align:right;">
7.04–7.05
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.050
</td>
<td style="text-align:right;">
0.954
</td>
<td style="text-align:right;">
0.916
</td>
<td style="text-align:right;">
0.038
</td>
<td style="text-align:right;">
39.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
789.1
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.63
</td>
<td style="text-align:right;">
83.7
</td>
<td style="text-align:right;">
9.29–9.86
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.089
</td>
<td style="text-align:right;">
0.917
</td>
<td style="text-align:right;">
0.833
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.82
</td>
<td style="text-align:right;">
76.0
</td>
<td style="text-align:right;">
7.62–8.29
</td>
<td style="text-align:right;">
666
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
0.933
</td>
<td style="text-align:right;">
0.858
</td>
<td style="text-align:right;">
0.074
</td>
<td style="text-align:right;">
71.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.58
</td>
<td style="text-align:right;">
57.3
</td>
<td style="text-align:right;">
11.46–11.82
</td>
<td style="text-align:right;">
351
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
0.880
</td>
<td style="text-align:right;">
0.826
</td>
<td style="text-align:right;">
0.054
</td>
<td style="text-align:right;">
42.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.62
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
6.51–6.62
</td>
<td style="text-align:right;">
111
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.069
</td>
<td style="text-align:right;">
0.929
</td>
<td style="text-align:right;">
0.859
</td>
<td style="text-align:right;">
0.070
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.91
</td>
<td style="text-align:right;">
45.8
</td>
<td style="text-align:right;">
13.77–14.15
</td>
<td style="text-align:right;">
376
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.707
</td>
<td style="text-align:right;">
0.615
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
35.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
8.62
</td>
<td style="text-align:right;">
63.9
</td>
<td style="text-align:right;">
8.52–8.71
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
0.343
</td>
<td style="text-align:right;">
0.447
</td>
<td style="text-align:right;">
-0.104
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.20
</td>
<td style="text-align:right;">
60.7
</td>
<td style="text-align:right;">
7.15–7.28
</td>
<td style="text-align:right;">
121
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.117
</td>
<td style="text-align:right;">
0.315
</td>
<td style="text-align:right;">
0.443
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
34.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.59
</td>
<td style="text-align:right;">
47.9
</td>
<td style="text-align:right;">
11.49–11.66
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.106
</td>
<td style="text-align:right;">
0.374
</td>
<td style="text-align:right;">
0.494
</td>
<td style="text-align:right;">
-0.120
</td>
<td style="text-align:right;">
22.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
5.94
</td>
<td style="text-align:right;">
45.7
</td>
<td style="text-align:right;">
5.84–6.04
</td>
<td style="text-align:right;">
191
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.104
</td>
<td style="text-align:right;">
0.347
</td>
<td style="text-align:right;">
0.463
</td>
<td style="text-align:right;">
-0.116
</td>
<td style="text-align:right;">
33.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.66
</td>
<td style="text-align:right;">
44.0
</td>
<td style="text-align:right;">
12.66–12.69
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.101
</td>
<td style="text-align:right;">
0.352
</td>
<td style="text-align:right;">
0.451
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
32.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
18.5
</td>
<td style="text-align:right;">
0.45–0.67
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.067
</td>
<td style="text-align:right;">
0.522
</td>
<td style="text-align:right;">
0.577
</td>
<td style="text-align:right;">
-0.056
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
2.13
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
2.06–2.21
</td>
<td style="text-align:right;">
150
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.033
</td>
<td style="text-align:right;">
0.669
</td>
<td style="text-align:right;">
0.681
</td>
<td style="text-align:right;">
-0.012
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.17
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
13.07–13.28
</td>
<td style="text-align:right;">
201
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.029
</td>
<td style="text-align:right;">
0.690
</td>
<td style="text-align:right;">
0.714
</td>
<td style="text-align:right;">
-0.024
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
3.21
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
3.21–3.25
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.022
</td>
<td style="text-align:right;">
0.840
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
-0.021
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.01
</td>
<td style="text-align:right;">
141.0
</td>
<td style="text-align:right;">
8.72–9.90
</td>
<td style="text-align:right;">
1177
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.257
</td>
<td style="text-align:right;">
-0.141
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.29
</td>
<td style="text-align:right;">
114.1
</td>
<td style="text-align:right;">
6.97–7.53
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.129
</td>
<td style="text-align:right;">
0.123
</td>
<td style="text-align:right;">
0.238
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
106.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.91
</td>
<td style="text-align:right;">
105.0
</td>
<td style="text-align:right;">
10.91–11.54
</td>
<td style="text-align:right;">
630
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.123
</td>
<td style="text-align:right;">
0.130
</td>
<td style="text-align:right;">
0.257
</td>
<td style="text-align:right;">
-0.126
</td>
<td style="text-align:right;">
105.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.96
</td>
<td style="text-align:right;">
87.1
</td>
<td style="text-align:right;">
4.84–5.08
</td>
<td style="text-align:right;">
236
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.116
</td>
<td style="text-align:right;">
0.132
</td>
<td style="text-align:right;">
0.267
</td>
<td style="text-align:right;">
-0.135
</td>
<td style="text-align:right;">
76.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.54
</td>
<td style="text-align:right;">
80.5
</td>
<td style="text-align:right;">
12.54–12.60
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.113
</td>
<td style="text-align:right;">
0.156
</td>
<td style="text-align:right;">
0.256
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
80.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.55
</td>
<td style="text-align:right;">
90.0
</td>
<td style="text-align:right;">
14.42–14.71
</td>
<td style="text-align:right;">
286
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
0.298
</td>
<td style="text-align:right;">
0.438
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
15.72
</td>
<td style="text-align:right;">
71.8
</td>
<td style="text-align:right;">
15.72–16.21
</td>
<td style="text-align:right;">
495
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
0.321
</td>
<td style="text-align:right;">
0.450
</td>
<td style="text-align:right;">
-0.130
</td>
<td style="text-align:right;">
71.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
11.76
</td>
<td style="text-align:right;">
70.4
</td>
<td style="text-align:right;">
11.65–11.98
</td>
<td style="text-align:right;">
321
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
0.339
</td>
<td style="text-align:right;">
0.458
</td>
<td style="text-align:right;">
-0.119
</td>
<td style="text-align:right;">
39.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.98
</td>
<td style="text-align:right;">
45.7
</td>
<td style="text-align:right;">
12.98–13.06
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.103
</td>
<td style="text-align:right;">
0.333
</td>
<td style="text-align:right;">
0.433
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
45.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.65
</td>
<td style="text-align:right;">
40.1
</td>
<td style="text-align:right;">
10.54–10.65
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.098
</td>
<td style="text-align:right;">
0.417
</td>
<td style="text-align:right;">
0.514
</td>
<td style="text-align:right;">
-0.098
</td>
<td style="text-align:right;">
40.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
15.73
</td>
<td style="text-align:right;">
18.4
</td>
<td style="text-align:right;">
15.50–15.85
</td>
<td style="text-align:right;">
345
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.058
</td>
<td style="text-align:right;">
0.780
</td>
<td style="text-align:right;">
0.750
</td>
<td style="text-align:right;">
0.030
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.81
</td>
<td style="text-align:right;">
15.6
</td>
<td style="text-align:right;">
5.50–5.94
</td>
<td style="text-align:right;">
435
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.025
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
0.942
</td>
<td style="text-align:right;">
0.027
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
12.43
</td>
<td style="text-align:right;">
14.5
</td>
<td style="text-align:right;">
12.23–12.65
</td>
<td style="text-align:right;">
411
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.033
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
0.909
</td>
<td style="text-align:right;">
0.025
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
4.50
</td>
<td style="text-align:right;">
13.7
</td>
<td style="text-align:right;">
4.45–4.50
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.024
</td>
<td style="text-align:right;">
0.962
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
0.017
</td>
<td style="text-align:right;">
10.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.72
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
8.57–8.98
</td>
<td style="text-align:right;">
401
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.026
</td>
<td style="text-align:right;">
0.950
</td>
<td style="text-align:right;">
0.918
</td>
<td style="text-align:right;">
0.032
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.47
</td>
<td style="text-align:right;">
29.9
</td>
<td style="text-align:right;">
13.27–13.55
</td>
<td style="text-align:right;">
275
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.072
</td>
<td style="text-align:right;">
0.731
</td>
<td style="text-align:right;">
0.799
</td>
<td style="text-align:right;">
-0.068
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
15.07
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
14.95–15.07
</td>
<td style="text-align:right;">
114
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.060
</td>
<td style="text-align:right;">
0.714
</td>
<td style="text-align:right;">
0.766
</td>
<td style="text-align:right;">
-0.053
</td>
<td style="text-align:right;">
13.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.08
</td>
<td style="text-align:right;">
12.5
</td>
<td style="text-align:right;">
10.96–11.14
</td>
<td style="text-align:right;">
177
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.032
</td>
<td style="text-align:right;">
0.896
</td>
<td style="text-align:right;">
0.931
</td>
<td style="text-align:right;">
-0.035
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
9.1
</td>
<td style="text-align:right;">
0.01–0.04
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.046
</td>
<td style="text-align:right;">
0.633
</td>
<td style="text-align:right;">
0.585
</td>
<td style="text-align:right;">
0.048
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
8.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.71
</td>
<td style="text-align:right;">
148.4
</td>
<td style="text-align:right;">
3.59–3.80
</td>
<td style="text-align:right;">
209
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.186
</td>
<td style="text-align:right;">
0.661
</td>
<td style="text-align:right;">
0.421
</td>
<td style="text-align:right;">
0.239
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.77
</td>
<td style="text-align:right;">
66.2
</td>
<td style="text-align:right;">
10.71–10.83
</td>
<td style="text-align:right;">
116
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:right;">
0.370
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
0.109
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
1.57
</td>
<td style="text-align:right;">
50.8
</td>
<td style="text-align:right;">
1.52–1.63
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
0.717
</td>
<td style="text-align:right;">
0.640
</td>
<td style="text-align:right;">
0.076
</td>
<td style="text-align:right;">
16.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.92
</td>
<td style="text-align:right;">
47.7
</td>
<td style="text-align:right;">
12.88–12.97
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
0.388
</td>
<td style="text-align:right;">
0.283
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
13.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
8.12
</td>
<td style="text-align:right;">
47.1
</td>
<td style="text-align:right;">
8.04–8.21
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.097
</td>
<td style="text-align:right;">
0.346
</td>
<td style="text-align:right;">
0.242
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:right;">
23.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
709.9
</td>
<td style="text-align:right;">
13.11–13.50
</td>
<td style="text-align:right;">
385
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.288
</td>
<td style="text-align:right;">
0.681
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
-0.287
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
88.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.11
</td>
<td style="text-align:right;">
247.4
</td>
<td style="text-align:right;">
12.08–12.11
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.153
</td>
<td style="text-align:right;">
0.812
</td>
<td style="text-align:right;">
0.965
</td>
<td style="text-align:right;">
-0.154
</td>
<td style="text-align:right;">
247.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
28.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.56
</td>
<td style="text-align:right;">
75.4
</td>
<td style="text-align:right;">
12.39–12.72
</td>
<td style="text-align:right;">
330
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.109
</td>
<td style="text-align:right;">
0.156
</td>
<td style="text-align:right;">
0.249
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
17.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.85
</td>
<td style="text-align:right;">
48.7
</td>
<td style="text-align:right;">
15.75–15.95
</td>
<td style="text-align:right;">
196
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.106
</td>
<td style="text-align:right;">
0.322
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:right;">
-0.089
</td>
<td style="text-align:right;">
15.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
10.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
13.73
</td>
<td style="text-align:right;">
38.8
</td>
<td style="text-align:right;">
13.73–13.88
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
0.233
</td>
<td style="text-align:right;">
0.308
</td>
<td style="text-align:right;">
-0.074
</td>
<td style="text-align:right;">
38.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
7.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
11.38
</td>
<td style="text-align:right;">
33.7
</td>
<td style="text-align:right;">
11.31–11.38
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.066
</td>
<td style="text-align:right;">
0.167
</td>
<td style="text-align:right;">
0.214
</td>
<td style="text-align:right;">
-0.046
</td>
<td style="text-align:right;">
33.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
23.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.32
</td>
<td style="text-align:right;">
24.6
</td>
<td style="text-align:right;">
8.99–9.54
</td>
<td style="text-align:right;">
546
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.053
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
22.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
34.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.88
</td>
<td style="text-align:right;">
74.7
</td>
<td style="text-align:right;">
10.74–11.08
</td>
<td style="text-align:right;">
341
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.133
</td>
<td style="text-align:right;">
0.390
</td>
<td style="text-align:right;">
0.518
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
65.2
</td>
<td style="text-align:right;">
1.18–1.44
</td>
<td style="text-align:right;">
255
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:right;">
0.720
</td>
<td style="text-align:right;">
0.585
</td>
<td style="text-align:right;">
0.135
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
65.0
</td>
<td style="text-align:right;">
0.01–0.18
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.117
</td>
<td style="text-align:right;">
0.735
</td>
<td style="text-align:right;">
0.634
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
48.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
45.3
</td>
<td style="text-align:right;">
2.49–2.74
</td>
<td style="text-align:right;">
249
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.885
</td>
<td style="text-align:right;">
0.872
</td>
<td style="text-align:right;">
0.014
</td>
<td style="text-align:right;">
45.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.74
</td>
<td style="text-align:right;">
34.7
</td>
<td style="text-align:right;">
3.74–3.90
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.085
</td>
<td style="text-align:right;">
0.729
</td>
<td style="text-align:right;">
0.639
</td>
<td style="text-align:right;">
0.090
</td>
<td style="text-align:right;">
34.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.38
</td>
<td style="text-align:right;">
183.4
</td>
<td style="text-align:right;">
11.13–11.62
</td>
<td style="text-align:right;">
490
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.145
</td>
<td style="text-align:right;">
0.927
</td>
<td style="text-align:right;">
0.782
</td>
<td style="text-align:right;">
0.145
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.55
</td>
<td style="text-align:right;">
159.3
</td>
<td style="text-align:right;">
9.33–9.97
</td>
<td style="text-align:right;">
637
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.126
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
0.821
</td>
<td style="text-align:right;">
0.125
</td>
<td style="text-align:right;">
125.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.61
</td>
<td style="text-align:right;">
119.8
</td>
<td style="text-align:right;">
7.50–7.91
</td>
<td style="text-align:right;">
401
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
0.946
</td>
<td style="text-align:right;">
0.852
</td>
<td style="text-align:right;">
0.093
</td>
<td style="text-align:right;">
109.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
7.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.00
</td>
<td style="text-align:right;">
109.1
</td>
<td style="text-align:right;">
12.79–13.18
</td>
<td style="text-align:right;">
386
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.113
</td>
<td style="text-align:right;">
0.916
</td>
<td style="text-align:right;">
0.802
</td>
<td style="text-align:right;">
0.114
</td>
<td style="text-align:right;">
103.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs pos1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
14.19
</td>
<td style="text-align:right;">
38.9
</td>
<td style="text-align:right;">
14.19–14.24
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.080
</td>
<td style="text-align:right;">
0.831
</td>
<td style="text-align:right;">
0.758
</td>
<td style="text-align:right;">
0.073
</td>
<td style="text-align:right;">
36.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.86
</td>
<td style="text-align:right;">
17.1
</td>
<td style="text-align:right;">
9.70–10.02
</td>
<td style="text-align:right;">
320
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.036
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.882
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.78
</td>
<td style="text-align:right;">
16.8
</td>
<td style="text-align:right;">
7.67–7.87
</td>
<td style="text-align:right;">
195
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.035
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.900
</td>
<td style="text-align:right;">
0.035
</td>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
14.08
</td>
<td style="text-align:right;">
14.3
</td>
<td style="text-align:right;">
13.95–14.18
</td>
<td style="text-align:right;">
221
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.055
</td>
<td style="text-align:right;">
0.720
</td>
<td style="text-align:right;">
0.659
</td>
<td style="text-align:right;">
0.061
</td>
<td style="text-align:right;">
7.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
12.95
</td>
<td style="text-align:right;">
13.7
</td>
<td style="text-align:right;">
12.69–12.95
</td>
<td style="text-align:right;">
260
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.050
</td>
<td style="text-align:right;">
0.796
</td>
<td style="text-align:right;">
0.743
</td>
<td style="text-align:right;">
0.053
</td>
<td style="text-align:right;">
8.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.33
</td>
<td style="text-align:right;">
12.8
</td>
<td style="text-align:right;">
11.24–11.44
</td>
<td style="text-align:right;">
200
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.037
</td>
<td style="text-align:right;">
0.893
</td>
<td style="text-align:right;">
0.858
</td>
<td style="text-align:right;">
0.035
</td>
<td style="text-align:right;">
7.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.08
</td>
<td style="text-align:right;">
27.6
</td>
<td style="text-align:right;">
12.01–12.15
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.082
</td>
<td style="text-align:right;">
0.431
</td>
<td style="text-align:right;">
0.545
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.21
</td>
<td style="text-align:right;">
26.3
</td>
<td style="text-align:right;">
13.15–13.29
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.079
</td>
<td style="text-align:right;">
0.348
</td>
<td style="text-align:right;">
0.451
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
16.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.28
</td>
<td style="text-align:right;">
18.8
</td>
<td style="text-align:right;">
15.27–15.28
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.067
</td>
<td style="text-align:right;">
0.408
</td>
<td style="text-align:right;">
0.462
</td>
<td style="text-align:right;">
-0.054
</td>
<td style="text-align:right;">
8.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.81
</td>
<td style="text-align:right;">
18.4
</td>
<td style="text-align:right;">
3.72–3.87
</td>
<td style="text-align:right;">
147
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
0.622
</td>
<td style="text-align:right;">
0.568
</td>
<td style="text-align:right;">
0.054
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
5.42
</td>
<td style="text-align:right;">
11.0
</td>
<td style="text-align:right;">
5.38–5.46
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.052
</td>
<td style="text-align:right;">
0.373
</td>
<td style="text-align:right;">
0.436
</td>
<td style="text-align:right;">
-0.063
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.29
</td>
<td style="text-align:right;">
82.4
</td>
<td style="text-align:right;">
13.03–13.54
</td>
<td style="text-align:right;">
511
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.122
</td>
<td style="text-align:right;">
0.680
</td>
<td style="text-align:right;">
0.793
</td>
<td style="text-align:right;">
-0.113
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
8.01
</td>
<td style="text-align:right;">
72.4
</td>
<td style="text-align:right;">
7.96–8.08
</td>
<td style="text-align:right;">
121
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.045
</td>
<td style="text-align:right;">
0.946
</td>
<td style="text-align:right;">
0.990
</td>
<td style="text-align:right;">
-0.043
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
6.08
</td>
<td style="text-align:right;">
69.8
</td>
<td style="text-align:right;">
5.95–6.30
</td>
<td style="text-align:right;">
345
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
0.933
</td>
<td style="text-align:right;">
0.983
</td>
<td style="text-align:right;">
-0.050
</td>
<td style="text-align:right;">
24.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
4.95
</td>
<td style="text-align:right;">
30.4
</td>
<td style="text-align:right;">
4.74–4.95
</td>
<td style="text-align:right;">
205
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.041
</td>
<td style="text-align:right;">
0.921
</td>
<td style="text-align:right;">
0.953
</td>
<td style="text-align:right;">
-0.032
</td>
<td style="text-align:right;">
30.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.03
</td>
<td style="text-align:right;">
22.0
</td>
<td style="text-align:right;">
11.99–12.03
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.053
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
0.853
</td>
<td style="text-align:right;">
-0.036
</td>
<td style="text-align:right;">
22.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.21
</td>
<td style="text-align:right;">
13.9
</td>
<td style="text-align:right;">
16.99–17.39
</td>
<td style="text-align:right;">
391
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.057
</td>
<td style="text-align:right;">
0.422
</td>
<td style="text-align:right;">
0.340
</td>
<td style="text-align:right;">
0.082
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
8.0
</td>
<td style="text-align:right;">
0.02–0.08
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.045
</td>
<td style="text-align:right;">
0.494
</td>
<td style="text-align:right;">
0.553
</td>
<td style="text-align:right;">
-0.059
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.14
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
4.07–4.20
</td>
<td style="text-align:right;">
126
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.030
</td>
<td style="text-align:right;">
0.452
</td>
<td style="text-align:right;">
0.386
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
13.73
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
13.64–13.86
</td>
<td style="text-align:right;">
221
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.031
</td>
<td style="text-align:right;">
0.237
</td>
<td style="text-align:right;">
0.218
</td>
<td style="text-align:right;">
0.019
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.26
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
6.19–6.40
</td>
<td style="text-align:right;">
206
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.023
</td>
<td style="text-align:right;">
0.182
</td>
<td style="text-align:right;">
0.147
</td>
<td style="text-align:right;">
0.035
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.76
</td>
<td style="text-align:right;">
97.7
</td>
<td style="text-align:right;">
10.58–10.94
</td>
<td style="text-align:right;">
362
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.153
</td>
<td style="text-align:right;">
0.398
</td>
<td style="text-align:right;">
0.546
</td>
<td style="text-align:right;">
-0.148
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
11.94
</td>
<td style="text-align:right;">
86.4
</td>
<td style="text-align:right;">
11.94–12.11
</td>
<td style="text-align:right;">
163
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.142
</td>
<td style="text-align:right;">
0.337
</td>
<td style="text-align:right;">
0.494
</td>
<td style="text-align:right;">
-0.157
</td>
<td style="text-align:right;">
73.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
15.85
</td>
<td style="text-align:right;">
75.7
</td>
<td style="text-align:right;">
15.60–16.15
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.132
</td>
<td style="text-align:right;">
0.339
</td>
<td style="text-align:right;">
0.465
</td>
<td style="text-align:right;">
-0.126
</td>
<td style="text-align:right;">
49.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.60
</td>
<td style="text-align:right;">
68.0
</td>
<td style="text-align:right;">
14.16–14.60
</td>
<td style="text-align:right;">
438
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
0.300
</td>
<td style="text-align:right;">
0.457
</td>
<td style="text-align:right;">
-0.157
</td>
<td style="text-align:right;">
59.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.11
</td>
<td style="text-align:right;">
53.4
</td>
<td style="text-align:right;">
13.11–13.16
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.112
</td>
<td style="text-align:right;">
0.373
</td>
<td style="text-align:right;">
0.478
</td>
<td style="text-align:right;">
-0.105
</td>
<td style="text-align:right;">
53.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.62
</td>
<td style="text-align:right;">
10.4
</td>
<td style="text-align:right;">
8.50–8.70
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.024
</td>
<td style="text-align:right;">
0.958
</td>
<td style="text-align:right;">
0.924
</td>
<td style="text-align:right;">
0.034
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.64
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
11.50–11.78
</td>
<td style="text-align:right;">
281
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.021
</td>
<td style="text-align:right;">
0.932
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
0.022
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
2.19–2.39
</td>
<td style="text-align:right;">
196
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.029
</td>
<td style="text-align:right;">
0.826
</td>
<td style="text-align:right;">
0.786
</td>
<td style="text-align:right;">
0.041
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.70
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
9.70–9.72
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.017
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.009
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.50
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
7.48–7.50
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.016
</td>
<td style="text-align:right;">
0.950
</td>
<td style="text-align:right;">
0.933
</td>
<td style="text-align:right;">
0.017
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
12.54
</td>
<td style="text-align:right;">
48.5
</td>
<td style="text-align:right;">
12.45–12.66
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.122
</td>
<td style="text-align:right;">
0.836
</td>
<td style="text-align:right;">
0.699
</td>
<td style="text-align:right;">
0.136
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.85
</td>
<td style="text-align:right;">
45.9
</td>
<td style="text-align:right;">
7.62–8.13
</td>
<td style="text-align:right;">
511
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.094
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.839
</td>
<td style="text-align:right;">
0.097
</td>
<td style="text-align:right;">
15.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.37
</td>
<td style="text-align:right;">
45.6
</td>
<td style="text-align:right;">
2.25–2.46
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.126
</td>
<td style="text-align:right;">
0.793
</td>
<td style="text-align:right;">
0.667
</td>
<td style="text-align:right;">
0.126
</td>
<td style="text-align:right;">
18.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.61
</td>
<td style="text-align:right;">
39.0
</td>
<td style="text-align:right;">
6.44–6.61
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.088
</td>
<td style="text-align:right;">
0.929
</td>
<td style="text-align:right;">
0.830
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
37.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
14.06
</td>
<td style="text-align:right;">
38.0
</td>
<td style="text-align:right;">
13.80–14.19
</td>
<td style="text-align:right;">
386
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.122
</td>
<td style="text-align:right;">
0.716
</td>
<td style="text-align:right;">
0.644
</td>
<td style="text-align:right;">
0.072
</td>
<td style="text-align:right;">
29.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.44
</td>
<td style="text-align:right;">
38.1
</td>
<td style="text-align:right;">
7.27–7.75
</td>
<td style="text-align:right;">
470
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.126
</td>
<td style="text-align:right;">
0.354
</td>
<td style="text-align:right;">
0.476
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.48
</td>
<td style="text-align:right;">
33.9
</td>
<td style="text-align:right;">
9.29–9.56
</td>
<td style="text-align:right;">
271
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.119
</td>
<td style="text-align:right;">
0.353
</td>
<td style="text-align:right;">
0.469
</td>
<td style="text-align:right;">
-0.116
</td>
<td style="text-align:right;">
26.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
6.27
</td>
<td style="text-align:right;">
31.4
</td>
<td style="text-align:right;">
6.24–6.27
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
0.354
</td>
<td style="text-align:right;">
0.491
</td>
<td style="text-align:right;">
-0.137
</td>
<td style="text-align:right;">
24.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.56
</td>
<td style="text-align:right;">
27.9
</td>
<td style="text-align:right;">
10.56–10.59
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.109
</td>
<td style="text-align:right;">
0.362
</td>
<td style="text-align:right;">
0.474
</td>
<td style="text-align:right;">
-0.112
</td>
<td style="text-align:right;">
19.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.79
</td>
<td style="text-align:right;">
21.4
</td>
<td style="text-align:right;">
11.59–11.85
</td>
<td style="text-align:right;">
265
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.095
</td>
<td style="text-align:right;">
0.381
</td>
<td style="text-align:right;">
0.450
</td>
<td style="text-align:right;">
-0.069
</td>
<td style="text-align:right;">
13.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
138.1
</td>
<td style="text-align:right;">
13.10–13.53
</td>
<td style="text-align:right;">
421
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.184
</td>
<td style="text-align:right;">
0.681
</td>
<td style="text-align:right;">
0.886
</td>
<td style="text-align:right;">
-0.204
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.10
</td>
<td style="text-align:right;">
43.5
</td>
<td style="text-align:right;">
12.07–12.10
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.090
</td>
<td style="text-align:right;">
0.815
</td>
<td style="text-align:right;">
0.893
</td>
<td style="text-align:right;">
-0.078
</td>
<td style="text-align:right;">
43.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
37.3
</td>
<td style="text-align:right;">
0.60–0.92
</td>
<td style="text-align:right;">
316
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.120
</td>
<td style="text-align:right;">
0.545
</td>
<td style="text-align:right;">
0.664
</td>
<td style="text-align:right;">
-0.119
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
15.6
</td>
<td style="text-align:right;">
0.29–0.70
</td>
<td style="text-align:right;">
410
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.082
</td>
<td style="text-align:right;">
0.466
</td>
<td style="text-align:right;">
0.581
</td>
<td style="text-align:right;">
-0.115
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.01
</td>
<td style="text-align:right;">
7.1
</td>
<td style="text-align:right;">
8.89–9.18
</td>
<td style="text-align:right;">
285
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.038
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.140
</td>
<td style="text-align:right;">
-0.024
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.98
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
5.87–6.07
</td>
<td style="text-align:right;">
191
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.039
</td>
<td style="text-align:right;">
0.137
</td>
<td style="text-align:right;">
0.184
</td>
<td style="text-align:right;">
-0.047
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.27
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
15.17–15.35
</td>
<td style="text-align:right;">
176
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.045
</td>
<td style="text-align:right;">
0.318
</td>
<td style="text-align:right;">
0.291
</td>
<td style="text-align:right;">
0.027
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.18
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
16.85–17.38
</td>
<td style="text-align:right;">
526
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.044
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
0.340
</td>
<td style="text-align:right;">
0.066
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
11.90
</td>
<td style="text-align:right;">
135.4
</td>
<td style="text-align:right;">
11.74–12.17
</td>
<td style="text-align:right;">
431
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.235
</td>
<td style="text-align:right;">
0.340
</td>
<td style="text-align:right;">
0.592
</td>
<td style="text-align:right;">
-0.252
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.51
</td>
<td style="text-align:right;">
132.9
</td>
<td style="text-align:right;">
10.28–10.70
</td>
<td style="text-align:right;">
416
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.230
</td>
<td style="text-align:right;">
0.427
</td>
<td style="text-align:right;">
0.668
</td>
<td style="text-align:right;">
-0.241
</td>
<td style="text-align:right;">
118.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.17
</td>
<td style="text-align:right;">
95.3
</td>
<td style="text-align:right;">
13.17–14.06
</td>
<td style="text-align:right;">
885
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.199
</td>
<td style="text-align:right;">
0.347
</td>
<td style="text-align:right;">
0.542
</td>
<td style="text-align:right;">
-0.195
</td>
<td style="text-align:right;">
95.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
9.28
</td>
<td style="text-align:right;">
74.3
</td>
<td style="text-align:right;">
9.23–9.28
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.171
</td>
<td style="text-align:right;">
0.491
</td>
<td style="text-align:right;">
0.665
</td>
<td style="text-align:right;">
-0.174
</td>
<td style="text-align:right;">
74.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
15.06
</td>
<td style="text-align:right;">
65.0
</td>
<td style="text-align:right;">
15.06–15.13
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
0.418
</td>
<td style="text-align:right;">
-0.163
</td>
<td style="text-align:right;">
65.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
15.52
</td>
<td style="text-align:right;">
14.7
</td>
<td style="text-align:right;">
15.43–15.66
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.070
</td>
<td style="text-align:right;">
0.784
</td>
<td style="text-align:right;">
0.692
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.00
</td>
<td style="text-align:right;">
10.0
</td>
<td style="text-align:right;">
12.92–13.07
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.043
</td>
<td style="text-align:right;">
0.914
</td>
<td style="text-align:right;">
0.833
</td>
<td style="text-align:right;">
0.080
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.92
</td>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
11.81–11.92
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.030
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.897
</td>
<td style="text-align:right;">
0.038
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
ht115 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
14.43
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
14.39–14.43
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.038
</td>
<td style="text-align:right;">
0.805
</td>
<td style="text-align:right;">
0.742
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
493.0
</td>
<td style="text-align:right;">
1.78–1.96
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.327
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:right;">
0.646
</td>
<td style="text-align:right;">
-0.332
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
223.4
</td>
<td style="text-align:right;">
0.23–0.78
</td>
<td style="text-align:right;">
553
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.227
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
0.613
</td>
<td style="text-align:right;">
-0.208
</td>
<td style="text-align:right;">
221.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.96
</td>
<td style="text-align:right;">
94.7
</td>
<td style="text-align:right;">
2.96–2.98
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.131
</td>
<td style="text-align:right;">
0.685
</td>
<td style="text-align:right;">
0.794
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
94.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
10.33
</td>
<td style="text-align:right;">
69.0
</td>
<td style="text-align:right;">
10.20–10.41
</td>
<td style="text-align:right;">
206
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.081
</td>
<td style="text-align:right;">
0.936
</td>
<td style="text-align:right;">
0.865
</td>
<td style="text-align:right;">
0.071
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.83
</td>
<td style="text-align:right;">
56.6
</td>
<td style="text-align:right;">
7.62–8.03
</td>
<td style="text-align:right;">
405
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.073
</td>
<td style="text-align:right;">
0.923
</td>
<td style="text-align:right;">
0.858
</td>
<td style="text-align:right;">
0.065
</td>
<td style="text-align:right;">
49.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.38
</td>
<td style="text-align:right;">
204.3
</td>
<td style="text-align:right;">
13.28–13.48
</td>
<td style="text-align:right;">
201
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.210
</td>
<td style="text-align:right;">
0.273
</td>
<td style="text-align:right;">
0.472
</td>
<td style="text-align:right;">
-0.199
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.49
</td>
<td style="text-align:right;">
176.9
</td>
<td style="text-align:right;">
11.42–11.63
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.197
</td>
<td style="text-align:right;">
0.273
</td>
<td style="text-align:right;">
0.476
</td>
<td style="text-align:right;">
-0.203
</td>
<td style="text-align:right;">
83.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.06
</td>
<td style="text-align:right;">
116.2
</td>
<td style="text-align:right;">
9.93–10.20
</td>
<td style="text-align:right;">
266
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.160
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:right;">
0.453
</td>
<td style="text-align:right;">
-0.163
</td>
<td style="text-align:right;">
105.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.69
</td>
<td style="text-align:right;">
106.2
</td>
<td style="text-align:right;">
4.55–4.76
</td>
<td style="text-align:right;">
206
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.157
</td>
<td style="text-align:right;">
0.386
</td>
<td style="text-align:right;">
0.537
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
57.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
6.70
</td>
<td style="text-align:right;">
101.7
</td>
<td style="text-align:right;">
6.63–6.76
</td>
<td style="text-align:right;">
131
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.148
</td>
<td style="text-align:right;">
0.287
</td>
<td style="text-align:right;">
0.441
</td>
<td style="text-align:right;">
-0.154
</td>
<td style="text-align:right;">
68.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.37
</td>
<td style="text-align:right;">
825.0
</td>
<td style="text-align:right;">
13.11–13.62
</td>
<td style="text-align:right;">
505
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.287
</td>
<td style="text-align:right;">
0.996
</td>
<td style="text-align:right;">
0.712
</td>
<td style="text-align:right;">
0.284
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.11
</td>
<td style="text-align:right;">
259.5
</td>
<td style="text-align:right;">
12.06–12.11
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.154
</td>
<td style="text-align:right;">
0.958
</td>
<td style="text-align:right;">
0.798
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
259.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
93.3
</td>
<td style="text-align:right;">
17.46–17.49
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.137
</td>
<td style="text-align:right;">
0.266
</td>
<td style="text-align:right;">
0.402
</td>
<td style="text-align:right;">
-0.137
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
47.2
</td>
<td style="text-align:right;">
0.02–0.07
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.106
</td>
<td style="text-align:right;">
0.589
</td>
<td style="text-align:right;">
0.489
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.64
</td>
<td style="text-align:right;">
42.0
</td>
<td style="text-align:right;">
15.53–15.91
</td>
<td style="text-align:right;">
371
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.092
</td>
<td style="text-align:right;">
0.246
</td>
<td style="text-align:right;">
0.348
</td>
<td style="text-align:right;">
-0.101
</td>
<td style="text-align:right;">
35.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
14.45
</td>
<td style="text-align:right;">
37.9
</td>
<td style="text-align:right;">
14.23–14.53
</td>
<td style="text-align:right;">
297
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
0.192
</td>
<td style="text-align:right;">
0.295
</td>
<td style="text-align:right;">
-0.103
</td>
<td style="text-align:right;">
23.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
3.24
</td>
<td style="text-align:right;">
14.1
</td>
<td style="text-align:right;">
2.96–3.34
</td>
<td style="text-align:right;">
375
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.056
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:right;">
0.385
</td>
<td style="text-align:right;">
-0.071
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.68
</td>
<td style="text-align:right;">
440.5
</td>
<td style="text-align:right;">
10.47–10.98
</td>
<td style="text-align:right;">
507
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.297
</td>
<td style="text-align:right;">
0.783
</td>
<td style="text-align:right;">
0.513
</td>
<td style="text-align:right;">
0.270
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.34
</td>
<td style="text-align:right;">
329.7
</td>
<td style="text-align:right;">
12.06–12.61
</td>
<td style="text-align:right;">
550
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.267
</td>
<td style="text-align:right;">
0.741
</td>
<td style="text-align:right;">
0.475
</td>
<td style="text-align:right;">
0.266
</td>
<td style="text-align:right;">
293.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.80
</td>
<td style="text-align:right;">
299.7
</td>
<td style="text-align:right;">
13.66–13.90
</td>
<td style="text-align:right;">
235
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.257
</td>
<td style="text-align:right;">
0.733
</td>
<td style="text-align:right;">
0.464
</td>
<td style="text-align:right;">
0.269
</td>
<td style="text-align:right;">
257.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
9.47
</td>
<td style="text-align:right;">
223.0
</td>
<td style="text-align:right;">
9.38–9.47
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.215
</td>
<td style="text-align:right;">
0.778
</td>
<td style="text-align:right;">
0.536
</td>
<td style="text-align:right;">
0.242
</td>
<td style="text-align:right;">
223.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.90
</td>
<td style="text-align:right;">
161.9
</td>
<td style="text-align:right;">
14.90–15.13
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.194
</td>
<td style="text-align:right;">
0.600
</td>
<td style="text-align:right;">
0.404
</td>
<td style="text-align:right;">
0.197
</td>
<td style="text-align:right;">
132.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.72
</td>
<td style="text-align:right;">
91.9
</td>
<td style="text-align:right;">
9.42–10.13
</td>
<td style="text-align:right;">
711
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
0.808
</td>
<td style="text-align:right;">
0.930
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.42
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:right;">
8.27–8.42
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
0.846
</td>
<td style="text-align:right;">
0.927
</td>
<td style="text-align:right;">
-0.081
</td>
<td style="text-align:right;">
71.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
16.06
</td>
<td style="text-align:right;">
46.5
</td>
<td style="text-align:right;">
15.93–16.17
</td>
<td style="text-align:right;">
231
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
0.812
</td>
<td style="text-align:right;">
0.740
</td>
<td style="text-align:right;">
0.072
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.27
</td>
<td style="text-align:right;">
36.2
</td>
<td style="text-align:right;">
7.25–7.27
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.053
</td>
<td style="text-align:right;">
0.875
</td>
<td style="text-align:right;">
0.939
</td>
<td style="text-align:right;">
-0.063
</td>
<td style="text-align:right;">
36.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.58
</td>
<td style="text-align:right;">
27.7
</td>
<td style="text-align:right;">
11.29–11.68
</td>
<td style="text-align:right;">
385
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.052
</td>
<td style="text-align:right;">
0.849
</td>
<td style="text-align:right;">
0.915
</td>
<td style="text-align:right;">
-0.066
</td>
<td style="text-align:right;">
25.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
671.7
</td>
<td style="text-align:right;">
1.79–1.96
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.375
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:right;">
0.715
</td>
<td style="text-align:right;">
-0.401
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
221.4
</td>
<td style="text-align:right;">
0.64–0.78
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.226
</td>
<td style="text-align:right;">
0.398
</td>
<td style="text-align:right;">
0.628
</td>
<td style="text-align:right;">
-0.230
</td>
<td style="text-align:right;">
221.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.96
</td>
<td style="text-align:right;">
137.1
</td>
<td style="text-align:right;">
2.96–2.99
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.155
</td>
<td style="text-align:right;">
0.684
</td>
<td style="text-align:right;">
0.827
</td>
<td style="text-align:right;">
-0.143
</td>
<td style="text-align:right;">
137.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.38
</td>
<td style="text-align:right;">
199.1
</td>
<td style="text-align:right;">
13.31–13.44
</td>
<td style="text-align:right;">
131
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.208
</td>
<td style="text-align:right;">
0.273
</td>
<td style="text-align:right;">
0.478
</td>
<td style="text-align:right;">
-0.205
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.81
</td>
<td style="text-align:right;">
91.5
</td>
<td style="text-align:right;">
11.45–11.97
</td>
<td style="text-align:right;">
510
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.143
</td>
<td style="text-align:right;">
0.307
</td>
<td style="text-align:right;">
0.433
</td>
<td style="text-align:right;">
-0.127
</td>
<td style="text-align:right;">
59.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
15.28
</td>
<td style="text-align:right;">
76.8
</td>
<td style="text-align:right;">
15.26–15.28
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.133
</td>
<td style="text-align:right;">
0.338
</td>
<td style="text-align:right;">
0.462
</td>
<td style="text-align:right;">
-0.123
</td>
<td style="text-align:right;">
49.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.05
</td>
<td style="text-align:right;">
55.6
</td>
<td style="text-align:right;">
9.96–10.13
</td>
<td style="text-align:right;">
170
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.110
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:right;">
0.397
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
42.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.54
</td>
<td style="text-align:right;">
37.0
</td>
<td style="text-align:right;">
3.48–3.65
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.095
</td>
<td style="text-align:right;">
0.419
</td>
<td style="text-align:right;">
0.539
</td>
<td style="text-align:right;">
-0.120
</td>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
477.8
</td>
<td style="text-align:right;">
13.10–13.54
</td>
<td style="text-align:right;">
435
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.192
</td>
<td style="text-align:right;">
0.997
</td>
<td style="text-align:right;">
0.794
</td>
<td style="text-align:right;">
0.203
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.09
</td>
<td style="text-align:right;">
130.9
</td>
<td style="text-align:right;">
12.05–12.09
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.959
</td>
<td style="text-align:right;">
0.859
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
130.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
48.6
</td>
<td style="text-align:right;">
0.21–0.64
</td>
<td style="text-align:right;">
421
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:right;">
0.681
</td>
<td style="text-align:right;">
0.599
</td>
<td style="text-align:right;">
0.082
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.14
</td>
<td style="text-align:right;">
196.2
</td>
<td style="text-align:right;">
6.82–7.32
</td>
<td style="text-align:right;">
495
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.170
</td>
<td style="text-align:right;">
0.281
</td>
<td style="text-align:right;">
0.107
</td>
<td style="text-align:right;">
0.175
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.32
</td>
<td style="text-align:right;">
178.9
</td>
<td style="text-align:right;">
8.32–8.86
</td>
<td style="text-align:right;">
535
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.258
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.141
</td>
<td style="text-align:right;">
178.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.86
</td>
<td style="text-align:right;">
132.5
</td>
<td style="text-align:right;">
9.86–9.98
</td>
<td style="text-align:right;">
125
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.135
</td>
<td style="text-align:right;">
0.260
</td>
<td style="text-align:right;">
0.113
</td>
<td style="text-align:right;">
0.146
</td>
<td style="text-align:right;">
132.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
11.12
</td>
<td style="text-align:right;">
123.3
</td>
<td style="text-align:right;">
10.98–11.33
</td>
<td style="text-align:right;">
341
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.133
</td>
<td style="text-align:right;">
0.241
</td>
<td style="text-align:right;">
0.122
</td>
<td style="text-align:right;">
0.119
</td>
<td style="text-align:right;">
112.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.82
</td>
<td style="text-align:right;">
76.9
</td>
<td style="text-align:right;">
5.71–5.82
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
0.226
</td>
<td style="text-align:right;">
0.114
</td>
<td style="text-align:right;">
0.112
</td>
<td style="text-align:right;">
76.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.48
</td>
<td style="text-align:right;">
310.7
</td>
<td style="text-align:right;">
12.25–12.63
</td>
<td style="text-align:right;">
371
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
0.729
</td>
<td style="text-align:right;">
0.460
</td>
<td style="text-align:right;">
0.269
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.78
</td>
<td style="text-align:right;">
296.5
</td>
<td style="text-align:right;">
10.40–11.25
</td>
<td style="text-align:right;">
846
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.244
</td>
<td style="text-align:right;">
0.783
</td>
<td style="text-align:right;">
0.532
</td>
<td style="text-align:right;">
0.250
</td>
<td style="text-align:right;">
259.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.78
</td>
<td style="text-align:right;">
259.4
</td>
<td style="text-align:right;">
13.63–13.89
</td>
<td style="text-align:right;">
255
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.240
</td>
<td style="text-align:right;">
0.723
</td>
<td style="text-align:right;">
0.487
</td>
<td style="text-align:right;">
0.237
</td>
<td style="text-align:right;">
230.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
15.09
</td>
<td style="text-align:right;">
171.7
</td>
<td style="text-align:right;">
14.89–15.36
</td>
<td style="text-align:right;">
469
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.201
</td>
<td style="text-align:right;">
0.573
</td>
<td style="text-align:right;">
0.366
</td>
<td style="text-align:right;">
0.207
</td>
<td style="text-align:right;">
158.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
9.40
</td>
<td style="text-align:right;">
147.1
</td>
<td style="text-align:right;">
9.35–9.40
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.175
</td>
<td style="text-align:right;">
0.768
</td>
<td style="text-align:right;">
0.583
</td>
<td style="text-align:right;">
0.185
</td>
<td style="text-align:right;">
147.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.98
</td>
<td style="text-align:right;">
118.4
</td>
<td style="text-align:right;">
9.83–10.14
</td>
<td style="text-align:right;">
311
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.112
</td>
<td style="text-align:right;">
0.818
</td>
<td style="text-align:right;">
0.929
</td>
<td style="text-align:right;">
-0.111
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.82
</td>
<td style="text-align:right;">
85.8
</td>
<td style="text-align:right;">
8.70–8.82
</td>
<td style="text-align:right;">
127
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.093
</td>
<td style="text-align:right;">
0.823
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:right;">
-0.118
</td>
<td style="text-align:right;">
85.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.69
</td>
<td style="text-align:right;">
56.2
</td>
<td style="text-align:right;">
7.55–7.69
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.071
</td>
<td style="text-align:right;">
0.855
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
-0.079
</td>
<td style="text-align:right;">
56.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
6.40
</td>
<td style="text-align:right;">
33.1
</td>
<td style="text-align:right;">
6.24–6.50
</td>
<td style="text-align:right;">
251
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.044
</td>
<td style="text-align:right;">
0.908
</td>
<td style="text-align:right;">
0.953
</td>
<td style="text-align:right;">
-0.044
</td>
<td style="text-align:right;">
23.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.14
</td>
<td style="text-align:right;">
27.9
</td>
<td style="text-align:right;">
11.14–11.19
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.052
</td>
<td style="text-align:right;">
0.860
</td>
<td style="text-align:right;">
0.912
</td>
<td style="text-align:right;">
-0.052
</td>
<td style="text-align:right;">
27.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
197.7
</td>
<td style="text-align:right;">
1.77–1.94
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.280
</td>
<td style="text-align:right;">
0.314
</td>
<td style="text-align:right;">
0.597
</td>
<td style="text-align:right;">
-0.283
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
50.1
</td>
<td style="text-align:right;">
0.62–0.77
</td>
<td style="text-align:right;">
150
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.145
</td>
<td style="text-align:right;">
0.394
</td>
<td style="text-align:right;">
0.543
</td>
<td style="text-align:right;">
-0.149
</td>
<td style="text-align:right;">
50.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.85
</td>
<td style="text-align:right;">
36.1
</td>
<td style="text-align:right;">
7.62–8.04
</td>
<td style="text-align:right;">
416
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
0.926
</td>
<td style="text-align:right;">
0.839
</td>
<td style="text-align:right;">
0.088
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
5.92
</td>
<td style="text-align:right;">
35.2
</td>
<td style="text-align:right;">
5.73–6.15
</td>
<td style="text-align:right;">
416
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
0.917
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
27.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
10.24
</td>
<td style="text-align:right;">
29.7
</td>
<td style="text-align:right;">
9.97–10.36
</td>
<td style="text-align:right;">
385
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.074
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.881
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
22.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.37
</td>
<td style="text-align:right;">
90.2
</td>
<td style="text-align:right;">
13.26–13.44
</td>
<td style="text-align:right;">
177
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.189
</td>
<td style="text-align:right;">
0.273
</td>
<td style="text-align:right;">
0.434
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.52
</td>
<td style="text-align:right;">
90.2
</td>
<td style="text-align:right;">
10.40–10.65
</td>
<td style="text-align:right;">
251
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.192
</td>
<td style="text-align:right;">
0.282
</td>
<td style="text-align:right;">
0.472
</td>
<td style="text-align:right;">
-0.190
</td>
<td style="text-align:right;">
26.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.65
</td>
<td style="text-align:right;">
80.7
</td>
<td style="text-align:right;">
11.65–11.76
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.181
</td>
<td style="text-align:right;">
0.291
</td>
<td style="text-align:right;">
0.481
</td>
<td style="text-align:right;">
-0.190
</td>
<td style="text-align:right;">
68.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
6.40
</td>
<td style="text-align:right;">
64.3
</td>
<td style="text-align:right;">
6.28–6.58
</td>
<td style="text-align:right;">
301
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.162
</td>
<td style="text-align:right;">
0.306
</td>
<td style="text-align:right;">
0.447
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
37.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.91
</td>
<td style="text-align:right;">
53.0
</td>
<td style="text-align:right;">
4.68–5.05
</td>
<td style="text-align:right;">
367
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.149
</td>
<td style="text-align:right;">
0.361
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
30.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.29
</td>
<td style="text-align:right;">
121.5
</td>
<td style="text-align:right;">
13.05–13.56
</td>
<td style="text-align:right;">
506
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.131
</td>
<td style="text-align:right;">
0.996
</td>
<td style="text-align:right;">
0.891
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.05
</td>
<td style="text-align:right;">
32.5
</td>
<td style="text-align:right;">
11.87–12.05
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.064
</td>
<td style="text-align:right;">
0.960
</td>
<td style="text-align:right;">
0.896
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:right;">
32.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.14
</td>
<td style="text-align:right;">
68.4
</td>
<td style="text-align:right;">
6.97–7.27
</td>
<td style="text-align:right;">
301
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.132
</td>
<td style="text-align:right;">
0.280
</td>
<td style="text-align:right;">
0.146
</td>
<td style="text-align:right;">
0.133
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.27
</td>
<td style="text-align:right;">
51.8
</td>
<td style="text-align:right;">
8.27–8.87
</td>
<td style="text-align:right;">
594
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.240
</td>
<td style="text-align:right;">
0.145
</td>
<td style="text-align:right;">
0.096
</td>
<td style="text-align:right;">
51.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.90
</td>
<td style="text-align:right;">
48.1
</td>
<td style="text-align:right;">
10.73–11.15
</td>
<td style="text-align:right;">
421
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.107
</td>
<td style="text-align:right;">
0.251
</td>
<td style="text-align:right;">
0.141
</td>
<td style="text-align:right;">
0.109
</td>
<td style="text-align:right;">
41.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.27
</td>
<td style="text-align:right;">
28.1
</td>
<td style="text-align:right;">
12.15–12.55
</td>
<td style="text-align:right;">
397
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.085
</td>
<td style="text-align:right;">
0.238
</td>
<td style="text-align:right;">
0.139
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
26.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
25.7
</td>
<td style="text-align:right;">
17.46–17.49
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.096
</td>
<td style="text-align:right;">
0.266
</td>
<td style="text-align:right;">
0.354
</td>
<td style="text-align:right;">
-0.089
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.81
</td>
<td style="text-align:right;">
81.9
</td>
<td style="text-align:right;">
10.68–10.99
</td>
<td style="text-align:right;">
306
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.174
</td>
<td style="text-align:right;">
0.780
</td>
<td style="text-align:right;">
0.616
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.53
</td>
<td style="text-align:right;">
70.1
</td>
<td style="text-align:right;">
12.33–12.71
</td>
<td style="text-align:right;">
375
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.168
</td>
<td style="text-align:right;">
0.728
</td>
<td style="text-align:right;">
0.579
</td>
<td style="text-align:right;">
0.149
</td>
<td style="text-align:right;">
53.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
16.24
</td>
<td style="text-align:right;">
69.3
</td>
<td style="text-align:right;">
16.06–16.42
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.170
</td>
<td style="text-align:right;">
0.612
</td>
<td style="text-align:right;">
0.417
</td>
<td style="text-align:right;">
0.195
</td>
<td style="text-align:right;">
41.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.81
</td>
<td style="text-align:right;">
67.0
</td>
<td style="text-align:right;">
13.71–13.89
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:right;">
0.732
</td>
<td style="text-align:right;">
0.558
</td>
<td style="text-align:right;">
0.174
</td>
<td style="text-align:right;">
47.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
20.49
</td>
<td style="text-align:right;">
65.0
</td>
<td style="text-align:right;">
20.45–20.49
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.031
</td>
<td style="text-align:right;">
0.013
</td>
<td style="text-align:right;">
0.013
</td>
<td style="text-align:right;">
-0.000
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.77
</td>
<td style="text-align:right;">
92.6
</td>
<td style="text-align:right;">
9.61–10.06
</td>
<td style="text-align:right;">
440
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.118
</td>
<td style="text-align:right;">
0.818
</td>
<td style="text-align:right;">
0.942
</td>
<td style="text-align:right;">
-0.124
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.61
</td>
<td style="text-align:right;">
65.4
</td>
<td style="text-align:right;">
8.40–8.61
</td>
<td style="text-align:right;">
217
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.096
</td>
<td style="text-align:right;">
0.839
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
-0.106
</td>
<td style="text-align:right;">
63.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.39
</td>
<td style="text-align:right;">
32.2
</td>
<td style="text-align:right;">
7.37–7.39
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.063
</td>
<td style="text-align:right;">
0.887
</td>
<td style="text-align:right;">
0.936
</td>
<td style="text-align:right;">
-0.049
</td>
<td style="text-align:right;">
32.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
15.88
</td>
<td style="text-align:right;">
25.0
</td>
<td style="text-align:right;">
15.69–16.10
</td>
<td style="text-align:right;">
405
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
0.798
</td>
<td style="text-align:right;">
0.724
</td>
<td style="text-align:right;">
0.074
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
mig6 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.06
</td>
<td style="text-align:right;">
17.0
</td>
<td style="text-align:right;">
11.06–11.09
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.052
</td>
<td style="text-align:right;">
0.864
</td>
<td style="text-align:right;">
0.915
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
17.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
8.58
</td>
<td style="text-align:right;">
37.7
</td>
<td style="text-align:right;">
8.42–8.96
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.063
</td>
<td style="text-align:right;">
0.849
</td>
<td style="text-align:right;">
0.904
</td>
<td style="text-align:right;">
-0.055
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.01
</td>
<td style="text-align:right;">
24.8
</td>
<td style="text-align:right;">
6.64–7.35
</td>
<td style="text-align:right;">
712
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.050
</td>
<td style="text-align:right;">
0.859
</td>
<td style="text-align:right;">
0.904
</td>
<td style="text-align:right;">
-0.045
</td>
<td style="text-align:right;">
20.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
10.36
</td>
<td style="text-align:right;">
22.7
</td>
<td style="text-align:right;">
10.27–10.62
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.050
</td>
<td style="text-align:right;">
0.860
</td>
<td style="text-align:right;">
0.905
</td>
<td style="text-align:right;">
-0.045
</td>
<td style="text-align:right;">
20.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
5.56
</td>
<td style="text-align:right;">
18.2
</td>
<td style="text-align:right;">
5.42–5.63
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.044
</td>
<td style="text-align:right;">
0.858
</td>
<td style="text-align:right;">
0.909
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
15.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.65
</td>
<td style="text-align:right;">
17.2
</td>
<td style="text-align:right;">
11.62–11.82
</td>
<td style="text-align:right;">
193
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.049
</td>
<td style="text-align:right;">
0.829
</td>
<td style="text-align:right;">
0.869
</td>
<td style="text-align:right;">
-0.040
</td>
<td style="text-align:right;">
10.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.18
</td>
<td style="text-align:right;">
39.2
</td>
<td style="text-align:right;">
7.15–7.23
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.095
</td>
<td style="text-align:right;">
0.441
</td>
<td style="text-align:right;">
0.320
</td>
<td style="text-align:right;">
0.121
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.06
</td>
<td style="text-align:right;">
29.3
</td>
<td style="text-align:right;">
8.93–9.14
</td>
<td style="text-align:right;">
206
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
0.459
</td>
<td style="text-align:right;">
0.355
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
13.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.54
</td>
<td style="text-align:right;">
24.1
</td>
<td style="text-align:right;">
4.42–4.61
</td>
<td style="text-align:right;">
185
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.077
</td>
<td style="text-align:right;">
0.501
</td>
<td style="text-align:right;">
0.426
</td>
<td style="text-align:right;">
0.075
</td>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
5.88
</td>
<td style="text-align:right;">
23.1
</td>
<td style="text-align:right;">
5.83–5.92
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.075
</td>
<td style="text-align:right;">
0.507
</td>
<td style="text-align:right;">
0.443
</td>
<td style="text-align:right;">
0.064
</td>
<td style="text-align:right;">
7.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.52
</td>
<td style="text-align:right;">
17.7
</td>
<td style="text-align:right;">
10.45–10.62
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.066
</td>
<td style="text-align:right;">
0.444
</td>
<td style="text-align:right;">
0.382
</td>
<td style="text-align:right;">
0.062
</td>
<td style="text-align:right;">
9.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
8.01
</td>
<td style="text-align:right;">
80.1
</td>
<td style="text-align:right;">
7.95–8.06
</td>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.049
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.989
</td>
<td style="text-align:right;">
-0.054
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
6.25
</td>
<td style="text-align:right;">
62.4
</td>
<td style="text-align:right;">
6.05–6.39
</td>
<td style="text-align:right;">
330
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.047
</td>
<td style="text-align:right;">
0.936
</td>
<td style="text-align:right;">
0.983
</td>
<td style="text-align:right;">
-0.047
</td>
<td style="text-align:right;">
22.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.78
</td>
<td style="text-align:right;">
61.3
</td>
<td style="text-align:right;">
13.69–13.78
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.101
</td>
<td style="text-align:right;">
0.801
</td>
<td style="text-align:right;">
0.882
</td>
<td style="text-align:right;">
-0.081
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.68
</td>
<td style="text-align:right;">
38.7
</td>
<td style="text-align:right;">
12.61–12.68
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.079
</td>
<td style="text-align:right;">
0.761
</td>
<td style="text-align:right;">
0.842
</td>
<td style="text-align:right;">
-0.081
</td>
<td style="text-align:right;">
38.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
4.78
</td>
<td style="text-align:right;">
22.8
</td>
<td style="text-align:right;">
4.65–4.89
</td>
<td style="text-align:right;">
235
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.034
</td>
<td style="text-align:right;">
0.924
</td>
<td style="text-align:right;">
0.960
</td>
<td style="text-align:right;">
-0.037
</td>
<td style="text-align:right;">
21.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.87
</td>
<td style="text-align:right;">
151.5
</td>
<td style="text-align:right;">
8.37–9.67
</td>
<td style="text-align:right;">
1291
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.145
</td>
<td style="text-align:right;">
0.232
</td>
<td style="text-align:right;">
0.109
</td>
<td style="text-align:right;">
0.123
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.10
</td>
<td style="text-align:right;">
149.3
</td>
<td style="text-align:right;">
6.58–7.37
</td>
<td style="text-align:right;">
791
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.146
</td>
<td style="text-align:right;">
0.256
</td>
<td style="text-align:right;">
0.106
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:right;">
135.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.67
</td>
<td style="text-align:right;">
124.0
</td>
<td style="text-align:right;">
10.67–10.96
</td>
<td style="text-align:right;">
286
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.132
</td>
<td style="text-align:right;">
0.239
</td>
<td style="text-align:right;">
0.112
</td>
<td style="text-align:right;">
0.127
</td>
<td style="text-align:right;">
124.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.57
</td>
<td style="text-align:right;">
108.3
</td>
<td style="text-align:right;">
5.09–5.57
</td>
<td style="text-align:right;">
485
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.126
</td>
<td style="text-align:right;">
0.253
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.137
</td>
<td style="text-align:right;">
108.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
11.96
</td>
<td style="text-align:right;">
92.3
</td>
<td style="text-align:right;">
11.96–12.09
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.117
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
0.120
</td>
<td style="text-align:right;">
0.135
</td>
<td style="text-align:right;">
92.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
86.5
</td>
<td style="text-align:right;">
1.20–1.48
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.132
</td>
<td style="text-align:right;">
0.625
</td>
<td style="text-align:right;">
0.742
</td>
<td style="text-align:right;">
-0.116
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
67.9
</td>
<td style="text-align:right;">
0.01–0.20
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.117
</td>
<td style="text-align:right;">
0.642
</td>
<td style="text-align:right;">
0.756
</td>
<td style="text-align:right;">
-0.114
</td>
<td style="text-align:right;">
60.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.11
</td>
<td style="text-align:right;">
50.6
</td>
<td style="text-align:right;">
2.93–3.31
</td>
<td style="text-align:right;">
376
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.099
</td>
<td style="text-align:right;">
0.656
</td>
<td style="text-align:right;">
0.754
</td>
<td style="text-align:right;">
-0.098
</td>
<td style="text-align:right;">
38.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
4.35
</td>
<td style="text-align:right;">
46.9
</td>
<td style="text-align:right;">
4.31–4.54
</td>
<td style="text-align:right;">
220
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.097
</td>
<td style="text-align:right;">
0.635
</td>
<td style="text-align:right;">
0.731
</td>
<td style="text-align:right;">
-0.095
</td>
<td style="text-align:right;">
45.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
6.13
</td>
<td style="text-align:right;">
24.3
</td>
<td style="text-align:right;">
5.97–6.21
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.073
</td>
<td style="text-align:right;">
0.653
</td>
<td style="text-align:right;">
0.721
</td>
<td style="text-align:right;">
-0.068
</td>
<td style="text-align:right;">
19.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
4.51
</td>
<td style="text-align:right;">
17.6
</td>
<td style="text-align:right;">
4.45–4.56
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.027
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
0.975
</td>
<td style="text-align:right;">
-0.030
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
2.47
</td>
<td style="text-align:right;">
12.4
</td>
<td style="text-align:right;">
2.38–2.53
</td>
<td style="text-align:right;">
152
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.043
</td>
<td style="text-align:right;">
0.846
</td>
<td style="text-align:right;">
0.809
</td>
<td style="text-align:right;">
0.037
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
16.01
</td>
<td style="text-align:right;">
9.8
</td>
<td style="text-align:right;">
15.91–16.10
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.044
</td>
<td style="text-align:right;">
0.743
</td>
<td style="text-align:right;">
0.748
</td>
<td style="text-align:right;">
-0.005
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
5.56
</td>
<td style="text-align:right;">
9.4
</td>
<td style="text-align:right;">
5.56–5.98
</td>
<td style="text-align:right;">
415
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.020
</td>
<td style="text-align:right;">
0.946
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
-0.026
</td>
<td style="text-align:right;">
8.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
8.2
</td>
<td style="text-align:right;">
0.57–0.75
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.037
</td>
<td style="text-align:right;">
0.807
</td>
<td style="text-align:right;">
0.795
</td>
<td style="text-align:right;">
0.012
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
21.1
</td>
<td style="text-align:right;">
0.01–0.18
</td>
<td style="text-align:right;">
170
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.095
</td>
<td style="text-align:right;">
0.632
</td>
<td style="text-align:right;">
0.532
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:right;">
13.2
</td>
<td style="text-align:right;">
2.31–2.51
</td>
<td style="text-align:right;">
195
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.070
</td>
<td style="text-align:right;">
0.781
</td>
<td style="text-align:right;">
0.671
</td>
<td style="text-align:right;">
0.110
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.18
</td>
<td style="text-align:right;">
11.3
</td>
<td style="text-align:right;">
1.18–1.21
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.069
</td>
<td style="text-align:right;">
0.665
</td>
<td style="text-align:right;">
0.609
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
11.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
3.52
</td>
<td style="text-align:right;">
8.1
</td>
<td style="text-align:right;">
3.52–3.54
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.045
</td>
<td style="text-align:right;">
0.876
</td>
<td style="text-align:right;">
0.832
</td>
<td style="text-align:right;">
0.044
</td>
<td style="text-align:right;">
8.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
12.52
</td>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
12.47–12.57
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.047
</td>
<td style="text-align:right;">
0.758
</td>
<td style="text-align:right;">
0.699
</td>
<td style="text-align:right;">
0.059
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:right;">
13.1
</td>
<td style="text-align:right;">
3.28–3.39
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.074
</td>
<td style="text-align:right;">
0.680
</td>
<td style="text-align:right;">
0.609
</td>
<td style="text-align:right;">
0.071
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
1.92
</td>
<td style="text-align:right;">
8.4
</td>
<td style="text-align:right;">
1.90–1.94
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.920
</td>
<td style="text-align:right;">
0.018
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
14.23
</td>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
14.11–14.27
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.053
</td>
<td style="text-align:right;">
0.446
</td>
<td style="text-align:right;">
0.496
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.50
</td>
<td style="text-align:right;">
103.8
</td>
<td style="text-align:right;">
13.20–13.77
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.157
</td>
<td style="text-align:right;">
0.702
</td>
<td style="text-align:right;">
0.866
</td>
<td style="text-align:right;">
-0.165
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.20
</td>
<td style="text-align:right;">
46.4
</td>
<td style="text-align:right;">
12.15–12.20
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.096
</td>
<td style="text-align:right;">
0.804
</td>
<td style="text-align:right;">
0.884
</td>
<td style="text-align:right;">
-0.081
</td>
<td style="text-align:right;">
46.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
1.11
</td>
<td style="text-align:right;">
14.3
</td>
<td style="text-align:right;">
0.95–1.22
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.070
</td>
<td style="text-align:right;">
0.648
</td>
<td style="text-align:right;">
0.737
</td>
<td style="text-align:right;">
-0.090
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.43
</td>
<td style="text-align:right;">
50.8
</td>
<td style="text-align:right;">
10.13–10.77
</td>
<td style="text-align:right;">
637
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.110
</td>
<td style="text-align:right;">
0.230
</td>
<td style="text-align:right;">
0.139
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
7.14
</td>
<td style="text-align:right;">
46.3
</td>
<td style="text-align:right;">
6.95–7.32
</td>
<td style="text-align:right;">
368
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.108
</td>
<td style="text-align:right;">
0.260
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:right;">
0.111
</td>
<td style="text-align:right;">
35.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.12
</td>
<td style="text-align:right;">
41.7
</td>
<td style="text-align:right;">
8.59–9.12
</td>
<td style="text-align:right;">
535
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
0.240
</td>
<td style="text-align:right;">
0.131
</td>
<td style="text-align:right;">
0.109
</td>
<td style="text-align:right;">
41.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.43
</td>
<td style="text-align:right;">
36.7
</td>
<td style="text-align:right;">
11.95–12.64
</td>
<td style="text-align:right;">
683
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.283
</td>
<td style="text-align:right;">
0.164
</td>
<td style="text-align:right;">
0.119
</td>
<td style="text-align:right;">
33.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.65
</td>
<td style="text-align:right;">
26.3
</td>
<td style="text-align:right;">
4.05–5.07
</td>
<td style="text-align:right;">
1020
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.087
</td>
<td style="text-align:right;">
0.247
</td>
<td style="text-align:right;">
0.191
</td>
<td style="text-align:right;">
0.056
</td>
<td style="text-align:right;">
17.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.42
</td>
<td style="text-align:right;">
45.9
</td>
<td style="text-align:right;">
10.23–10.61
</td>
<td style="text-align:right;">
371
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.136
</td>
<td style="text-align:right;">
0.519
</td>
<td style="text-align:right;">
0.650
</td>
<td style="text-align:right;">
-0.131
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
6.46
</td>
<td style="text-align:right;">
33.1
</td>
<td style="text-align:right;">
6.29–6.61
</td>
<td style="text-align:right;">
316
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
0.630
</td>
<td style="text-align:right;">
0.730
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
13.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
8.30
</td>
<td style="text-align:right;">
29.5
</td>
<td style="text-align:right;">
8.08–8.52
</td>
<td style="text-align:right;">
432
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.106
</td>
<td style="text-align:right;">
0.567
</td>
<td style="text-align:right;">
0.661
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:right;">
13.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.04
</td>
<td style="text-align:right;">
28.8
</td>
<td style="text-align:right;">
11.84–12.23
</td>
<td style="text-align:right;">
388
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.110
</td>
<td style="text-align:right;">
0.459
</td>
<td style="text-align:right;">
0.563
</td>
<td style="text-align:right;">
-0.103
</td>
<td style="text-align:right;">
24.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
4.57
</td>
<td style="text-align:right;">
26.9
</td>
<td style="text-align:right;">
4.44–4.74
</td>
<td style="text-align:right;">
296
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.098
</td>
<td style="text-align:right;">
0.642
</td>
<td style="text-align:right;">
0.736
</td>
<td style="text-align:right;">
-0.094
</td>
<td style="text-align:right;">
21.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
par1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
4.78
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
4.72–4.83
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.017
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:right;">
0.956
</td>
<td style="text-align:right;">
-0.016
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
1.88
</td>
<td style="text-align:right;">
751.1
</td>
<td style="text-align:right;">
1.80–1.96
</td>
<td style="text-align:right;">
162
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.391
</td>
<td style="text-align:right;">
0.725
</td>
<td style="text-align:right;">
0.315
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.96
</td>
<td style="text-align:right;">
199.4
</td>
<td style="text-align:right;">
2.96–2.98
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.181
</td>
<td style="text-align:right;">
0.855
</td>
<td style="text-align:right;">
0.684
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
199.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
195.8
</td>
<td style="text-align:right;">
0.62–0.79
</td>
<td style="text-align:right;">
169
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.213
</td>
<td style="text-align:right;">
0.609
</td>
<td style="text-align:right;">
0.403
</td>
<td style="text-align:right;">
0.206
</td>
<td style="text-align:right;">
195.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.23
</td>
<td style="text-align:right;">
166.8
</td>
<td style="text-align:right;">
13.05–13.38
</td>
<td style="text-align:right;">
326
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.176
</td>
<td style="text-align:right;">
0.799
</td>
<td style="text-align:right;">
0.596
</td>
<td style="text-align:right;">
0.203
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
14.38
</td>
<td style="text-align:right;">
144.0
</td>
<td style="text-align:right;">
14.38–14.42
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.174
</td>
<td style="text-align:right;">
0.765
</td>
<td style="text-align:right;">
0.597
</td>
<td style="text-align:right;">
0.168
</td>
<td style="text-align:right;">
108.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
43.7
</td>
<td style="text-align:right;">
0.79–0.98
</td>
<td style="text-align:right;">
186
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.099
</td>
<td style="text-align:right;">
0.599
</td>
<td style="text-align:right;">
0.707
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
13.62
</td>
<td style="text-align:right;">
36.6
</td>
<td style="text-align:right;">
13.56–13.69
</td>
<td style="text-align:right;">
121
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.090
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
0.285
</td>
<td style="text-align:right;">
0.121
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
11.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.25
</td>
<td style="text-align:right;">
28.0
</td>
<td style="text-align:right;">
9.17–9.34
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.074
</td>
<td style="text-align:right;">
0.254
</td>
<td style="text-align:right;">
0.309
</td>
<td style="text-align:right;">
-0.055
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
10.7
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.18
</td>
<td style="text-align:right;">
27.7
</td>
<td style="text-align:right;">
7.07–7.23
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.073
</td>
<td style="text-align:right;">
0.251
</td>
<td style="text-align:right;">
0.287
</td>
<td style="text-align:right;">
-0.036
</td>
<td style="text-align:right;">
14.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
9.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
3.09
</td>
<td style="text-align:right;">
27.5
</td>
<td style="text-align:right;">
3.01–3.16
</td>
<td style="text-align:right;">
145
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.082
</td>
<td style="text-align:right;">
0.420
</td>
<td style="text-align:right;">
0.521
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
165.3
</td>
<td style="text-align:right;">
0.28–0.59
</td>
<td style="text-align:right;">
301
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.192
</td>
<td style="text-align:right;">
0.513
</td>
<td style="text-align:right;">
0.675
</td>
<td style="text-align:right;">
-0.162
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
67.3
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
1.86
</td>
<td style="text-align:right;">
74.8
</td>
<td style="text-align:right;">
1.74–2.00
</td>
<td style="text-align:right;">
260
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
0.573
</td>
<td style="text-align:right;">
0.698
</td>
<td style="text-align:right;">
-0.126
</td>
<td style="text-align:right;">
22.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
17.0
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.40
</td>
<td style="text-align:right;">
51.3
</td>
<td style="text-align:right;">
13.18–13.58
</td>
<td style="text-align:right;">
400
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.026
</td>
<td style="text-align:right;">
0.981
</td>
<td style="text-align:right;">
0.997
</td>
<td style="text-align:right;">
-0.016
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
31.0
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
9.35
</td>
<td style="text-align:right;">
37.6
</td>
<td style="text-align:right;">
9.26–9.53
</td>
<td style="text-align:right;">
265
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.045
</td>
<td style="text-align:right;">
0.961
</td>
<td style="text-align:right;">
0.904
</td>
<td style="text-align:right;">
0.057
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
31.7
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
10.53
</td>
<td style="text-align:right;">
22.9
</td>
<td style="text-align:right;">
10.53–10.55
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.032
</td>
<td style="text-align:right;">
0.960
</td>
<td style="text-align:right;">
0.914
</td>
<td style="text-align:right;">
0.046
</td>
<td style="text-align:right;">
22.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
11.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
167.3
</td>
<td style="text-align:right;">
17.44–17.49
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.185
</td>
<td style="text-align:right;">
0.420
</td>
<td style="text-align:right;">
0.266
</td>
<td style="text-align:right;">
0.154
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.88
</td>
<td style="text-align:right;">
135.5
</td>
<td style="text-align:right;">
15.70–15.98
</td>
<td style="text-align:right;">
275
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
0.402
</td>
<td style="text-align:right;">
0.215
</td>
<td style="text-align:right;">
0.188
</td>
<td style="text-align:right;">
104.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.89
</td>
<td style="text-align:right;">
91.0
</td>
<td style="text-align:right;">
6.73–7.17
</td>
<td style="text-align:right;">
437
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.123
</td>
<td style="text-align:right;">
0.150
</td>
<td style="text-align:right;">
0.267
</td>
<td style="text-align:right;">
-0.117
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
14.41
</td>
<td style="text-align:right;">
78.0
</td>
<td style="text-align:right;">
14.18–14.57
</td>
<td style="text-align:right;">
381
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.122
</td>
<td style="text-align:right;">
0.357
</td>
<td style="text-align:right;">
0.193
</td>
<td style="text-align:right;">
0.163
</td>
<td style="text-align:right;">
60.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.17
</td>
<td style="text-align:right;">
72.9
</td>
<td style="text-align:right;">
8.17–8.43
</td>
<td style="text-align:right;">
259
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
0.158
</td>
<td style="text-align:right;">
0.276
</td>
<td style="text-align:right;">
-0.118
</td>
<td style="text-align:right;">
72.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.81
</td>
<td style="text-align:right;">
730.3
</td>
<td style="text-align:right;">
13.63–13.91
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.388
</td>
<td style="text-align:right;">
0.336
</td>
<td style="text-align:right;">
0.732
</td>
<td style="text-align:right;">
-0.395
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.55
</td>
<td style="text-align:right;">
565.4
</td>
<td style="text-align:right;">
12.15–12.63
</td>
<td style="text-align:right;">
481
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.346
</td>
<td style="text-align:right;">
0.377
</td>
<td style="text-align:right;">
0.718
</td>
<td style="text-align:right;">
-0.342
</td>
<td style="text-align:right;">
543.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.91
</td>
<td style="text-align:right;">
526.4
</td>
<td style="text-align:right;">
14.91–15.28
</td>
<td style="text-align:right;">
370
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.333
</td>
<td style="text-align:right;">
0.242
</td>
<td style="text-align:right;">
0.599
</td>
<td style="text-align:right;">
-0.356
</td>
<td style="text-align:right;">
447.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
16.33
</td>
<td style="text-align:right;">
513.6
</td>
<td style="text-align:right;">
16.28–16.46
</td>
<td style="text-align:right;">
182
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.330
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:right;">
0.620
</td>
<td style="text-align:right;">
-0.330
</td>
<td style="text-align:right;">
466.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.75
</td>
<td style="text-align:right;">
355.4
</td>
<td style="text-align:right;">
10.51–10.89
</td>
<td style="text-align:right;">
372
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.266
</td>
<td style="text-align:right;">
0.551
</td>
<td style="text-align:right;">
0.787
</td>
<td style="text-align:right;">
-0.235
</td>
<td style="text-align:right;">
322.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
7.0
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.58
</td>
<td style="text-align:right;">
84.6
</td>
<td style="text-align:right;">
13.36–13.69
</td>
<td style="text-align:right;">
320
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
0.784
</td>
<td style="text-align:right;">
0.877
</td>
<td style="text-align:right;">
-0.093
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
166.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
63.8
</td>
<td style="text-align:right;">
2.76–3.13
</td>
<td style="text-align:right;">
366
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.086
</td>
<td style="text-align:right;">
0.901
</td>
<td style="text-align:right;">
0.802
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
30.1
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
16.35
</td>
<td style="text-align:right;">
45.4
</td>
<td style="text-align:right;">
16.25–16.41
</td>
<td style="text-align:right;">
156
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.092
</td>
<td style="text-align:right;">
0.689
</td>
<td style="text-align:right;">
0.772
</td>
<td style="text-align:right;">
-0.084
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
106.3
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.39
</td>
<td style="text-align:right;">
41.5
</td>
<td style="text-align:right;">
11.20–11.48
</td>
<td style="text-align:right;">
281
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.077
</td>
<td style="text-align:right;">
0.782
</td>
<td style="text-align:right;">
0.863
</td>
<td style="text-align:right;">
-0.081
</td>
<td style="text-align:right;">
18.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
183.2
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs mig6
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
4.14
</td>
<td style="text-align:right;">
30.1
</td>
<td style="text-align:right;">
4.14–4.16
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
0.958
</td>
<td style="text-align:right;">
0.930
</td>
<td style="text-align:right;">
0.029
</td>
<td style="text-align:right;">
30.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
83.7
</td>
<td style="text-align:right;">
FALSE
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.54
</td>
<td style="text-align:right;">
135.7
</td>
<td style="text-align:right;">
13.21–13.86
</td>
<td style="text-align:right;">
645
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
0.646
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
10.89
</td>
<td style="text-align:right;">
102.1
</td>
<td style="text-align:right;">
10.56–11.24
</td>
<td style="text-align:right;">
681
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
0.918
</td>
<td style="text-align:right;">
0.825
</td>
<td style="text-align:right;">
0.093
</td>
<td style="text-align:right;">
40.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.56
</td>
<td style="text-align:right;">
87.1
</td>
<td style="text-align:right;">
9.43–9.56
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.838
</td>
<td style="text-align:right;">
0.097
</td>
<td style="text-align:right;">
83.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.59
</td>
<td style="text-align:right;">
82.2
</td>
<td style="text-align:right;">
7.32–7.95
</td>
<td style="text-align:right;">
632
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.085
</td>
<td style="text-align:right;">
0.943
</td>
<td style="text-align:right;">
0.854
</td>
<td style="text-align:right;">
0.089
</td>
<td style="text-align:right;">
65.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
14.86
</td>
<td style="text-align:right;">
58.6
</td>
<td style="text-align:right;">
14.86–14.91
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.663
</td>
<td style="text-align:right;">
0.545
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:right;">
58.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.76
</td>
<td style="text-align:right;">
226.0
</td>
<td style="text-align:right;">
10.70–10.81
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.221
</td>
<td style="text-align:right;">
0.251
</td>
<td style="text-align:right;">
0.475
</td>
<td style="text-align:right;">
-0.224
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.20
</td>
<td style="text-align:right;">
198.4
</td>
<td style="text-align:right;">
7.12–7.29
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.203
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
0.443
</td>
<td style="text-align:right;">
-0.188
</td>
<td style="text-align:right;">
132.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.22
</td>
<td style="text-align:right;">
197.6
</td>
<td style="text-align:right;">
8.97–9.49
</td>
<td style="text-align:right;">
515
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.203
</td>
<td style="text-align:right;">
0.251
</td>
<td style="text-align:right;">
0.435
</td>
<td style="text-align:right;">
-0.184
</td>
<td style="text-align:right;">
132.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.81
</td>
<td style="text-align:right;">
176.7
</td>
<td style="text-align:right;">
11.81–11.90
</td>
<td style="text-align:right;">
90
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.196
</td>
<td style="text-align:right;">
0.268
</td>
<td style="text-align:right;">
0.467
</td>
<td style="text-align:right;">
-0.199
</td>
<td style="text-align:right;">
146.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.96
</td>
<td style="text-align:right;">
162.4
</td>
<td style="text-align:right;">
12.90–13.06
</td>
<td style="text-align:right;">
150
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.186
</td>
<td style="text-align:right;">
0.250
</td>
<td style="text-align:right;">
0.431
</td>
<td style="text-align:right;">
-0.181
</td>
<td style="text-align:right;">
112.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.34
</td>
<td style="text-align:right;">
609.3
</td>
<td style="text-align:right;">
13.13–13.58
</td>
<td style="text-align:right;">
440
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
0.971
</td>
<td style="text-align:right;">
0.702
</td>
<td style="text-align:right;">
0.269
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.13
</td>
<td style="text-align:right;">
250.1
</td>
<td style="text-align:right;">
12.10–12.13
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.154
</td>
<td style="text-align:right;">
0.964
</td>
<td style="text-align:right;">
0.800
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:right;">
250.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.60
</td>
<td style="text-align:right;">
65.8
</td>
<td style="text-align:right;">
6.37–6.87
</td>
<td style="text-align:right;">
496
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
0.165
</td>
<td style="text-align:right;">
0.283
</td>
<td style="text-align:right;">
-0.118
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
5.37
</td>
<td style="text-align:right;">
58.7
</td>
<td style="text-align:right;">
5.09–5.37
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.097
</td>
<td style="text-align:right;">
0.166
</td>
<td style="text-align:right;">
0.255
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
54.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.70
</td>
<td style="text-align:right;">
50.9
</td>
<td style="text-align:right;">
8.13–9.04
</td>
<td style="text-align:right;">
906
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.090
</td>
<td style="text-align:right;">
0.159
</td>
<td style="text-align:right;">
0.249
</td>
<td style="text-align:right;">
-0.090
</td>
<td style="text-align:right;">
46.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
10.05
</td>
<td style="text-align:right;">
43.3
</td>
<td style="text-align:right;">
10.05–10.13
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
0.149
</td>
<td style="text-align:right;">
0.265
</td>
<td style="text-align:right;">
-0.117
</td>
<td style="text-align:right;">
43.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.08
</td>
<td style="text-align:right;">
31.2
</td>
<td style="text-align:right;">
4.06–4.08
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.078
</td>
<td style="text-align:right;">
0.188
</td>
<td style="text-align:right;">
0.295
</td>
<td style="text-align:right;">
-0.107
</td>
<td style="text-align:right;">
31.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
16.36
</td>
<td style="text-align:right;">
114.9
</td>
<td style="text-align:right;">
16.16–16.49
</td>
<td style="text-align:right;">
330
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.157
</td>
<td style="text-align:right;">
0.263
</td>
<td style="text-align:right;">
0.425
</td>
<td style="text-align:right;">
-0.162
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.65
</td>
<td style="text-align:right;">
97.9
</td>
<td style="text-align:right;">
14.55–14.82
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.145
</td>
<td style="text-align:right;">
0.253
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
-0.154
</td>
<td style="text-align:right;">
86.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.55
</td>
<td style="text-align:right;">
63.1
</td>
<td style="text-align:right;">
13.52–13.55
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.120
</td>
<td style="text-align:right;">
0.339
</td>
<td style="text-align:right;">
0.439
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
63.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
17.54
</td>
<td style="text-align:right;">
41.2
</td>
<td style="text-align:right;">
17.49–17.60
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.098
</td>
<td style="text-align:right;">
0.306
</td>
<td style="text-align:right;">
0.406
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
29.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.51
</td>
<td style="text-align:right;">
27.1
</td>
<td style="text-align:right;">
12.45–12.51
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.080
</td>
<td style="text-align:right;">
0.373
</td>
<td style="text-align:right;">
0.453
</td>
<td style="text-align:right;">
-0.080
</td>
<td style="text-align:right;">
27.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.37
</td>
<td style="text-align:right;">
134.5
</td>
<td style="text-align:right;">
11.19–11.55
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.128
</td>
<td style="text-align:right;">
0.778
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
-0.131
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.52
</td>
<td style="text-align:right;">
96.8
</td>
<td style="text-align:right;">
9.37–10.13
</td>
<td style="text-align:right;">
756
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.103
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
0.917
</td>
<td style="text-align:right;">
-0.100
</td>
<td style="text-align:right;">
77.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.62
</td>
<td style="text-align:right;">
84.9
</td>
<td style="text-align:right;">
7.51–7.87
</td>
<td style="text-align:right;">
351
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
0.852
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
67.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.21
</td>
<td style="text-align:right;">
49.1
</td>
<td style="text-align:right;">
12.93–13.33
</td>
<td style="text-align:right;">
400
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
0.787
</td>
<td style="text-align:right;">
0.871
</td>
<td style="text-align:right;">
-0.084
</td>
<td style="text-align:right;">
41.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs par1
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
29.4
</td>
<td style="text-align:right;">
0.02–0.04
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.061
</td>
<td style="text-align:right;">
0.876
</td>
<td style="text-align:right;">
0.833
</td>
<td style="text-align:right;">
0.043
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.23
</td>
<td style="text-align:right;">
75.4
</td>
<td style="text-align:right;">
13.13–13.37
</td>
<td style="text-align:right;">
235
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.799
</td>
<td style="text-align:right;">
0.686
</td>
<td style="text-align:right;">
0.112
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
14.37
</td>
<td style="text-align:right;">
48.3
</td>
<td style="text-align:right;">
14.37–14.42
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
0.763
</td>
<td style="text-align:right;">
0.659
</td>
<td style="text-align:right;">
0.104
</td>
<td style="text-align:right;">
48.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.10
</td>
<td style="text-align:right;">
46.2
</td>
<td style="text-align:right;">
10.99–11.21
</td>
<td style="text-align:right;">
217
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.066
</td>
<td style="text-align:right;">
0.918
</td>
<td style="text-align:right;">
0.849
</td>
<td style="text-align:right;">
0.069
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
9.94
</td>
<td style="text-align:right;">
22.2
</td>
<td style="text-align:right;">
9.76–9.99
</td>
<td style="text-align:right;">
226
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.041
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
0.891
</td>
<td style="text-align:right;">
0.043
</td>
<td style="text-align:right;">
21.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.60
</td>
<td style="text-align:right;">
20.2
</td>
<td style="text-align:right;">
7.48–7.78
</td>
<td style="text-align:right;">
295
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.038
</td>
<td style="text-align:right;">
0.946
</td>
<td style="text-align:right;">
0.896
</td>
<td style="text-align:right;">
0.050
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.77
</td>
<td style="text-align:right;">
127.4
</td>
<td style="text-align:right;">
10.72–10.83
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.166
</td>
<td style="text-align:right;">
0.261
</td>
<td style="text-align:right;">
0.423
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.96
</td>
<td style="text-align:right;">
125.1
</td>
<td style="text-align:right;">
12.88–13.11
</td>
<td style="text-align:right;">
231
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.164
</td>
<td style="text-align:right;">
0.247
</td>
<td style="text-align:right;">
0.411
</td>
<td style="text-align:right;">
-0.164
</td>
<td style="text-align:right;">
71.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.88
</td>
<td style="text-align:right;">
123.9
</td>
<td style="text-align:right;">
11.83–11.88
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.166
</td>
<td style="text-align:right;">
0.274
</td>
<td style="text-align:right;">
0.454
</td>
<td style="text-align:right;">
-0.180
</td>
<td style="text-align:right;">
71.6
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.34
</td>
<td style="text-align:right;">
84.4
</td>
<td style="text-align:right;">
9.26–9.48
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.131
</td>
<td style="text-align:right;">
0.224
</td>
<td style="text-align:right;">
0.359
</td>
<td style="text-align:right;">
-0.135
</td>
<td style="text-align:right;">
63.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
8.02
</td>
<td style="text-align:right;">
72.4
</td>
<td style="text-align:right;">
7.88–8.18
</td>
<td style="text-align:right;">
296
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
0.242
</td>
<td style="text-align:right;">
0.367
</td>
<td style="text-align:right;">
-0.125
</td>
<td style="text-align:right;">
56.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
309.4
</td>
<td style="text-align:right;">
13.13–13.50
</td>
<td style="text-align:right;">
360
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.166
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
0.794
</td>
<td style="text-align:right;">
0.175
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.13
</td>
<td style="text-align:right;">
123.7
</td>
<td style="text-align:right;">
12.11–12.13
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
0.964
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
123.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
32.4
</td>
<td style="text-align:right;">
0.32–0.56
</td>
<td style="text-align:right;">
231
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
0.516
</td>
<td style="text-align:right;">
0.592
</td>
<td style="text-align:right;">
-0.076
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
13.59
</td>
<td style="text-align:right;">
74.1
</td>
<td style="text-align:right;">
13.22–13.85
</td>
<td style="text-align:right;">
622
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:right;">
0.329
</td>
<td style="text-align:right;">
0.213
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.89
</td>
<td style="text-align:right;">
57.3
</td>
<td style="text-align:right;">
15.75–16.10
</td>
<td style="text-align:right;">
346
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.397
</td>
<td style="text-align:right;">
0.296
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
35.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.22
</td>
<td style="text-align:right;">
56.9
</td>
<td style="text-align:right;">
12.10–12.22
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
0.222
</td>
<td style="text-align:right;">
0.145
</td>
<td style="text-align:right;">
0.077
</td>
<td style="text-align:right;">
56.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
17.49
</td>
<td style="text-align:right;">
45.1
</td>
<td style="text-align:right;">
17.10–17.49
</td>
<td style="text-align:right;">
396
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
0.420
</td>
<td style="text-align:right;">
0.337
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
44.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
11.10
</td>
<td style="text-align:right;">
35.2
</td>
<td style="text-align:right;">
10.91–11.10
</td>
<td style="text-align:right;">
185
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
0.198
</td>
<td style="text-align:right;">
0.119
</td>
<td style="text-align:right;">
0.079
</td>
<td style="text-align:right;">
35.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:right;">
127.0
</td>
<td style="text-align:right;">
1.22–1.61
</td>
<td style="text-align:right;">
386
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
0.586
</td>
<td style="text-align:right;">
0.768
</td>
<td style="text-align:right;">
-0.182
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
16.33
</td>
<td style="text-align:right;">
114.9
</td>
<td style="text-align:right;">
16.11–16.52
</td>
<td style="text-align:right;">
405
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.158
</td>
<td style="text-align:right;">
0.289
</td>
<td style="text-align:right;">
0.440
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
2.65
</td>
<td style="text-align:right;">
109.5
</td>
<td style="text-align:right;">
2.61–2.94
</td>
<td style="text-align:right;">
326
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.148
</td>
<td style="text-align:right;">
0.621
</td>
<td style="text-align:right;">
0.787
</td>
<td style="text-align:right;">
-0.166
</td>
<td style="text-align:right;">
108.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
3.94
</td>
<td style="text-align:right;">
106.3
</td>
<td style="text-align:right;">
3.94–4.09
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.146
</td>
<td style="text-align:right;">
0.625
</td>
<td style="text-align:right;">
0.775
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
91.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
106.0
</td>
<td style="text-align:right;">
0.01–0.22
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.146
</td>
<td style="text-align:right;">
0.633
</td>
<td style="text-align:right;">
0.777
</td>
<td style="text-align:right;">
-0.144
</td>
<td style="text-align:right;">
91.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.21
</td>
<td style="text-align:right;">
130.0
</td>
<td style="text-align:right;">
11.04–11.51
</td>
<td style="text-align:right;">
466
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.126
</td>
<td style="text-align:right;">
0.797
</td>
<td style="text-align:right;">
0.912
</td>
<td style="text-align:right;">
-0.115
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.99
</td>
<td style="text-align:right;">
124.1
</td>
<td style="text-align:right;">
9.86–10.04
</td>
<td style="text-align:right;">
176
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.116
</td>
<td style="text-align:right;">
0.807
</td>
<td style="text-align:right;">
0.928
</td>
<td style="text-align:right;">
-0.120
</td>
<td style="text-align:right;">
89.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.02
</td>
<td style="text-align:right;">
89.9
</td>
<td style="text-align:right;">
12.87–13.18
</td>
<td style="text-align:right;">
305
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.105
</td>
<td style="text-align:right;">
0.802
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
-0.108
</td>
<td style="text-align:right;">
70.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
7.59
</td>
<td style="text-align:right;">
76.8
</td>
<td style="text-align:right;">
7.49–7.84
</td>
<td style="text-align:right;">
341
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.084
</td>
<td style="text-align:right;">
0.848
</td>
<td style="text-align:right;">
0.919
</td>
<td style="text-align:right;">
-0.072
</td>
<td style="text-align:right;">
64.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs rpn12
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.85
</td>
<td style="text-align:right;">
75.5
</td>
<td style="text-align:right;">
8.84–8.85
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.086
</td>
<td style="text-align:right;">
0.853
</td>
<td style="text-align:right;">
0.948
</td>
<td style="text-align:right;">
-0.095
</td>
<td style="text-align:right;">
75.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
13.60
</td>
<td style="text-align:right;">
88.7
</td>
<td style="text-align:right;">
13.23–14.04
</td>
<td style="text-align:right;">
805
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.180
</td>
<td style="text-align:right;">
0.807
</td>
<td style="text-align:right;">
0.631
</td>
<td style="text-align:right;">
0.176
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
7.60
</td>
<td style="text-align:right;">
48.5
</td>
<td style="text-align:right;">
7.26–8.04
</td>
<td style="text-align:right;">
776
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.096
</td>
<td style="text-align:right;">
0.943
</td>
<td style="text-align:right;">
0.822
</td>
<td style="text-align:right;">
0.121
</td>
<td style="text-align:right;">
31.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
11.40
</td>
<td style="text-align:right;">
45.4
</td>
<td style="text-align:right;">
11.16–11.98
</td>
<td style="text-align:right;">
810
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
0.920
</td>
<td style="text-align:right;">
0.845
</td>
<td style="text-align:right;">
0.076
</td>
<td style="text-align:right;">
41.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.46
</td>
<td style="text-align:right;">
43.2
</td>
<td style="text-align:right;">
2.26–2.76
</td>
<td style="text-align:right;">
502
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.122
</td>
<td style="text-align:right;">
0.816
</td>
<td style="text-align:right;">
0.688
</td>
<td style="text-align:right;">
0.128
</td>
<td style="text-align:right;">
25.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.26
</td>
<td style="text-align:right;">
39.5
</td>
<td style="text-align:right;">
6.05–6.26
</td>
<td style="text-align:right;">
205
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.088
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
0.834
</td>
<td style="text-align:right;">
0.100
</td>
<td style="text-align:right;">
39.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
9.47
</td>
<td style="text-align:right;">
116.4
</td>
<td style="text-align:right;">
9.24–9.55
</td>
<td style="text-align:right;">
310
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.214
</td>
<td style="text-align:right;">
0.244
</td>
<td style="text-align:right;">
0.489
</td>
<td style="text-align:right;">
-0.245
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.43
</td>
<td style="text-align:right;">
115.2
</td>
<td style="text-align:right;">
7.29–7.55
</td>
<td style="text-align:right;">
260
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.214
</td>
<td style="text-align:right;">
0.260
</td>
<td style="text-align:right;">
0.479
</td>
<td style="text-align:right;">
-0.219
</td>
<td style="text-align:right;">
84.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
11.36
</td>
<td style="text-align:right;">
110.5
</td>
<td style="text-align:right;">
11.22–11.50
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.210
</td>
<td style="text-align:right;">
0.236
</td>
<td style="text-align:right;">
0.427
</td>
<td style="text-align:right;">
-0.191
</td>
<td style="text-align:right;">
80.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
12.92
</td>
<td style="text-align:right;">
82.1
</td>
<td style="text-align:right;">
12.86–13.00
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.180
</td>
<td style="text-align:right;">
0.283
</td>
<td style="text-align:right;">
0.453
</td>
<td style="text-align:right;">
-0.170
</td>
<td style="text-align:right;">
35.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.45
</td>
<td style="text-align:right;">
60.2
</td>
<td style="text-align:right;">
4.35–4.55
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.158
</td>
<td style="text-align:right;">
0.326
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:right;">
-0.174
</td>
<td style="text-align:right;">
27.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
71.9
</td>
<td style="text-align:right;">
13.12–13.50
</td>
<td style="text-align:right;">
375
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
0.886
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
57.5
</td>
<td style="text-align:right;">
0.42–0.69
</td>
<td style="text-align:right;">
261
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
0.491
</td>
<td style="text-align:right;">
0.645
</td>
<td style="text-align:right;">
-0.155
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
1.81
</td>
<td style="text-align:right;">
39.3
</td>
<td style="text-align:right;">
1.70–1.95
</td>
<td style="text-align:right;">
245
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.121
</td>
<td style="text-align:right;">
0.571
</td>
<td style="text-align:right;">
0.711
</td>
<td style="text-align:right;">
-0.139
</td>
<td style="text-align:right;">
12.1
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
12.12
</td>
<td style="text-align:right;">
29.0
</td>
<td style="text-align:right;">
12.09–12.12
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.062
</td>
<td style="text-align:right;">
0.965
</td>
<td style="text-align:right;">
0.897
</td>
<td style="text-align:right;">
0.068
</td>
<td style="text-align:right;">
29.0
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
9.88
</td>
<td style="text-align:right;">
26.8
</td>
<td style="text-align:right;">
9.67–10.03
</td>
<td style="text-align:right;">
356
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.058
</td>
<td style="text-align:right;">
0.965
</td>
<td style="text-align:right;">
0.911
</td>
<td style="text-align:right;">
0.055
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
15.83
</td>
<td style="text-align:right;">
46.0
</td>
<td style="text-align:right;">
15.68–15.96
</td>
<td style="text-align:right;">
271
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.133
</td>
<td style="text-align:right;">
0.414
</td>
<td style="text-align:right;">
0.308
</td>
<td style="text-align:right;">
0.106
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
13.80
</td>
<td style="text-align:right;">
42.3
</td>
<td style="text-align:right;">
13.53–14.21
</td>
<td style="text-align:right;">
670
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.326
</td>
<td style="text-align:right;">
0.202
</td>
<td style="text-align:right;">
0.124
</td>
<td style="text-align:right;">
24.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
12.53
</td>
<td style="text-align:right;">
33.3
</td>
<td style="text-align:right;">
12.49–12.53
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.094
</td>
<td style="text-align:right;">
0.248
</td>
<td style="text-align:right;">
0.155
</td>
<td style="text-align:right;">
0.093
</td>
<td style="text-align:right;">
33.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
16.96
</td>
<td style="text-align:right;">
20.3
</td>
<td style="text-align:right;">
16.96–17.49
</td>
<td style="text-align:right;">
537
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
0.435
</td>
<td style="text-align:right;">
0.348
</td>
<td style="text-align:right;">
0.086
</td>
<td style="text-align:right;">
20.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
7.7
</td>
<td style="text-align:right;">
0.45–0.79
</td>
<td style="text-align:right;">
343
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.058
</td>
<td style="text-align:right;">
0.502
</td>
<td style="text-align:right;">
0.539
</td>
<td style="text-align:right;">
-0.036
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.74
</td>
<td style="text-align:right;">
121.6
</td>
<td style="text-align:right;">
13.59–13.97
</td>
<td style="text-align:right;">
382
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.223
</td>
<td style="text-align:right;">
0.349
</td>
<td style="text-align:right;">
0.562
</td>
<td style="text-align:right;">
-0.213
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
14.98
</td>
<td style="text-align:right;">
81.2
</td>
<td style="text-align:right;">
14.98–15.05
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.179
</td>
<td style="text-align:right;">
0.257
</td>
<td style="text-align:right;">
0.421
</td>
<td style="text-align:right;">
-0.165
</td>
<td style="text-align:right;">
75.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.08
</td>
<td style="text-align:right;">
80.6
</td>
<td style="text-align:right;">
11.85–12.36
</td>
<td style="text-align:right;">
506
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.183
</td>
<td style="text-align:right;">
0.402
</td>
<td style="text-align:right;">
0.569
</td>
<td style="text-align:right;">
-0.167
</td>
<td style="text-align:right;">
75.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
6.48
</td>
<td style="text-align:right;">
75.1
</td>
<td style="text-align:right;">
6.36–6.59
</td>
<td style="text-align:right;">
226
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.163
</td>
<td style="text-align:right;">
0.574
</td>
<td style="text-align:right;">
0.725
</td>
<td style="text-align:right;">
-0.151
</td>
<td style="text-align:right;">
18.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
4.88
</td>
<td style="text-align:right;">
73.3
</td>
<td style="text-align:right;">
4.80–5.10
</td>
<td style="text-align:right;">
301
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.160
</td>
<td style="text-align:right;">
0.640
</td>
<td style="text-align:right;">
0.730
</td>
<td style="text-align:right;">
-0.089
</td>
<td style="text-align:right;">
63.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
9.81
</td>
<td style="text-align:right;">
95.0
</td>
<td style="text-align:right;">
9.55–10.12
</td>
<td style="text-align:right;">
565
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.119
</td>
<td style="text-align:right;">
0.824
</td>
<td style="text-align:right;">
0.939
</td>
<td style="text-align:right;">
-0.115
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
11.18
</td>
<td style="text-align:right;">
80.9
</td>
<td style="text-align:right;">
11.12–11.42
</td>
<td style="text-align:right;">
296
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.122
</td>
<td style="text-align:right;">
0.803
</td>
<td style="text-align:right;">
0.906
</td>
<td style="text-align:right;">
-0.102
</td>
<td style="text-align:right;">
62.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
8.55
</td>
<td style="text-align:right;">
59.6
</td>
<td style="text-align:right;">
7.55–8.55
</td>
<td style="text-align:right;">
1004
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.091
</td>
<td style="text-align:right;">
0.856
</td>
<td style="text-align:right;">
0.939
</td>
<td style="text-align:right;">
-0.083
</td>
<td style="text-align:right;">
55.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
12.42
</td>
<td style="text-align:right;">
30.8
</td>
<td style="text-align:right;">
12.42–12.53
</td>
<td style="text-align:right;">
106
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.073
</td>
<td style="text-align:right;">
0.827
</td>
<td style="text-align:right;">
0.912
</td>
<td style="text-align:right;">
-0.085
</td>
<td style="text-align:right;">
30.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
pos1 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.53
</td>
<td style="text-align:right;">
23.8
</td>
<td style="text-align:right;">
13.53–13.59
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.076
</td>
<td style="text-align:right;">
0.760
</td>
<td style="text-align:right;">
0.829
</td>
<td style="text-align:right;">
-0.069
</td>
<td style="text-align:right;">
22.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:right;">
29.7
</td>
<td style="text-align:right;">
2.06–2.40
</td>
<td style="text-align:right;">
336
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
0.792
</td>
<td style="text-align:right;">
0.657
</td>
<td style="text-align:right;">
0.134
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
29.0
</td>
<td style="text-align:right;">
0.01–0.03
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.111
</td>
<td style="text-align:right;">
0.648
</td>
<td style="text-align:right;">
0.532
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
14.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
3.41
</td>
<td style="text-align:right;">
25.0
</td>
<td style="text-align:right;">
3.41–3.48
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.078
</td>
<td style="text-align:right;">
0.890
</td>
<td style="text-align:right;">
0.819
</td>
<td style="text-align:right;">
0.071
</td>
<td style="text-align:right;">
20.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
6.61
</td>
<td style="text-align:right;">
22.3
</td>
<td style="text-align:right;">
6.42–6.87
</td>
<td style="text-align:right;">
450
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.068
</td>
<td style="text-align:right;">
0.929
</td>
<td style="text-align:right;">
0.830
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
8.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
I
</td>
<td style="text-align:right;">
12.51
</td>
<td style="text-align:right;">
21.1
</td>
<td style="text-align:right;">
12.43–12.58
</td>
<td style="text-align:right;">
150
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.082
</td>
<td style="text-align:right;">
0.819
</td>
<td style="text-align:right;">
0.714
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
7.20
</td>
<td style="text-align:right;">
23.3
</td>
<td style="text-align:right;">
7.07–7.29
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.099
</td>
<td style="text-align:right;">
0.323
</td>
<td style="text-align:right;">
0.428
</td>
<td style="text-align:right;">
-0.105
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
10.42
</td>
<td style="text-align:right;">
19.1
</td>
<td style="text-align:right;">
10.36–10.49
</td>
<td style="text-align:right;">
125
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.091
</td>
<td style="text-align:right;">
0.385
</td>
<td style="text-align:right;">
0.525
</td>
<td style="text-align:right;">
-0.140
</td>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
8.98
</td>
<td style="text-align:right;">
18.5
</td>
<td style="text-align:right;">
8.93–9.07
</td>
<td style="text-align:right;">
136
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.089
</td>
<td style="text-align:right;">
0.379
</td>
<td style="text-align:right;">
0.434
</td>
<td style="text-align:right;">
-0.054
</td>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
4.93
</td>
<td style="text-align:right;">
8.4
</td>
<td style="text-align:right;">
4.86–5.03
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.061
</td>
<td style="text-align:right;">
0.469
</td>
<td style="text-align:right;">
0.504
</td>
<td style="text-align:right;">
-0.035
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
II
</td>
<td style="text-align:right;">
1.78
</td>
<td style="text-align:right;">
7.9
</td>
<td style="text-align:right;">
1.76–1.79
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.046
</td>
<td style="text-align:right;">
0.896
</td>
<td style="text-align:right;">
0.853
</td>
<td style="text-align:right;">
0.043
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
8.04
</td>
<td style="text-align:right;">
41.2
</td>
<td style="text-align:right;">
7.97–8.16
</td>
<td style="text-align:right;">
186
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.057
</td>
<td style="text-align:right;">
0.992
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
0.059
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
6.88
</td>
<td style="text-align:right;">
30.6
</td>
<td style="text-align:right;">
6.80–6.93
</td>
<td style="text-align:right;">
126
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.049
</td>
<td style="text-align:right;">
0.989
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:right;">
0.049
</td>
<td style="text-align:right;">
12.4
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
5.80
</td>
<td style="text-align:right;">
25.5
</td>
<td style="text-align:right;">
5.75–5.80
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.047
</td>
<td style="text-align:right;">
0.984
</td>
<td style="text-align:right;">
0.930
</td>
<td style="text-align:right;">
0.054
</td>
<td style="text-align:right;">
25.5
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
13.33
</td>
<td style="text-align:right;">
17.2
</td>
<td style="text-align:right;">
13.16–13.51
</td>
<td style="text-align:right;">
346
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.061
</td>
<td style="text-align:right;">
0.803
</td>
<td style="text-align:right;">
0.890
</td>
<td style="text-align:right;">
-0.086
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
III
</td>
<td style="text-align:right;">
4.75
</td>
<td style="text-align:right;">
13.3
</td>
<td style="text-align:right;">
4.57–4.75
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.038
</td>
<td style="text-align:right;">
0.958
</td>
<td style="text-align:right;">
0.927
</td>
<td style="text-align:right;">
0.031
</td>
<td style="text-align:right;">
13.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
6.03
</td>
<td style="text-align:right;">
14.6
</td>
<td style="text-align:right;">
5.87–6.17
</td>
<td style="text-align:right;">
296
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.059
</td>
<td style="text-align:right;">
0.132
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
-0.039
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
8.30
</td>
<td style="text-align:right;">
9.6
</td>
<td style="text-align:right;">
8.05–8.70
</td>
<td style="text-align:right;">
650
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.045
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.144
</td>
<td style="text-align:right;">
-0.029
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
9.70
</td>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
9.70–9.76
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.036
</td>
<td style="text-align:right;">
0.115
</td>
<td style="text-align:right;">
0.152
</td>
<td style="text-align:right;">
-0.037
</td>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
0.59–0.95
</td>
<td style="text-align:right;">
356
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
-0.051
</td>
<td style="text-align:right;">
0.521
</td>
<td style="text-align:right;">
0.523
</td>
<td style="text-align:right;">
-0.002
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
IV
</td>
<td style="text-align:right;">
4.87
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
4.83–4.87
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.038
</td>
<td style="text-align:right;">
0.141
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
-0.019
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
20.49
</td>
<td style="text-align:right;">
22.2
</td>
<td style="text-align:right;">
20.48–20.49
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.011
</td>
<td style="text-align:right;">
0.011
</td>
<td style="text-align:right;">
0.012
</td>
<td style="text-align:right;">
-0.002
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
12.34
</td>
<td style="text-align:right;">
21.3
</td>
<td style="text-align:right;">
12.06–12.53
</td>
<td style="text-align:right;">
468
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-0.096
</td>
<td style="text-align:right;">
0.504
</td>
<td style="text-align:right;">
0.584
</td>
<td style="text-align:right;">
-0.080
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
10.31
</td>
<td style="text-align:right;">
19.8
</td>
<td style="text-align:right;">
10.20–10.46
</td>
<td style="text-align:right;">
256
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-0.089
</td>
<td style="text-align:right;">
0.575
</td>
<td style="text-align:right;">
0.664
</td>
<td style="text-align:right;">
-0.088
</td>
<td style="text-align:right;">
11.8
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
1.28
</td>
<td style="text-align:right;">
19.3
</td>
<td style="text-align:right;">
1.18–1.43
</td>
<td style="text-align:right;">
245
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
0.738
</td>
<td style="text-align:right;">
0.654
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
V
</td>
<td style="text-align:right;">
13.54
</td>
<td style="text-align:right;">
16.3
</td>
<td style="text-align:right;">
13.54–13.60
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
-0.084
</td>
<td style="text-align:right;">
0.452
</td>
<td style="text-align:right;">
0.508
</td>
<td style="text-align:right;">
-0.056
</td>
<td style="text-align:right;">
16.3
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
15.79
</td>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
15.69–15.90
</td>
<td style="text-align:right;">
201
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.047
</td>
<td style="text-align:right;">
0.783
</td>
<td style="text-align:right;">
0.725
</td>
<td style="text-align:right;">
0.058
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
13.02
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
12.95–13.07
</td>
<td style="text-align:right;">
116
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.034
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
0.834
</td>
<td style="text-align:right;">
0.076
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
TRUE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516
</td>
<td style="text-align:left;">
rpn12 vs vha5
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;">
14.69
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
14.68–14.69
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.035
</td>
<td style="text-align:right;">
0.801
</td>
<td style="text-align:right;">
0.748
</td>
<td style="text-align:right;">
0.053
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
FALSE
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

</div>

# III. The interval

<div class="scale">

13.8 Mb → 37 kb

</div>

<div class="claim">

A near-isogenic line series on the right arm of chromosome III fine-maps
the general RNAi-response QTL to 37 kb, and gives an allelic series in
hatching.

</div>

## Figure 3 — NILs fine-map the chromosome III QTL to 37 kb

<div class="meta">

**Script** `scripts/Figure3_quad.R` · 50% *pos-1* RNAi<br> **Cross**
JU1793 × JU2466, HT115 vs *pos-1*<br> **chrIII peak** 13.784 Mb, LOD
139.7 · **resolved** 13.658–13.695 Mb (37 kb)

</div>

<div class="plate">

<img src="plots/Figure3_quad.png" alt="Three panels: pooled phenotype histogram with both cross parents marked, parental allele frequency along chromosome III, and the NIL series with introgression genotypes beside embryo hatching, rows unlabelled." width="100%" />
<p class="filecap">
plots/Figure3_quad.pdf · .png
</p>

</div>

Three panels, arranged so the figure reads in the order the argument
runs: two strains at opposite ends of the pooled panel, crossing them
drives one parent’s haplotype toward the QTL, and NILs carrying pieces
of it give an allelic series whose genotypes and phenotypes sit side by
side.

<div class="panel">

<span class="pl">A</span> The pooled *pos-1* phenotype distribution as a
histogram, with both cross parents marked by equal-height lollipops —
JU1793 vst 0.121, rank 83/84; JU2466 −0.051, rank 10/84. This makes
visible that the cross was built from opposite ends of the pooled panel.

</div>

<div class="panel">

<span class="pl">B</span> **Parental allele frequency** along the right
arm of chromosome III, 8 Mb to the telomere, in the JU1793×JU2466 F2
pool under *pos-1* RNAi. JU1793’s haplotype frequency is filled from
below in orange and JU2466’s above it in teal, so the two sum to 1 and
the panel reads as which parent occupies the pool at each position.
Frequencies are count-weighted within 50 kb bins, then shown as a
centred 250 kb rolling mean. The JU1793 fraction rises from `42%` at 8
Mb to `79%` at the telomere. The solid line is the HT115 control pool; a
dashed line marks 50%.

</div>

<div class="aside">

<span class="ch">Why a frequency trace instead of a LOD trace</span>

A LOD trace answers whether there is a QTL here, which Figure 2 already
answers. It says nothing about **which parent’s allele** the selection
favoured, and that is the claim this panel is making: *sid-2* sits at
13.68 Mb, and the resistant parent’s haplotype is what sweeps toward it.

The control is what makes this selection rather than a segregation
artefact. Over the same interval the HT115 pool runs the **other** way,
`42%` → `33%`, while the *pos-1* pool goes `42%` → `79%`.

Two details that matter for trusting it. Frequencies are
**count-weighted** — counts summed within each bin and the frequency
taken from the sums — so a marker with 200 reads counts for more than
one with 4; averaging per-marker frequencies would let the shallowest
markers pull the trace around. And the parent assignment is taken from
the cross export’s own README (`p1`/`p2` are the JU1793 and JU2466
haplotypes), not inferred: a swap would invert the panel and still look
plausible, so the script asserts that the *pos-1* pool ends up more
JU1793 than the control at the right end and stops if it does not.

</div>

<div class="panel">

<span class="pl">C</span> The NIL series with **genotype and phenotype
in one panel**, one row per strain. The rows carry **no strain names**;
from the bottom up they are JU1793, wSZ196, wSZ191, wSZ176, JU2466 — the
two parents being the two single-colour rows, and the per-strain values
listed below, so each row is identifiable from its genotype and hatching
together. Figure S10 draws the same series with the names on. Left, the
introgressions on the right arm of chromosome III, 13.635 Mb to the
telomere, JU1793 genotype in orange and JU2466 in teal — the smaller
half of the panel, on thin bars, because it carries two breakpoints
where the hatching carries a five-level series with intervals. Right,
embryos hatched under 50% *pos-1* RNAi on the same rows, one plate per
strain with Wilson 95% binomial intervals. The two halves share the row
axis and carry separate x scales, labelled beneath each.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Strain
</th>
<th style="text-align:right;">
Embryos
</th>
<th style="text-align:right;">
Hatched
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:right;">
179
</td>
<td style="text-align:right;">
99.4%
</td>
</tr>
<tr>
<td style="text-align:left;">
wSZ196
</td>
<td style="text-align:right;">
261
</td>
<td style="text-align:right;">
97.3%
</td>
</tr>
<tr>
<td style="text-align:left;">
wSZ191
</td>
<td style="text-align:right;">
252
</td>
<td style="text-align:right;">
79.4%
</td>
</tr>
<tr>
<td style="text-align:left;">
wSZ176
</td>
<td style="text-align:right;">
272
</td>
<td style="text-align:right;">
58.5%
</td>
</tr>
<tr>
<td style="text-align:left;">
JU2466
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
35.7%
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from supplemental_data/hatching_assays/nil_series_hatching.tsv

</div>

<div class="aside">

<span class="ch">The introgressions run to the end of the
chromosome</span>

This is what distinguishes the strains, so the panel says it rather than
leaving it to be inferred. Five of the six NILs carry a JU2466 segment
that begins at a breakpoint and continues to the chromosome III terminus
at **13,783,801 bp**, which the right-hand edge of the genotype track
marks. **wSZ191 is the exception** — an internal 13.658–13.695 Mb
segment that stops short.

That contrast is the fine-mapping argument: wSZ191 carries the resolved
interval and nothing distal to it, wSZ196 carries everything distal and
not the interval, and their hatching differs (`79.4%` against `97.3%`).
So compressing the axis to the breakpoints alone would crop exactly the
fact the panel exists to show; the window keeps its right edge at the
terminus, and the compression is in the left flank and the row height
instead.

The shaded band with dotted edges is the interval the series resolves,
**13.658–13.695 Mb** — the region wSZ191 carries and wSZ196 does not.

</div>

<div class="panel">

<span class="pl">Previously two panels</span> C (hatching) and D
(genotypes) drew the same five strains as rows twice, in two panels with
two x axes, so pairing a genotype with its phenotype meant carrying a
row position across a panel boundary. The two-panel functions are still
in `Figure3_common.R` — `Figure3_chrIII.R` and the supplements use them.

</div>

<div class="caveat">

<span class="ch">Two caveats, the second worth a sentence in the
text</span>

1.  One plate per strain per condition, so the hatching intervals in C
    describe *counting* uncertainty, not between-plate variability, and
    no strain is replicated. The ordering should be read as an allelic
    series, not as a set of tested contrasts. The HT115 control arm and
    four further NILs are in the supplement below.
2.  The chromosome III peak marker at 13.784 Mb **is the terminal marker
    of the chromosome**. A scan cannot place a peak past the chromosome
    end, so that position is a boundary artefact rather than a
    localisation, and the NIL interval at 13.658–13.695 Mb is the claim.
    “The peak coincides with the NIL interval” is not a sentence these
    data support.

</div>

## Figure S10 — the complete NIL hatching experiment

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_nil_hatching_full.R` · 50% *pos-1*
RNAi<br> **Supports** Figure 3C · 10 strains, both food conditions

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_nil_hatching_full.png" alt="Horizontal bar chart of embryos hatched for ten NIL strains on both HT115 control and pos-1 RNAi food, with Wilson intervals." width="100%" />
<p class="filecap">
SUPP_FIG_XX_nil_hatching_full
</p>

</div>

All 10 strains on both food conditions at 50% *pos-1* RNAi, strains
labelled and ordered by *pos-1* hatching, with embryo counts on each bar
and Wilson 95% binomial intervals. **Control hatching is 97–100% for
every strain**, which is what makes the *pos-1* differences in Figure 3C
attributable to the RNAi rather than to the introgressions themselves.

Adds four NILs that Figure 3 does not use — wSZ192, wSZ193, wSZ194,
wSZ195 — and wSZ153.

<div class="caveat">

<span class="ch">Why those four are not in Figure 3D</span>

None of the four has introgression coordinates in the NIL ranges file,
so their genotypes cannot be drawn. **wSZ192 and wSZ193 in particular
would sharpen the interval if their breakpoints can be recovered.** One
plate per strain per condition, as in Figure 3C.

</div>

# IV. The residue

<div class="scale">

one residue

</div>

<div class="claim">

Editing *sid-2* residue 96 moves RNAi sensitivity in three genetic
backgrounds and in both directions, and the residue sits on the lumenal
face of the ectodomain alongside residues already known to be required
for dsRNA uptake.

</div>

## Figure 4 — editing *sid-2* residue 96 moves sensitivity in three backgrounds

<div class="meta">

**Script** `scripts/Figure4_sid2.R`<br> **Assets**
`scripts/sid2_ribbon_render.py` · `scripts/sid2_zoom_render.py`<br>
**Variant** chrIII:13,680,248 C>A · T96K<br> **Dose** 50% *pos-1* (A) ·
**25%** *pos-1* (B)

</div>

<div class="plate">

<img src="plots/Figure4_sid2.png" alt="Four panels: residue-96 allele swaps in JU1793 and JU2466, the same swap in N2 at 25% dose, the SID-2 ectodomain ribbon with a zoom on T96, and sid-2 wild coding variation drawn along the protein." width="100%" />
<p class="filecap">
plots/Figure4_sid2.pdf · .png
</p>

</div>

Panels carry letters only; all descriptive text is here.

<div class="panel">

<span class="pl">A</span> Embryos hatched on **50%** *pos-1* RNAi for
reciprocal edits at residue 96 in the two cross parents. JU1793\[96T\],
the resistant parental allele, 94.8% (`n = 213` embryos); JU1793\[96K\]
53.1% (`n = 207`); JU2466\[96K\], the sensitive parental allele, 5.4%
(`n = 204`); JU2466\[96T\] 18.4% (`n = 217`). Parental alleles are in
the strain colours and edited alleles in grey. Fisher’s exact test, each
edit against the unedited allele in the same background: `p = 1.7e-24`
and `p = 3.9e-05`. Editing 96T→K removes most of JU1793’s resistance and
96K→T recovers part of JU2466’s, so residue 96 contributes in both
directions without accounting for the whole difference between the
parents.

</div>

<div class="panel">

<span class="pl">B</span> The same substitution in N2, at **25%**
*pos-1* RNAi. N2 carries 96T and hatches 32.3% (`n = 220`); two
independently derived lines edited to 96K hatch 4.4% (wSZ203, `n = 273`)
and 3.7% (wSZ204, `n = 295`), Fisher `p = 7.8e-17` and `p = 6.2e-19`.
The two edited lines are labelled by genotype and not distinguished on
the axis, because they are independent lines of the same edit. This is
the third genetic background to show the effect and the only one with
two independent edits, and the direction is the same throughout: **96K
sensitive, 96T resistant**.

</div>

<div class="panel">

<span class="pl">C</span> The SID-2 ectodomain coloured by **local net
charge**, AlphaFold3 model, rotated onto the membrane normal so the
intestinal lumen is up; the grey slab is the bilayer. Colour is the net
side-chain charge of every ectodomain residue with a Cα within 12 Å, by
Henderson–Hasselbalch at pH 4.4 — the gut-lumen pH at which SID-2
functions — on a diverging scale saturating at ±`2` e, red negative and
blue positive. Left, residues 21–188 of chain A: T96 (orange), the two
lysines that make its pocket, K93 and K132 (blue), and SID-2’s three
extracellular histidines H32, H168 and H175 (teal). Right, the pocket
enlarged, with Cα distances from T96 of `6.6` Å to K93 and `6.8` Å to
K132 (nearest heavy atoms `4.5` and `4.4` Å).

</div>

<div class="aside">

<span class="ch">Why this panel is about charge and not about
shape</span>

The panel used to be coloured by secondary structure, which showed the
fold — not the argument, and not anything in dispute. The argument is
electrostatic, and it rests on a precedent running the same direction as
our own result.

SID-1’s dsRNA recognition is **electrostatic and sequence-independent**:
basic side chains against the phosphate backbone, no sequence
preference. That transfers as physics without any fold relationship,
which matters because the fold relationship does *not* hold — an
unrelated immunoglobulin domain scores higher against the SID-2 model
(TM `0.427`) than SID-1’s BRD1 does (`0.404`), while the two genuine
BRDs of SID-1 score `0.525` against each other. The resemblance is the
background level for any compact β-sandwich of this size. Nothing in the
manuscript claims otherwise, and nothing should: do not superpose the
8XC1 dsRNA onto this model.

What the charge measurement says: the ectodomain is net **acidic**
overall (`−8.34` e at pH 7.4, `+0.47` e at pH 4.4 across the modelled
21–188), yet T96’s own neighbourhood is `+1.24` e — the **82nd
percentile** of the domain — and T96K takes it to `+2.24` e, the
**98th**. T96 is solvent-exposed at the domain median (SASA `62.4` Å²,
48th percentile), so the side chain is available. The precedent: the
triple His→Arg mutant of McEwan et al. internalised **more** dsRNA than
wild type — replacing pH-dependent positive charge with permanent
positive charge increased uptake. T96K adds permanent positive charge to
the same surface and increases RNAi sensitivity in all three
backgrounds.

Read it as a hypothesis with a direction, not a result. Whether SID-2
contacts dsRNA at all is untested here, charge counting ignores pKa
shifts and glycan shielding (three of the nine N-glycosylation sequons
are in the ectodomain), and the direction is a correlation across three
backgrounds. The discriminating experiment is charge-matched: T96R
should behave like T96K if charge is the mechanism, T96Q (isosteric,
neutral) like the wild type. See Figure S15 for the distribution and
Figure S16 for model confidence.

</div>

<div class="panel">

<span class="pl">D</span> *sid-2* coding variation in the wild
population: the protein drawn vertically with residue 1 at the top and
its topology as a filled bar, and each protein-altering variant labelled
to the right with its CeNDR allele frequency, then two columns giving
the JU1793 and JU2466 allele state. Callout rows are evenly spaced and
joined to the residue by a leader rather than sitting at the residue’s
own height: four of the eight fall between residues 141 and 153 and
would overlap completely at true scale.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Variant
</th>
<th style="text-align:right;">
Residue
</th>
<th style="text-align:right;">
CeNDR AF
</th>
<th style="text-align:right;">
Isotypes
</th>
<th style="text-align:center;">
JU1793
</th>
<th style="text-align:center;">
JU2466
</th>
<th style="text-align:center;">
N2
</th>
<th style="text-align:center;">
XZ1516
</th>
<th style="text-align:center;">
Parents differ
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
V5L
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.015
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:center;">
L
</td>
<td style="text-align:center;">
V
</td>
<td style="text-align:center;">
V
</td>
<td style="text-align:center;">
V
</td>
<td style="text-align:center;">
yes
</td>
</tr>
<tr>
<td style="text-align:left;">
D78A
</td>
<td style="text-align:right;">
78
</td>
<td style="text-align:right;">
0.009
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:center;">
D
</td>
<td style="text-align:center;">
D
</td>
<td style="text-align:center;">
D
</td>
<td style="text-align:center;">
A
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
T96K
</td>
<td style="text-align:right;">
96
</td>
<td style="text-align:right;">
0.453
</td>
<td style="text-align:right;">
537
</td>
<td style="text-align:center;">
T
</td>
<td style="text-align:center;">
K
</td>
<td style="text-align:center;">
T
</td>
<td style="text-align:center;">
K
</td>
<td style="text-align:center;">
yes
</td>
</tr>
<tr>
<td style="text-align:left;">
M141V
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
0.009
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:center;">
M
</td>
<td style="text-align:center;">
M
</td>
<td style="text-align:center;">
M
</td>
<td style="text-align:center;">
V
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
Q144P
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
0.009
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:center;">
Q
</td>
<td style="text-align:center;">
Q
</td>
<td style="text-align:center;">
Q
</td>
<td style="text-align:center;">
P
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
A151I/T
</td>
<td style="text-align:right;">
151
</td>
<td style="text-align:right;">
0.195
</td>
<td style="text-align:right;">
539
</td>
<td style="text-align:center;">
A
</td>
<td style="text-align:center;">
A
</td>
<td style="text-align:center;">
A
</td>
<td style="text-align:center;">
I/T
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
P153T
</td>
<td style="text-align:right;">
153
</td>
<td style="text-align:right;">
0.470
</td>
<td style="text-align:right;">
538
</td>
<td style="text-align:center;">
T
</td>
<td style="text-align:center;">
T
</td>
<td style="text-align:center;">
P
</td>
<td style="text-align:center;">
T
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
L209M
</td>
<td style="text-align:right;">
209
</td>
<td style="text-align:right;">
0.195
</td>
<td style="text-align:right;">
539
</td>
<td style="text-align:center;">
L
</td>
<td style="text-align:center;">
L
</td>
<td style="text-align:center;">
L
</td>
<td style="text-align:center;">
M
</td>
<td style="text-align:center;">
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from supplemental_data/structure/sid2_variants_cendr.tsv

</div>

<div class="tnote">

Panel D is not an exhaustive catalogue: the amino-acid annotation covers
variants segregating among the four cross parents. 81 variants of any
kind segregate in the 3 kb span in CeNDR, and the 8 protein-altering
ones above are what is annotated.

</div>

<div class="aside">

<span class="ch">Linkage, and why it does not confound the cross</span>

T96K and P153T are in near-complete linkage disequilibrium across the
wild population, r² = 0.935: 96K never occurs without 153T, though 153T
occurs without 96K in nine isotypes. **Both cross parents carry 153T** —
JU1793 is 96T/T153 and JU2466 is 96K/T153 — so P153T does not segregate
in the JU1793×JU2466 cross and cannot account for the mapped effect.

Of the 8 annotated coding variants, only 2 differ between the two
parents: **V5L and T96K**. That narrows the cross’s candidate coding
changes to two by inspection.

</div>

<div class="caveat">

<span class="ch">Caveats — all of which belong in the text</span>

1.  **One plate per strain per condition throughout A and B.** The
    intervals are Wilson binomial intervals on that plate’s embryo count
    and describe counting uncertainty, not between-plate variability.
    wSZ203 and wSZ204 in panel B are the closest thing to a biological
    replicate anywhere in the figure.
2.  **Doses differ between A and B** — 50% against 25%. The percentages
    are not comparable between the two panels.
3.  The proximity in panel C is **spatial context, not evidence of a
    shared site**. Three of the four annotated residues are nearer to
    T96 than the domain’s median Cα–Cα distance (23.8 Å), but 41% of the
    ectodomain lies within 20 Å, so three of four landing there gives
    binomial `p = 0.19` and a permutation test on their mean distance
    gives `p = 0.30`.
4.  The prediction is a dimer, but **ipTM is 0.45–0.46 across all five
    models** and 43–46% of the model is called disordered, so no dimer
    is shown and dimerisation is not claimed. Only the ectodomain
    reaches the confident pLDDT band (mean 72.3; T96 79.3).
5.  Secondary structure is assigned from the backbone with the Kabsch &
    Sander hydrogen-bond criterion, not DSSP, which is not installed. It
    is corroborated independently: the criterion calls residues 194–211
    100% helix, exactly the span DeepTMHMM calls the transmembrane
    helix.
6.  T96 is the threonine of an N94-C95-T96 sequon, so T96K removes a
    predicted N-glycosylation site — **but the N94A test of that
    hypothesis fails, and the figure makes no glycosylation claim.** See
    the allele-swaps supplement.
7.  **The variant has no marginal effect across the wild panel.** See
    the in-panel supplement.
8.  Published allele positions come from the UniProt annotation for
    `G5EEV9` as recorded in the earlier structure-modelling work.
    **Verify the residue numbers against UniProt before publication.**

</div>

## Effect sizes across the three chromosome III experiments

**Script** `scripts/effect_size_ladder.R`<br> **Table**
`plots/diagnostics/TABLE_effect_size_ladder.tsv`

Three experiments interrogate the same locus: the NIL series narrows it
to 37 kb, the allele swaps test one residue inside that interval, and
the JU1793 × JU2466 cross measures it as a selection response. Do their
effect sizes agree?

<div class="caveat">

<span class="ch">The two hatching assays do not share a scale</span>

**JU2466 hatches `0.357` under *pos-1* in the NIL series and `0.045` in
the allele-swap series — the same strain, the same nominal 50% dose, a
7.9-fold difference.** JU1793 is `0.994` against `0.948`. So a raw
percentage-point effect from one experiment cannot be set beside one
from the other, and every effect below is therefore given twice: as a
raw difference in hatched fraction, and normalised to the JU1793–JU2466
span measured **within that same experiment**, where 0 is JU2466 and 1
is JU1793.

This is also evidence on an open question: `METHODS.txt` carries a \[TO
FILL\] asking that both 50% dose figures be confirmed, and two assays
nominally at the same dose are not behaving the same way.

</div>

<table>
<thead>
<tr>
<th style="text-align:left;">
Experiment
</th>
<th style="text-align:left;">
Contrast
</th>
<th style="text-align:right;">
Δ hatched
</th>
<th style="text-align:right;">
% of parental span
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NIL series
</td>
<td style="text-align:left;">
JU1793 -> wSZ191 (JU2466 alleles at the 37 kb interval only)
</td>
<td style="text-align:right;">
+0.201 \[+0.149, +0.255\]
</td>
<td style="text-align:right;">
32% \[23, 40\]
</td>
</tr>
<tr>
<td style="text-align:left;">
NIL series
</td>
<td style="text-align:left;">
JU1793 -> wSZ196 (JU2466 alleles distal to the interval only)
</td>
<td style="text-align:right;">
+0.021 \[-0.008, +0.049\]
</td>
<td style="text-align:right;">
3% \[-1, 8\]
</td>
</tr>
<tr>
<td style="text-align:left;">
NIL series
</td>
<td style="text-align:left;">
JU1793 -> wSZ176 (JU2466 alleles across interval + distal)
</td>
<td style="text-align:right;">
+0.410 \[+0.347, +0.469\]
</td>
<td style="text-align:right;">
64% \[55, 74\]
</td>
</tr>
<tr>
<td style="text-align:left;">
NIL series
</td>
<td style="text-align:left;">
JU1793 -> JU2466 (whole genome)
</td>
<td style="text-align:right;">
+0.638 \[+0.569, +0.697\]
</td>
<td style="text-align:right;">
100% \[89, 109\]
</td>
</tr>
<tr>
<td style="text-align:left;">
JU allele swap
</td>
<td style="text-align:left;">
96T vs 96K in the JU1793 background (cost of the JU2466 residue)
</td>
<td style="text-align:right;">
+0.417 \[+0.340, +0.489\]
</td>
<td style="text-align:right;">
46% \[38, 54\]
</td>
</tr>
<tr>
<td style="text-align:left;">
JU allele swap
</td>
<td style="text-align:left;">
96T vs 96K in the JU2466 background (cost of the JU2466 residue)
</td>
<td style="text-align:right;">
+0.139 \[+0.087, +0.198\]
</td>
<td style="text-align:right;">
15% \[10, 22\]
</td>
</tr>
<tr>
<td style="text-align:left;">
JU allele swap
</td>
<td style="text-align:left;">
JU1793 -> JU2466 (whole genome)
</td>
<td style="text-align:right;">
+0.903 \[+0.857, +0.931\]
</td>
<td style="text-align:right;">
100% \[95, 103\]
</td>
</tr>
<tr>
<td style="text-align:left;">
N2 allele swap
</td>
<td style="text-align:left;">
96T vs 96K in the N2 background, 25% pos-1 food
</td>
<td style="text-align:right;">
+0.282 \[+0.221, +0.348\]
</td>
<td style="text-align:right;">
—
</td>
</tr>
<tr>
<td style="text-align:left;">
N2 allele swap
</td>
<td style="text-align:left;">
96T vs 96K in the N2 background, 50% pos-1 food
</td>
<td style="text-align:right;">
+0.026 \[-0.004, +0.068\]
</td>
<td style="text-align:right;">
—
</td>
</tr>
<tr>
<td style="text-align:left;">
N2 allele swap
</td>
<td style="text-align:left;">
96T vs 96K in the N2 background, no pos-1 food
</td>
<td style="text-align:right;">
-0.006 \[-0.034, +0.009\]
</td>
<td style="text-align:right;">
—
</td>
</tr>
</tbody>
</table>

<div class="derived">

**The interval and the residue agree in magnitude.** The 37 kb interval
alone is **32% \[23, 40\]** of the parental span. Residue 96 alone costs
**46% \[38, 54\]** of the span in the JU1793 background and **15% \[10,
22\]** in the JU2466 background — and 32% sits between them, the mean of
the two backgrounds being 30.5%. Every row of the table is stated the
same way, as the cost of carrying the JU2466 allele, so the swap rows
and the NIL rows are directly comparable. So residue 96 is of the right
size to account for the *whole* interval effect, with nothing left over
that requires a second causal variant. That is what Figure S18
independently implies: only two missense differences exist in the entire
37 kb and both are in *sid-2*.

**But this is a match in magnitude, not a statistical identity.** The
swap is **3-fold background-dependent**, so it has no single effect size
to test against, and neither background’s interval overlaps the NIL
estimate cleanly — the JU1793 estimate touches it only at the edge and
the JU2466 estimate misses it.

**Residue 96 matters in the N2 background too, and by a comparable
amount.** The N2 swap has no normalised column because that experiment
has no second parental genotype — one background, two alleles, across a
dose series — so only its raw effect is quotable. At **25% *pos-1* food
it is +0.282 \[0.221, 0.348\]**, which sits between the two JU
backgrounds’ raw effects (`+0.139` and `+0.417`). That is the third
background of Figure 4’s claim, and it lands mid-range rather than at
either extreme.

**The N2 window is open only at 25%.** Across the dose series the effect
is `−0.007` at 0%, **`+0.283` at 25%**, `+0.026` at 50%, `+0.007` at 75%
and `−0.007` at 100%: N2 carrying 96T is already at `0.060` hatching by
50% food, so above 25% there is no room left for the allele to matter
and the contrast collapses into the floor. The 0% row is the negative
control and behaves — no RNAi, no effect. Note this also means the N2
comparison is at a **different dose** from the JU experiments’ nominal
50%, so even its raw effect is not strictly beside theirs.

**The two chromosome III segments are not additive.** The region distal
to the interval does nothing measurable alone — `+0.021` \[−0.008,
0.049\], 3% of span, an interval spanning zero — yet it more than
doubles the interval’s effect when combined: additivity predicts `0.772`
hatched for wSZ176 and the observed value is `0.585`, an excess of
`0.188`. So the interval’s 32% is itself contingent on the chromosome
III background it sits in.

</div>

<div class="aside">

<span class="ch">The cross QTL agrees on direction and rank, and its
magnitude is not convertible</span>

**Direction.** Δfreq = `−0.416` at III:13.78 in the HT115-vs-*pos-1*
contrast. By the sign convention — parent 1 (JU1793) frequency in the
first pool minus the second — a negative value means the JU1793 allele
is **enriched in the *pos-1*-selected pool**, i.e. JU1793 confers
resistance. That is the same direction as both hatching experiments.

**Rank.** III:13.78 carries LOD `140`, against `57` for the next
separated locus in that contrast. By Δfreq alone chrIV:6.42 is
marginally larger (`0.428` against `0.416`).

**Magnitude, not attempted.** An allele-frequency shift in a selected
pool is a selection response. Converting it to a difference in hatched
fraction requires the selection intensity and the number of generations,
and `METHODS.txt` carries both as \[TO FILL\]. No conversion is made
here rather than a fabricated one.

**One coherence worth noting.** The chromosome III right arm accounts
for **64% \[55, 74\]** of the parental span in the NIL series, leaving
roughly a third elsewhere — and the cross independently finds a second
large locus at **chrIV:6.42 Mb, Δfreq `0.428`**, which the NIL series
could not have detected, because every introgression it tests lies on
chromosome III.

</div>

<div class="caveat">

<span class="ch">An asymmetry in the 94 swaps</span>

In the **JU2466** background the N94A swap gives `0.425` hatched against
JU2466’s `0.045` — **+0.38, or 42% of the span, nearly three times the
15% that residue 96 accounts for in the same background.** In the
**JU1793** background it does nothing at all (`0.993` against `0.948`).

This does not revive the glycosylation hypothesis, and the Figure 4
caveats are right that the N94A test of it fails: if the sequon
mattered, removing it should have made JU1793 sensitive, and it did not.
But the effect-size asymmetry is real and is not otherwise recorded —
the largest single-residue effect measured in the JU2466 background is
at 94, not at 96.

</div>

## Figure S11 — the full N2 dose series

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_n2_swap_dose.R`<br> **Supports** Figure
4B’s dose choice · doses 0, 25, 50, 75, 100%

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_n2_swap_dose.png" alt="Embryo hatching for N2 and two edited 96K lines across the full pos-1 dilution series, and the 25% dose alone." width="100%" />
<p class="filecap">
SUPP_FIG_XX_n2_swap_dose
</p>

</div>

<div class="panel">

<span class="pl">A</span> Every genotype at every dose: 0, 25, 50, 75
and 100% *pos-1* RNAi bacteria. N2 carries 96T; wSZ203 and wSZ204 are
independent lines edited to 96K. One plate per strain per dose with
Wilson 95% intervals.

</div>

<div class="panel">

<span class="pl">B</span> The 25% dose alone, as in Figure 4B.

</div>

<table>
<thead>
<tr>
<th style="text-align:right;">
% pos-1
</th>
<th style="text-align:right;">
96T
</th>
<th style="text-align:right;">
96K
</th>
<th style="text-align:right;">
Difference
</th>
<th style="text-align:right;">
Fisher p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.985
</td>
<td style="text-align:right;">
0.992
</td>
<td style="text-align:right;">
-0.006
</td>
<td style="text-align:right;">
0.43
</td>
</tr>
<tr>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.323
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
+0.282
</td>
<td style="text-align:right;">
8.0e-25
</td>
</tr>
<tr>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0.060
</td>
<td style="text-align:right;">
0.034
</td>
<td style="text-align:right;">
+0.026
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
0.034
</td>
<td style="text-align:right;">
0.027
</td>
<td style="text-align:right;">
+0.007
</td>
<td style="text-align:right;">
0.64
</td>
</tr>
<tr>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.007
</td>
<td style="text-align:right;">
-0.007
</td>
<td style="text-align:right;">
0.55
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from
supplemental_data/hatching_assays/n2_allele_swaps_hatching.tsv

</div>

25% is not merely where the difference is largest (+0.282, Fisher
8.0e-25) — it is the only dose that can carry one. Every other dilution
sits at the ceiling or the floor for all three genotypes. **A
single-dose experiment at full strength would have found nothing** —
worth stating, because it is also the likely reason a sub-maximal dose
is needed to see *sid-2* alleles in a resistant background at all.

## Figure S12 — everything held back from Figure 4A

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_sid2_allele_swaps_full.R` · 50% *pos-1*
RNAi<br> **Supports** Figure 4A · 7 constructs, both food conditions

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_allele_swaps_full.png" alt="Three panels: all seven allele-swap constructs on both food conditions, the N-glycosylation sequon map, and the ectodomain coloured by pLDDT." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_allele_swaps_full
</p>

</div>

<div class="panel">

<span class="pl">A</span> All seven constructs, HT115 and *pos-1*.
Control hatching is 97.7–100% for every construct, so the *pos-1*
differences are not a property of the edits.

</div>

<div class="panel">

<span class="pl">B</span> T96 is the threonine of an N94-C95-T96 sequon;
nine N-x-S/T motifs are marked, and the four on the cytoplasmic side
cannot be glycosylated.

</div>

<div class="panel">

<span class="pl">C</span> The ectodomain coloured by pLDDT, with
per-domain means and the dimer ipTM values that limit what the structure
can be used for.

</div>

<div class="caveat">

<span class="ch">Why the glycosylation hypothesis is not in Figure
4</span>

N94A removes the same sequon as T96K while leaving residue 96 alone, and
should phenocopy T96K if loss of the glycan is what matters. It does
not:

<table>
<thead>
<tr>
<th style="text-align:left;">
Background
</th>
<th style="text-align:left;">
Construct
</th>
<th style="text-align:left;">
Motif
</th>
<th style="text-align:right;">
Embryos
</th>
<th style="text-align:right;">
pos-1 hatching
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:left;">
96T (wt)
</td>
<td style="text-align:left;">
NxT
</td>
<td style="text-align:right;">
213
</td>
<td style="text-align:right;">
0.948
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:left;">
94A
</td>
<td style="text-align:left;">
AxT
</td>
<td style="text-align:right;">
267
</td>
<td style="text-align:right;">
0.993
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:left;">
96K
</td>
<td style="text-align:left;">
NxK
</td>
<td style="text-align:right;">
207
</td>
<td style="text-align:right;">
0.531
</td>
</tr>
<tr>
<td style="text-align:left;">
JU2466
</td>
<td style="text-align:left;">
96T
</td>
<td style="text-align:left;">
NxT
</td>
<td style="text-align:right;">
217
</td>
<td style="text-align:right;">
0.184
</td>
</tr>
<tr>
<td style="text-align:left;">
JU2466
</td>
<td style="text-align:left;">
94A
</td>
<td style="text-align:left;">
AxK
</td>
<td style="text-align:right;">
207
</td>
<td style="text-align:right;">
0.425
</td>
</tr>
<tr>
<td style="text-align:left;">
JU2466
</td>
<td style="text-align:left;">
96K (wt)
</td>
<td style="text-align:left;">
NxK
</td>
<td style="text-align:right;">
419
</td>
<td style="text-align:right;">
0.045
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from
supplemental_data/hatching_assays/ju_allele_swaps_hatching.csv

</div>

The AxT construct cannot be glycosylated at N94 and is fully resistant,
so the glycan is not required for resistance; and removing Asn94 gains
more in JU2466 (+0.38) than restoring Thr96 does (+0.14), while gaining
almost nothing in JU1793 (+0.045). Positions 94 and 96 interact. A
description that fits all six numbers is that Lys96 is deleterious and
an unglycosylated Asn94 makes it worse, either relievable independently
— **but that is epistasis between neighbouring residues, not a
mechanism, and each genotype is a single plate.**

</div>

## Figure S13 — the honest negative check

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_sid2_allele_in_panel.R`<br> **Supports**
the candidate, against the mapping panel<br> **n** 230 of 231 phenotyped
strains · site chrIII:13,680,248 C>A

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_allele_in_panel.png" alt="Two panels: the pooled pos-1 phenotype split by sid-2 residue 96 across the mapping panel, and the association scan 300 kb either side of the variant." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_allele_in_panel
</p>

</div>

The T96K allele has no marginal effect across the mapping panel. This is
a negative result and it is the honest check on the candidate.

<div class="panel">

<span class="pl">A</span> The 2023 pooled *pos-1* phenotype split by
*sid-2* residue 96 for every phenotyped strain with a genotype (230 of
231). 96K is at 36% frequency in this set. 96T: `n = 147`, mean −0.0147,
median −0.0220. 96K: `n = 83`, mean −0.0236, median −0.0164. Wilcoxon
`p = 0.99`; Welch `p = 0.16`; `r² = 0.007`. The means differ in the
direction the crosses predict but the medians do not, because the 96T
group carries the resistant tail — so the mean difference is a tail
effect rather than a shift. Open circles mark the four strains whose
allele is known independently: JU1793 and N2 are 96T, JU2466 and XZ1516
are 96K.

</div>

<div class="panel">

<span class="pl">B</span> The association scan 300 kb either side of the
variant. The T96K marker itself reaches `−log10 p = 0.62` (`p = 0.24`),
ranking **18,662 of 64,423** markers on chromosome III, and nothing
within 300 kb clears either threshold. The chromosome III peak for this
trait is 7.7 Mb away, at 5.97 Mb.

</div>

<div class="aside">

<span class="ch">How to say this in the text</span>

This is not evidence against the allele; it is what a
background-dependent effect looks like from a marginal test. The editing
experiments put the T96K contribution at roughly half of the
JU1793–JU2466 difference in one pair of backgrounds, and a 36%-frequency
variant with an effect that context-dependent would not be expected to
surface in a marginal scan of 231 strains. **It is worth stating plainly
that *sid-2* was found by the cross and the NILs, not by the GWAS.**

</div>

## Figure S14 — surface charge of the ectodomain

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_sid2_electrostatics.R`<br> **Supports**
Figure 4C, as hypothesis framing<br> **Scope** residues 21–193 · gut
lumen pH 4.4

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_electrostatics.png" alt="Three panels: the ectodomain coloured by residue class, net charge against pH, and the distribution of Ca distances from residue 96." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_electrostatics
</p>

</div>

A hypothesis-framing figure: it estimates no binding affinity and
contains no docking result.

<div class="panel">

<span class="pl">A</span> The ectodomain coloured by residue class, with
histidine separated because it is the class that titrates between
neutral pH and the acidic gut lumen. T96 in orange, the uptake-critical
residues in white.

</div>

<div class="panel">

<span class="pl">B</span> Net charge against pH, by
Henderson–Hasselbalch over the side chains of residues 21–193 (11 Asp, 7
Glu, 3 His, 4 Cys, 5 Tyr, 9 Lys, 1 Arg; termini omitted because the
ectodomain runs into the transmembrane helix). The domain carries
`−8.4 e` at pH 7.4 but only `−0.2 e` at the gut-lumen pH of 4.4, and its
isoelectric point, 4.38, is essentially the pH it works at. Against that
near-neutral background the +1 e from T96K takes the surface from −0.2 e
to `+0.8 e` and the isoelectric point to 4.50 — a small change in
absolute terms, but at lumenal pH it is the difference between a
slightly negative and a slightly positive face, and **96K is the allele
with efficient uptake**.

</div>

<div class="panel">

<span class="pl">C</span> Cα distances from residue 96 to every other
ectodomain residue, with the four published uptake-critical residues
marked against the null. Three of four are nearer than the median, but
41% of the domain is within 20 Å, so binomial `p = 0.19` and a
permutation test on their mean distance gives `p = 0.30`.

</div>

<div class="caveat">

<span class="ch">This corrects an earlier internal figure</span>

The stage4 analysis in the structure-modelling folder reported the
ectodomain flipping to *strongly* net positive at pH 4.4. It assigned
fixed fractional charges of −0.240 to Asp and −0.334 to Glu, which
correspond to pKa values near 3.4–3.7; with standard pKa values Asp is
still about 76% ionised at pH 4.4 and the domain comes out close to
electroneutral. **The numbers in panel B are the ones to quote.**

</div>

<div class="caveat">

<span class="ch">No potential map is shown, deliberately</span>

The one in the structure-modelling folder is a Coulomb potential with a
distance-dependent dielectric sampled on a plane 55 Å above the protein,
in arbitrary units, with no ionic screening — and its headline panel is
the difference between the T96 and K96 surfaces, necessarily a smooth
positive blob centred on residue 96, since adding +1 e there cannot
produce anything else. The HADDOCK dsRNA docking in the same folder is
also unusable: positive HADDOCK scores (+98 to +200, where favourable is
negative), interaction energies of order 1e4, restraint-violation
energies of 480–1380, 4–13 structures per cluster, and different
91-residue active-residue lists between the two arms being compared.

</div>

# Diagnostics

Not manuscript figures. These settle methodological questions that arose
while assembling the figures above, and they live in
`plots/diagnostics/`. They are tracked so this report reads from a
clone, but unlike the twenty-one they need the Dryad archive to rebuild.

## Leakage in the MIP-seq validation

<div class="meta">

**Script** `scripts/baugh_leakage_vs_similarity.R`<br> **Asks** whether
the strains the deconvolution resolves badly are the genetically similar
ones<br> **Predictor** `baugh_strain_similarity.tsv` — identity-by-state
within the Baugh design matrix, the reference actually solved

</div>

<img src="plots/diagnostics/baugh_leakage_vs_similarity.png" alt="Four panels relating identity-by-state to the NNLS minus MIP-seq discrepancy, with a control panel and a within-pair residual-correlation panel." width="100%" />

This ports the dilution experiment’s question to the one dataset with an
independent measurement of the same material. Relatedness predicts
NNLS–MIP-seq disagreement: ρ `+0.319` over all 98 strains
(`p = 0.0014`), `+0.375` over the 95 whose nearest neighbour is itself
measured here (`p = 0.00018`), surviving with abundance held constant.
Confusable pairs trade signal — the most-confusable 2% of pairs have
median residual correlation `−0.100` against `−0.004` for the rest
(Wilcoxon `p = 0.0079`).

<table>
<thead>
<tr>
<th style="text-align:left;">
Strain
</th>
<th style="text-align:left;">
Closest relative
</th>
<th style="text-align:right;">
IBS
</th>
<th style="text-align:right;">
NNLS slope
</th>
<th style="text-align:right;">
MIP slope
</th>
<th style="text-align:right;">
Ratio
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NIC256
</td>
<td style="text-align:left;">
NIC271
</td>
<td style="text-align:right;">
0.9948
</td>
<td style="text-align:right;">
-1.37e-04
</td>
<td style="text-align:right;">
-4.32e-04
</td>
<td style="text-align:right;">
0.32
</td>
</tr>
<tr>
<td style="text-align:left;">
NIC271
</td>
<td style="text-align:left;">
NIC256
</td>
<td style="text-align:right;">
0.9948
</td>
<td style="text-align:right;">
2.73e-05
</td>
<td style="text-align:right;">
-5.87e-05
</td>
<td style="text-align:right;">
-0.47
</td>
</tr>
<tr>
<td style="text-align:left;">
PS2025
</td>
<td style="text-align:left;">
ECA348
</td>
<td style="text-align:right;">
0.9841
</td>
<td style="text-align:right;">
7.67e-05
</td>
<td style="text-align:right;">
1.98e-04
</td>
<td style="text-align:right;">
0.39
</td>
</tr>
<tr>
<td style="text-align:left;">
JU782
</td>
<td style="text-align:left;">
NIC271
</td>
<td style="text-align:right;">
0.9816
</td>
<td style="text-align:right;">
-5.26e-04
</td>
<td style="text-align:right;">
-6.63e-04
</td>
<td style="text-align:right;">
0.79
</td>
</tr>
<tr>
<td style="text-align:left;">
NIC262
</td>
<td style="text-align:left;">
NIC271
</td>
<td style="text-align:right;">
0.9813
</td>
<td style="text-align:right;">
2.73e-04
</td>
<td style="text-align:right;">
4.88e-04
</td>
<td style="text-align:right;">
0.56
</td>
</tr>
<tr>
<td style="text-align:left;">
CX11264
</td>
<td style="text-align:left;">
CX11262
</td>
<td style="text-align:right;">
0.9729
</td>
<td style="text-align:right;">
-5.20e-04
</td>
<td style="text-align:right;">
-4.88e-04
</td>
<td style="text-align:right;">
1.07
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:left;">
JU2106
</td>
<td style="text-align:right;">
0.9726
</td>
<td style="text-align:right;">
-5.12e-04
</td>
<td style="text-align:right;">
-5.84e-04
</td>
<td style="text-align:right;">
0.88
</td>
</tr>
<tr>
<td style="text-align:left;">
JU2106
</td>
<td style="text-align:left;">
JU1793
</td>
<td style="text-align:right;">
0.9726
</td>
<td style="text-align:right;">
-4.79e-04
</td>
<td style="text-align:right;">
-3.29e-04
</td>
<td style="text-align:right;">
1.46
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from baugh_nnls_with_mipseq.RData + baugh_strain_similarity.tsv

</div>

<div class="aside">

<span class="ch">Does it reach the manuscript? Mostly not — and that is
the finding</span>

Two things protect the manuscript, and they stack.

**The reference here contains only pooled strains.** All 102 columns are
strains that were in the pool, so the set-level leakage that costs the
dilution experiment \~20% of each pure pool is **zero by construction**.
That 20% is not a property of NNLS — it is what happens when the solver
is offered 127 candidates that are not in the pool. The transferable
finding is a **design rule**: restrict the reference to strains actually
pooled.

**The phenotype is a change, and the error differences out.**
Confusability is systematic across timepoints, so a pair splitting mass
in constant proportion gets the level wrong and the change much less so:

| error measure        | ρ vs relatedness |          p |
|:---------------------|-----------------:|-----------:|
| per-sample frequency |         `+0.322` |   `0.0012` |
| fitted slope         |         `+0.130` | **`0.20`** |

Median slope error is `0.086` per mille for the 13 strains above IBS
0.97 against `0.075` for the other 85 — 1.15×, Mann-Whitney `p = 0.48`.
The aggregation says it from the other side too: slopes agree at ρ
`0.974` where per-sample frequencies agree at `0.835`. Every manuscript
phenotype is a difference, so this applies throughout.

**Caveat on that conclusion.** `p = 0.20` at n = 98 is absence of
evidence, and the effect falls from 0.322 to 0.130 rather than to zero.
Read it as “no longer detectable at this n”, not “provably absent”.

**What still stands.** Individual strains at the extreme remain
unreliable even in slopes — NIC256’s NNLS slope is 0.32× the MIP slope
and NIC271’s has the *wrong sign* — so an extreme phenotypic outlier
from that pair deserves a check before it is believed. And a
single-timepoint frequency should not be used as a phenotype: that is
where ρ +0.32 lives, and it is structured along kinship, the same axis a
GWAS kinship matrix models.

</div>

## GWAS interval admission

<div class="meta">

**Script** `scripts/diagnostic_gwas_intervals.R`<br> **Settles** how to
admit a QTL from the association scan so an unsupported lone marker is
excluded and a supported sub-Bonferroni region is kept

</div>

<img src="plots/diagnostics/gwas_interval_diagnostic.png" alt="Chromosome III association scan with both thresholds, zooms on an isolated marker and a supported cluster, and local support against significance genome-wide." width="100%" />

The marker at **5.966 Mb** clears Bonferroni (`8.68`) with **zero**
other eigen-passing markers within 100 kb, out of 628 present. The
cluster at **12.70–12.80 Mb** peaks *below* Bonferroni (`6.31`) with
**14**. A threshold cannot separate them; local support can, and
admission at ≥1 supporting marker excludes the first and keeps the
second, stable up to k = 5.

Corroboration arrived independently: genome-wide, the 14 isolated
markers have median allele frequency `0.082` against `0.394` for the 451
supported ones — the low-frequency signature of spurious association.

### Interval extent against the LD cutoff

<img src="plots/diagnostics/gwas_qtl_intervals_eigen.png" alt="Interval width against LD cutoff per locus, the intervals themselves faceted by cutoff, and the LD profile of the chromosome III peak." width="100%" />

Admission by local support gives **11 loci** at the eigen threshold.
Extent is then LD to the peak marker, and the cutoff matters enormously:

<table>
<thead>
<tr>
<th style="text-align:left;">
locus
</th>
<th style="text-align:right;">
-log10 p
</th>
<th style="text-align:right;">
r² 0.5
</th>
<th style="text-align:right;">
r² 0.6
</th>
<th style="text-align:right;">
r² 0.7
</th>
<th style="text-align:right;">
r² 0.8
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
IV:15.32
</td>
<td style="text-align:right;">
8.84
</td>
<td style="text-align:right;">
1046
</td>
<td style="text-align:right;">
242
</td>
<td style="text-align:right;">
241
</td>
<td style="text-align:right;">
241
</td>
</tr>
<tr>
<td style="text-align:left;">
X:4.88
</td>
<td style="text-align:right;">
7.83
</td>
<td style="text-align:right;">
1461
</td>
<td style="text-align:right;">
1355
</td>
<td style="text-align:right;">
1020
</td>
<td style="text-align:right;">
1020
</td>
</tr>
<tr>
<td style="text-align:left;">
IV:13.41
</td>
<td style="text-align:right;">
7.49
</td>
<td style="text-align:right;">
2033
</td>
<td style="text-align:right;">
2033
</td>
<td style="text-align:right;">
149
</td>
<td style="text-align:right;">
36
</td>
</tr>
<tr>
<td style="text-align:left;">
III:12.72
</td>
<td style="text-align:right;">
6.31
</td>
<td style="text-align:right;">
2347
</td>
<td style="text-align:right;">
2347
</td>
<td style="text-align:right;">
2347
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
X:5.80
</td>
<td style="text-align:right;">
5.62
</td>
<td style="text-align:right;">
3322
</td>
<td style="text-align:right;">
3322
</td>
<td style="text-align:right;">
3322
</td>
<td style="text-align:right;">
3322
</td>
</tr>
<tr>
<td style="text-align:left;">
V:0.61
</td>
<td style="text-align:right;">
5.40
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
373
</td>
<td style="text-align:right;">
201
</td>
<td style="text-align:right;">
179
</td>
</tr>
<tr>
<td style="text-align:left;">
IV:10.53
</td>
<td style="text-align:right;">
5.17
</td>
<td style="text-align:right;">
5885
</td>
<td style="text-align:right;">
5682
</td>
<td style="text-align:right;">
5125
</td>
<td style="text-align:right;">
4411
</td>
</tr>
<tr>
<td style="text-align:left;">
III:3.82
</td>
<td style="text-align:right;">
4.76
</td>
<td style="text-align:right;">
5453
</td>
<td style="text-align:right;">
5453
</td>
<td style="text-align:right;">
5453
</td>
<td style="text-align:right;">
3765
</td>
</tr>
<tr>
<td style="text-align:left;">
III:4.41
</td>
<td style="text-align:right;">
4.76
</td>
<td style="text-align:right;">
4868
</td>
<td style="text-align:right;">
4839
</td>
<td style="text-align:right;">
4279
</td>
<td style="text-align:right;">
3765
</td>
</tr>
<tr>
<td style="text-align:left;">
IV:17.10
</td>
<td style="text-align:right;">
4.70
</td>
<td style="text-align:right;">
2204
</td>
<td style="text-align:right;">
1111
</td>
<td style="text-align:right;">
1111
</td>
<td style="text-align:right;">
858
</td>
</tr>
<tr>
<td style="text-align:left;">
IV:12.08
</td>
<td style="text-align:right;">
4.60
</td>
<td style="text-align:right;">
5939
</td>
<td style="text-align:right;">
5568
</td>
<td style="text-align:right;">
5081
</td>
<td style="text-align:right;">
5081
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from plots/diagnostics/TABLE_gwas_interval_r2_sweep_eigen.tsv —
widths in kb

</div>

**Loci that localise** — interval under 10% of their chromosome — go
`3`, `4`, `5`, `6` of 11 as the cutoff rises 0.5 → 0.8. It tracks signal
strength: median peak −log₁₀p is `6.90` for the localising loci against
`4.76` for the rest, because a marginal peak’s LD partners are
scattered.

<div class="caveat">

<span class="ch">What this implies for the concordance claim, and it is
cutoff-dependent</span>

The chromosome III association locus at 12.718 Mb meets the NIL interval
at r² 0.5, 0.6 **and** 0.7 — all giving a 2,347 kb interval, 17% of the
chromosome — and collapses to **14 kb** at r² 0.8, sitting `0.93` Mb
away. The signal-drop interval (97 kb) is `0.86` Mb away.

Panel C shows why. The peak’s distant LD partners span r²
`0.505`–`0.775`, so any cutoff above 0.78 drops every one of them and
the interval falls back to the local block. The plateau at 0.5–0.7 is
held open by a marker at 13.784 Mb (r² `0.701`) — which is the
chromosome’s **terminal marker**, the same boundary artefact that makes
the cross scan’s peak position uninformative.

</div>

<div class="caveat">

<span class="ch">Where *sid-2* falls in that block, and what it does not
show</span>

Stated explicitly, because it is the question a reader arrives with. The
*sid-2* T96K variant sits at `III:13,680,248` = `13.680` Mb, `0.962` Mb
distal to the association peak. Against the intervals above:

| cutoff   | interval         |    width | contains *sid-2*? |
|:---------|:-----------------|---------:|:------------------|
| r² ≥ 0.5 | 11.436–13.784 Mb | 2,347 kb | **yes**           |
| r² ≥ 0.6 | 11.436–13.784 Mb | 2,347 kb | **yes**           |
| r² ≥ 0.7 | 11.436–13.784 Mb | 2,347 kb | **yes**           |
| r² ≥ 0.8 | 12.718–12.732 Mb |    14 kb | no                |

So at r² ≥ 0.7 the association’s LD block does contain *sid-2*, and the
“megabases apart” objection to the concordance is answered at that
cutoff: the pooled *pos-1* peak, both cross peaks and *sid-2* all fall
inside one block.

**Three limits, all of which belong with the claim.** It is
cutoff-dependent — one step to r² ≥ 0.8 and the interval is 14 kb and
excludes *sid-2*. The plateau that reaches *sid-2* rests on the terminal
marker noted above, so the distal edge is a boundary artefact rather
than a measured LD boundary. And containment is not evidence: the NIL
interval is flat in this scan — 88 markers, best −log₁₀p `0.98`, with
T96K itself at `0.62`, ranking 18,662 of 64,423 on the chromosome.

What this supports is the negative statement, which is the useful one.
The association cannot resolve *within* this block, so *sid-2* lying
inside it is **consistent with** the NIL result rather than in tension
with it. The positive claim — that *sid-2* is the gene — is carried by
the NILs and the allele swaps, not by the scan.

One distinction worth keeping: at r² 0.5 the overlap is **genuine** — 8
markers inside the NIL window reach r² ≈ `0.527` with the peak. At 0.6
and 0.7 no qualifying marker falls inside it, and the interval spans it
only as a min–max hull. So the overlap is real only under the most
permissive cutoff.

**The defensible statement** is that the association signal and the NIL
interval are about 0.9 Mb apart, and that they coincide only under an
interval definition wide enough to cover a sixth of chromosome III. That
reads as a distinct locus, not concordance.

</div>

The supported cluster sits at 12.70–12.80 Mb. The NIL interval is
13.658–13.695 Mb, and **13.5–13.9 Mb has nothing above the eigen line at
all**. Defining intervals rigorously may therefore show the association
signal is a *distinct locus* rather than confirming concordance with the
crosses. Worth deciding how to present before the intervals are drawn.
:::

## Coverage against reference size

<div class="meta">

**Script** `scripts/diagnostic_reference_size.R`<br> **Asks** how slope
recovery degrades as coverage falls, and how much worse it gets as the
reference admits candidates that are not in the pool

</div>

<img src="plots/diagnostics/reference_size_vs_coverage.png" alt="Slope agreement, RMSE, discrepancy and relatedness correlation against fraction of full-depth reads, one line per reference size." width="100%" />

**A larger reference costs per-strain accuracy, and almost nothing in
the slopes.** The 102 pooled strains are always present; larger
references add candidates from the other 438 CeNDR isotypes. Depth is
the fraction of full-depth reads retained, by binomial thinning of the
archived alt counts (exact read thinning); full depth is 2.576 alt reads
per marker, so absolute coverage is not recoverable from this archive
and the axis is relative.

The reconstruction is validated rather than assumed: rebuilding the
reference from CeNDR genotypes and applying the recovered flip mask
reproduces the archive’s design matrix with **0 mismatching cells of
126,184,812**, and the deconvolution at reference 102 / full depth
reproduces the archived frequencies with a **maximum difference of
exactly 0** across 2,346 strain-samples.

<table>
<thead>
<tr>
<th style="text-align:left;">
Depth (fraction)
</th>
<th style="text-align:right;">
Alt reads/marker
</th>
<th style="text-align:right;">
R=102
</th>
<th style="text-align:right;">
R=150
</th>
<th style="text-align:right;">
R=250
</th>
<th style="text-align:right;">
R=400
</th>
<th style="text-align:right;">
R=540
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1/64
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.936
</td>
<td style="text-align:right;">
0.938
</td>
<td style="text-align:right;">
0.948
</td>
<td style="text-align:right;">
0.955
</td>
</tr>
<tr>
<td style="text-align:left;">
1/32
</td>
<td style="text-align:right;">
0.080
</td>
<td style="text-align:right;">
0.940
</td>
<td style="text-align:right;">
0.946
</td>
<td style="text-align:right;">
0.945
</td>
<td style="text-align:right;">
0.941
</td>
<td style="text-align:right;">
0.937
</td>
</tr>
<tr>
<td style="text-align:left;">
1/16
</td>
<td style="text-align:right;">
0.161
</td>
<td style="text-align:right;">
0.952
</td>
<td style="text-align:right;">
0.955
</td>
<td style="text-align:right;">
0.955
</td>
<td style="text-align:right;">
0.954
</td>
<td style="text-align:right;">
0.952
</td>
</tr>
<tr>
<td style="text-align:left;">
1/8
</td>
<td style="text-align:right;">
0.322
</td>
<td style="text-align:right;">
0.970
</td>
<td style="text-align:right;">
0.971
</td>
<td style="text-align:right;">
0.971
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
0.966
</td>
</tr>
<tr>
<td style="text-align:left;">
1/4
</td>
<td style="text-align:right;">
0.644
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.975
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0.972
</td>
</tr>
<tr>
<td style="text-align:left;">
1/2
</td>
<td style="text-align:right;">
1.288
</td>
<td style="text-align:right;">
0.976
</td>
<td style="text-align:right;">
0.976
</td>
<td style="text-align:right;">
0.974
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.974
</td>
</tr>
<tr>
<td style="text-align:left;">
1/1
</td>
<td style="text-align:right;">
2.576
</td>
<td style="text-align:right;">
0.974
</td>
<td style="text-align:right;">
0.976
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.974
</td>
<td style="text-align:right;">
0.975
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from plots/diagnostics/TABLE_reference_size_vs_coverage.tsv

</div>

<div class="aside">

<span class="ch">Three things the surface says</span>

**Slope recovery is flat in reference size.** At full depth ρ goes
`0.974` at 102 candidates to `0.975` at 540; at 1/64 depth it is `0.937`
against `0.955` — if anything slightly *better* with the larger
reference. So a reference carrying five times as many absent candidates
costs essentially nothing in the quantity the manuscript uses.

**Per-strain accuracy does pay.** Median \|NNLS − MIP\| rises from
`2.30` to `2.85` per mille at full depth, a 24% increase, and from
`2.92` to `3.54` at 1/64 depth. The cost is real but it lands on
individual frequencies, not on fitted slopes — the same split as the
leakage result above.

**Relatedness bites harder as coverage falls, and *less* as the
reference grows.** ρ(nn_ibs, discrepancy) at reference 102 runs `0.322`
at full depth up to `0.504` at 1/64 — confusability needs reads to
resolve, as predicted. But at fixed depth it *falls* with reference size
(`0.322` → `0.178` at full depth), because when every strain has a close
relative among 540 candidates, nearest-neighbour identity stops
discriminating. That is a property of the predictor losing range, not of
the deconvolution improving.

</div>

## Off-pool leakage against coverage

**Script** `scripts/diagnostic_downsample_leakage.R`<br> **Figure**
`plots/diagnostics/downsample_leakage.png`

The section above measures *discrepancy* — how far a pool member’s
estimated frequency sits from its MIP-seq measurement. This one measures
*leakage*: frequency handed to reference columns that are not in the
pool at all. They are different failure modes, and only the second is
what a real experiment risks when it cannot name its own members. At 102
candidates leakage is zero by construction, because every column is a
pool member; at 540 there are 438 columns that should carry nothing.

No deconvolution is rerun. The reference-size diagnostic saved every
fitted frequency for all 540 columns and its own tables then dropped the
absent ones in a join against the MIP-seq measurements, so these numbers
come out of the same fit and cannot drift from what is reported above.

<img src="plots/diagnostics/downsample_leakage.png" alt="Off-pool leakage against sequencing depth by reference size, frequency absorbed against identity-by-state to the closest pool member, the correlation between the two against depth, and slope agreement with MIP-seq as fitted versus restricted to the true pool." width="100%" />

<table>
<caption>
Median across samples of the pool frequency assigned to strains absent
from the pool.
</caption>
<thead>
<tr>
<th style="text-align:left;">
Depth
</th>
<th style="text-align:right;">
R=150
</th>
<th style="text-align:right;">
R=250
</th>
<th style="text-align:right;">
R=400
</th>
<th style="text-align:right;">
R=540
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1/64
</td>
<td style="text-align:right;">
5.03
</td>
<td style="text-align:right;">
10.99
</td>
<td style="text-align:right;">
16.67
</td>
<td style="text-align:right;">
20.11
</td>
</tr>
<tr>
<td style="text-align:left;">
1/32
</td>
<td style="text-align:right;">
4.36
</td>
<td style="text-align:right;">
8.65
</td>
<td style="text-align:right;">
14.64
</td>
<td style="text-align:right;">
17.77
</td>
</tr>
<tr>
<td style="text-align:left;">
1/16
</td>
<td style="text-align:right;">
4.05
</td>
<td style="text-align:right;">
8.49
</td>
<td style="text-align:right;">
13.00
</td>
<td style="text-align:right;">
15.96
</td>
</tr>
<tr>
<td style="text-align:left;">
1/8
</td>
<td style="text-align:right;">
3.96
</td>
<td style="text-align:right;">
8.24
</td>
<td style="text-align:right;">
12.11
</td>
<td style="text-align:right;">
15.24
</td>
</tr>
<tr>
<td style="text-align:left;">
1/4
</td>
<td style="text-align:right;">
3.89
</td>
<td style="text-align:right;">
7.78
</td>
<td style="text-align:right;">
12.14
</td>
<td style="text-align:right;">
15.17
</td>
</tr>
<tr>
<td style="text-align:left;">
1/2
</td>
<td style="text-align:right;">
3.78
</td>
<td style="text-align:right;">
7.54
</td>
<td style="text-align:right;">
12.01
</td>
<td style="text-align:right;">
15.27
</td>
</tr>
<tr>
<td style="text-align:left;">
1/1
</td>
<td style="text-align:right;">
3.62
</td>
<td style="text-align:right;">
7.47
</td>
<td style="text-align:right;">
11.81
</td>
<td style="text-align:right;">
14.87
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from plots/diagnostics/TABLE_downsample_leakage.tsv and
TABLE_downsample_leakage_slopes.tsv

</div>

<div class="aside">

<span class="ch">Both things are true at once</span>

**Leakage is large, and it scales with the count of absent candidates.**
At the full 540-strain reference and full depth, `14.9%` of the pool
goes to strains that are not in it; total pool frequency retained is
`0.8505`. Across reference sizes it runs `3.6%`, `7.5%`, `11.8%`,
`14.9%` for 48, 148, 298 and 438 absent candidates — about `0.34` per
mille each throughout. Nothing subtler than the number of extra columns
is needed to predict it.

**Coverage is almost irrelevant to it.** Over a 64-fold depth reduction
leakage at 540 candidates rises only from `14.9%` to `20.1%`, most of
that in the last two halvings. Low coverage is not what makes an unknown
membership list expensive, which was the thing worth checking and is not
what I would have guessed.

**It is not spread evenly, and individual strains are displaced
outright.** Among strain-samples above 0.5% of the pool the median
retained fraction is `0.898`, but the 5–95% range is `0.396`–`1.04`. The
cleanest case: ECA36 falls from `4.17` per mille at 102 candidates to
`0.00` at 540, while JU3226 — absent from the pool, IBS `0.9898` to
ECA36 — absorbs `4.39`. One strain is swapped for its look-alike almost
exactly.

**And none of it reaches the phenotype.** Spearman agreement between
NNLS growth slopes and MIP-seq slopes is `0.975` at 540 candidates
against `0.974` at 102, with RMSE `0.0001` for both. Restricting the
reference to the true pool and renormalising — what an experiment with a
known membership list gets — moves the third decimal at most. The cost
of not knowing the membership, as loss in ρ against the 102-candidate
reference, is between `-0.001` and `+0.003` at every depth above 1/64.

So the answer to “is this a big deal” is: not for anything the
manuscript claims. Per-strain identity degrades badly at the extreme of
relatedness, and the slope phenotype behind Figure 1 and the pooled GWAS
does not notice, because displacement is consistent across the samples
of a replicate and so shifts a strain’s whole trajectory rather than its
trend.

</div>

<div class="aside">

<span class="ch">Relatedness grades the magnitude, not the
membership</span>

Panel C plots a whole-set rank correlation, and read alone it is
misleading — it falls from `+0.304` at 150 candidates to `+0.051` at 540
and looks like a relationship dissolving. It is not. NNLS is
non-negative, so its solutions are sparse: at 540 candidates and full
depth `56.8%` of absent columns absorb *exactly* zero, and those ties
flatten any Spearman.

Split in two, it resolves. IBS does **not** predict which candidates
enter the solution’s support at all — logistic P(absorbs \> 0) on IBS
has slope `-1.13`, `p = 0.71`, and the rank correlation with the
indicator is `-0.021`. Among the candidates that do absorb something,
IBS predicts how much, and does so consistently at every reference size
and depth: `+0.46`, `+0.39`, `+0.34`, `+0.35` at 150, 250, 400 and 540
at full depth, and between `+0.28` and `+0.49` across the whole grid.
The whole-set decline is entirely the zero fraction climbing from `25%`
to `57%`.

Panel B is the same fact in the form that matters: mean frequency
absorbed rises `9.6`-fold across IBS deciles, from `0.081` per mille in
the least related tenth of absent candidates to `0.785` in the most
related. Decile *means* are plotted rather than medians for exactly the
reason above — with 57% exact zeros every decile median is 0, which
would draw a flat line through data that has a ten-fold gradient in it.

One incidental result worth keeping: sparsity increases with depth. The
zero fraction at 540 candidates goes from `20.5%` at 1/64 depth to
`56.8%` at full depth, so shallow sequencing does not just leak more in
total — it spreads the leak across more wrong strains.

</div>

<table>
<caption>
Full depth. The whole-set correlation tracks the zero fraction; the
among-absorbers correlation does not move.
</caption>
<thead>
<tr>
<th style="text-align:right;">
Reference
</th>
<th style="text-align:right;">
Absent candidates
</th>
<th style="text-align:right;">
Absorbing exactly 0
</th>
<th style="text-align:right;">
ρ, whole set
</th>
<th style="text-align:right;">
ρ, among absorbers
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
150
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
25.0%
</td>
<td style="text-align:right;">
+0.304
</td>
<td style="text-align:right;">
+0.464
</td>
</tr>
<tr>
<td style="text-align:right;">
250
</td>
<td style="text-align:right;">
148
</td>
<td style="text-align:right;">
41.9%
</td>
<td style="text-align:right;">
+0.136
</td>
<td style="text-align:right;">
+0.392
</td>
</tr>
<tr>
<td style="text-align:right;">
400
</td>
<td style="text-align:right;">
298
</td>
<td style="text-align:right;">
47.7%
</td>
<td style="text-align:right;">
+0.126
</td>
<td style="text-align:right;">
+0.340
</td>
</tr>
<tr>
<td style="text-align:right;">
540
</td>
<td style="text-align:right;">
438
</td>
<td style="text-align:right;">
56.8%
</td>
<td style="text-align:right;">
+0.051
</td>
<td style="text-align:right;">
+0.353
</td>
</tr>
</tbody>
</table>

<div class="caveat">

<span class="ch">What is still not established</span>

Which reference columns the solver admits into its support in the first
place. That is a property of the collinearity of the whole design
matrix, not of any one pairwise distance, and nothing measured here
explains it — a strain’s IBS to its closest pool member carries no
information about whether it absorbs anything (`p = 0.71`). The claim
this section supports is about how much an admitted candidate takes,
plus the wholesale displacement seen in the near-duplicate pairs above
IBS ≈ 0.985. It does not support a story about which strains the
deconvolution will choose to be wrong about.

</div>

## Gene content at the four strongest mig-6 QTL

**Scripts** `scripts/make_mig6_locus_tables.R`,
`scripts/diagnostic_mig6_locus_genes.R`<br> **Figures**
`plots/diagnostics/DIAG_mig6_locus_*.png`<br> **Table**
`plots/diagnostics/TABLE_mig6_locus_census.tsv`

The census Figure S18 runs on the 37 kb NIL interval, applied to every
independent HT115-vs-*mig-6* cross QTL above LOD 500, in a 100 kb window
on each peak. Colour is the parent carrying the alternate allele: parent
1 (N2, JU1793) pink, parent 2 (XZ1516, JU2466) green, the same p1/p2
ordering the cross allele-frequency tables use.

Eight peaks clear LOD 500, but four are `peak.rank > 1` with
`separated = FALSE` — shoulders of one sweep rather than independent
QTL, by the trough test in `cross_qtl_full_summary.R` — and a window on
a shoulder is a window on the same locus twice. The four independent
loci are the ones below.

<table>
<thead>
<tr>
<th style="text-align:left;">
Locus
</th>
<th style="text-align:right;">
LOD
</th>
<th style="text-align:right;">
Δfreq
</th>
<th style="text-align:right;">
Coding genes
</th>
<th style="text-align:right;">
Sites
</th>
<th style="text-align:right;">
Differ
</th>
<th style="text-align:right;">
Protein-altering
</th>
<th style="text-align:right;">
HIGH
</th>
<th style="text-align:right;">
No-call p1
</th>
<th style="text-align:right;">
No-call p2
</th>
<th style="text-align:right;">
Divergent bp
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31 Mb
</td>
<td style="text-align:right;">
940
</td>
<td style="text-align:right;">
-0.320
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
2664
</td>
<td style="text-align:right;">
229
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
193
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76 Mb
</td>
<td style="text-align:right;">
806
</td>
<td style="text-align:right;">
-0.384
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
1876
</td>
<td style="text-align:right;">
217
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15 Mb
</td>
<td style="text-align:right;">
798
</td>
<td style="text-align:right;">
-0.836
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2565
</td>
<td style="text-align:right;">
134
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:1.88 Mb
</td>
<td style="text-align:right;">
728
</td>
<td style="text-align:right;">
0.410
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
4464
</td>
<td style="text-align:right;">
235
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
985
</td>
<td style="text-align:right;">
133
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

<div class="derived">

**No HIGH-impact difference at any of the four loci.** No stop gained or
lost, no frameshift, no splice-acceptor or splice-donor change, no start
lost — across 86 protein-coding genes in 400 kb. The protein-altering
differences are all missense or inframe indels: 21 at V:10.76, 12 at
III:13.31, 11 at I:1.88 and 7 at X:7.15.

**None of these windows resolves to a candidate the way the NIL interval
does.** That interval held 2 protein-altering differences in 37 kb, both
in one gene. These hold 7–21 spread across many genes, which is what an
unfine-mapped peak looks like — the contrast is the point, not a defect
in these panels.

**The chromosome III *mig-6* peak is not the *sid-2* locus.** It sits at
13.31 Mb, and *sid-2* is at 13.679–13.682 Mb, so a 100 kb window on this
peak excludes *sid-2* by roughly 270 kb. The genes it does contain are a
different set.

**Two of the JU X:7.15 hits are worth a second look** on identity rather
than statistics: *bar-1* (β-catenin, Wnt signalling) at 209A→209T and
*ist-1* (insulin-receptor substrate) at 339T→339K. Neither is nominated
here as a candidate — this panel ranks nothing — but they are the two
whose function would make a reader stop.

</div>

<div class="caveat">

<span class="ch">Divergent regions, and what actually limits these
censuses</span>

**No divergent region overlaps any of the four windows** for either
parent of the relevant cross. The panels shade them where they occur, so
the shading is absent here rather than omitted. That was not the
expected result.

**What limits the census instead is per-parent missingness, and it is
concentrated in N2.** At I:1.88, 985 of 4,464 sites (22%) are a no-call
in **N2** against 133 in XZ1516; at III:13.31 it is 193 against 30. At
X:7.15 both JU parents are called at every one of 2,565 sites. So the
absence of HIGH-impact differences is near-airtight on chromosome X and
weakest on chromosome I, and the divergent-region track — the instrument
that would normally flag this — says nothing about it.

**Colour is near-degenerate for the N2 cross.** N2 is the reference
genome, so at almost every site where the parents differ it is XZ1516
that carries the alternate allele, and the pink key is essentially
unused in those three panels. It is informative in the JU cross, where
neither parent is the reference: two of the seven X:7.15 hits are JU1793
and five are JU2466.

</div>

<img src="plots/diagnostics/DIAG_mig6_locus_III_13-31.png" alt="" width="100%" />

<img src="plots/diagnostics/DIAG_mig6_locus_V_10-76.png" alt="" width="100%" />

<img src="plots/diagnostics/DIAG_mig6_locus_X_7-15.png" alt="" width="100%" />

<img src="plots/diagnostics/DIAG_mig6_locus_I_1-88.png" alt="" width="100%" />

## Figure S19 — SID-2 across two species, with elegans variation on top

<div class="meta">

**Scripts** `scripts/make_sid2_alignment_tables.py`,
`scripts/SUPP_FIG_XX_sid2_briggsae_alignment.R`<br> **Supports** Figure
4, by asking which SID-2 residues a second species conserves<br>
**Scope** *C. elegans* SID-2 (311 aa) against *C. briggsae* CBR-SID-2
(314 aa) · 17 population missense variants

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_briggsae_alignment.png" alt="Two panels: a per-residue conservation strip comparing C. elegans and C. briggsae SID-2 with the membrane topology beneath and annotated residues marked above, and a lollipop plot of the 17 missense variants segregating in the C. elegans population by allele frequency." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_briggsae_alignment
</p>

</div>

*C. briggsae* is insensitive to environmental RNAi, and a *C. elegans*
*sid-2* transgene confers sensitivity on it (Winston et al. 2007). The
two species therefore bracket a functional difference that *sid-2* alone
is sufficient to explain, which makes their SID-2 proteins the natural
comparison for asking which residues the protein needs.

<div class="derived">

**The two proteins are only 47.3% identical** — 140 identical, 156
different and 15 gapped over 311 *C. elegans* positions — so
conservation here is informative precisely because it is rare.

**All three histidines implicated in dsRNA uptake differ in *C.
briggsae*: H32→R, H168→S, H175→R.** Two of the three go to **arginine**,
which is exactly the substitution McEwan et al. 2012 made. That is why
the panel labels them *implicated* rather than critical: **their triple
His→Arg mutant internalised more dsRNA than wild type, not less**, so
the *briggsae* state at these positions cannot by itself explain why
*briggsae* fails to take dsRNA up. The observation is striking and it
cuts against the simple reading.

**N94 and T96 are both identical in *C. briggsae*.** The focal residue
of Figure 4 is the one thing in this neighbourhood that a species 47%
identical still conserves — and T96 is the conserved state, while the
*C. elegans* population carries 96K at **46%** allele frequency. So K is
the derived, common variant and T is the ancestral one, which is the
opposite of how a rare loss-of-function allele would look.

**Of the 17 missense variants segregating in the population, only 2 sit
at a position *C. briggsae* conserves** — and T96K is one of them. The
other 15 fall at positions already divergent between the species.

**XZ1516 differs from N2 at 7 of the 8 curated sites** (78, 96, 141,
144, 151, 153, 209; only residue 5 is shared), so it carries a strongly
diverged *sid-2* haplotype rather than a single variant of interest.

</div>

<div class="caveat">

<span class="ch">What the alignment cannot carry</span>

**At 47% identity the alignment has real uncertainty**, particularly
across the low-complexity stretches, and a single global
Needleman–Wunsch alignment (BLOSUM62, gap −11/−1) is one hypothesis
about correspondence rather than a fact. The staged table records the
aligned *briggsae* residue and a gap flag per position so a reader can
see where the comparison is weak instead of trusting it uniformly.
Positions 32, 34, 94, 96, 168 and 175 all sit in well-aligned blocks, so
the conclusions above do not rest on the ambiguous regions.

**Two species is not a conservation analysis.** Identity or difference
against one outgroup says nothing about the rate at a site. A proper
test would need an alignment across the *Caenorhabditis* genus, which
this figure does not attempt.

**Sequences are UniProt G5EEV9 and A8XSB8**, fetched rather than derived
here, and the *briggsae* entry is unreviewed (PE=4, predicted). Its gene
model has not been checked against *briggsae* RNA-seq.

</div>

## eQTL and parental expression at the censused loci

**Script** `scripts/candidate_eqtl_expression.R`<br> **Table**
`plots/diagnostics/TABLE_candidate_eqtl.tsv`<br> **Source** the
207-isolate expression matrix, its eQTL table and feature table (outside
the repository)

The coding-variant censuses above ask whether a gene’s protein differs
between the parents. This asks the complementary question — whether its
*expression* does, and whether it has a mapped eQTL in the 207-isolate
study. All four parents of the two crosses (N2, XZ1516, JU1793, JU2466)
are among the 207, so the comparison is direct rather than inferred.

Of **75 coding genes** across the five censused windows, **9 carry a
mapped eQTL** and **5 of those are local**. Expression differences are
reported as a z-score against the 207-strain distribution: how unusual
the parental gap is, not whether it is significant.

<table>
<thead>
<tr>
<th style="text-align:left;">
Locus
</th>
<th style="text-align:left;">
Gene
</th>
<th style="text-align:right;">
p1
</th>
<th style="text-align:right;">
p2
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
H²
</th>
<th style="text-align:right;">
eQTL
</th>
<th style="text-align:right;">
var. exp.
</th>
<th style="text-align:right;">
coding Δ
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15
</td>
<td style="text-align:left;">
ist-1
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-3.45
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15
</td>
<td style="text-align:left;">
C55B6.1
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
3.19
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15
</td>
<td style="text-align:left;">
C54D1.7
</td>
<td style="text-align:right;">
-1.00
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:right;">
-2.86
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15
</td>
<td style="text-align:left;">
alh-10
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:right;">
3.72
</td>
<td style="text-align:right;">
-2.69
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
JU1793xJU2466 X:7.15
</td>
<td style="text-align:left;">
clec-86
</td>
<td style="text-align:right;">
1.55
</td>
<td style="text-align:right;">
3.02
</td>
<td style="text-align:right;">
-2.07
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:1.88
</td>
<td style="text-align:left;">
tub-2
</td>
<td style="text-align:right;">
3.16
</td>
<td style="text-align:right;">
4.19
</td>
<td style="text-align:right;">
-3.35
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
Distant eQTL
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:1.88
</td>
<td style="text-align:left;">
egl-30
</td>
<td style="text-align:right;">
3.61
</td>
<td style="text-align:right;">
4.79
</td>
<td style="text-align:right;">
-2.16
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
Distant eQTL/Local eQTL
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:1.88
</td>
<td style="text-align:left;">
scm-1
</td>
<td style="text-align:right;">
5.51
</td>
<td style="text-align:right;">
5.78
</td>
<td style="text-align:right;">
-2.06
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
Local eQTL
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 I:1.88
</td>
<td style="text-align:left;">
tag-96
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-1.92
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
Local eQTL
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:left;">
Y43F4B.10
</td>
<td style="text-align:right;">
-0.90
</td>
<td style="text-align:right;">
4.30
</td>
<td style="text-align:right;">
-4.13
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:left;">
cpf-2
</td>
<td style="text-align:right;">
5.81
</td>
<td style="text-align:right;">
4.62
</td>
<td style="text-align:right;">
3.75
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:left;">
dph-7
</td>
<td style="text-align:right;">
-1.00
</td>
<td style="text-align:right;">
1.74
</td>
<td style="text-align:right;">
-2.95
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
Local eQTL
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:left;">
dro-1
</td>
<td style="text-align:right;">
6.26
</td>
<td style="text-align:right;">
5.22
</td>
<td style="text-align:right;">
2.48
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 III:13.31
</td>
<td style="text-align:left;">
F56A8.3
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:right;">
3.38
</td>
<td style="text-align:right;">
-2.04
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
pas-2
</td>
<td style="text-align:right;">
6.91
</td>
<td style="text-align:right;">
6.20
</td>
<td style="text-align:right;">
3.62
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
Distant eQTL
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
secs-1
</td>
<td style="text-align:right;">
-1.00
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
-3.30
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
F28H7.2
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-2.46
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
D1054.8
</td>
<td style="text-align:right;">
3.76
</td>
<td style="text-align:right;">
5.07
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
—
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
D1054.1
</td>
<td style="text-align:right;">
3.20
</td>
<td style="text-align:right;">
2.37
</td>
<td style="text-align:right;">
2.01
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
Local eQTL
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N2xXZ1516 V:10.76
</td>
<td style="text-align:left;">
D1054.9
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
Distant eQTL
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
NIL interval III:13.66-13.70
</td>
<td style="text-align:left;">
dyf-2
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-1.27
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
Distant eQTL
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

<div class="derived">

**This kills the two candidates that needed expression to work.**
`set-25` (H3K9 methyltransferase, nuclear RNAi) and `lam-2` (laminin γ,
and it sits on the JU chromosome X peak) were nominated on position and
function precisely because neither carries a protein-altering difference
— so expression was the only mechanism left. Neither differs between the
parents (`z` = +0.29 and −0.71) and neither has an eQTL. That was the
testable form of the hypothesis and it fails.

**`ist-1` becomes the strongest single candidate at X:7.15.** It is the
only gene in the set with both a coding difference and a large
expression difference — `z` = **−3.45**, JU1793 more than three
population SDs below JU2466 — and it lies 11–24 kb from the peak.

**`sid-2` is a positive control for the method, and it behaves.** `z` =
+0.36 with no eQTL: no expression difference between JU1793 and JU2466
at all. That is what it should look like, because the *sid-2* effect is
coding — T96K, confirmed by the allele swap — so this analysis is not
manufacturing signal where the answer is already known.

**`mtm-6` and `wdr-5.3` keep their coding evidence and lose their
regulatory story.** `wdr-5.3` matters here: its expression is highly
heritable across the 207 (H² = 0.70) and the parents still do not
differ, which makes that a meaningful null rather than a noisy one.

**Genes the coding census had not nominated** now come forward on
expression: `Y43F4B.10` (`z` = −4.13, the largest in any window),
`cpf-2` (+3.75) and `dph-7` (−2.95, local eQTL) at III:13.31; `pas-2`
(+3.62), `secs-1` (−3.30) and `D1054.1` — whose local eQTL explains
**69%** of expression variance, the strongest cis-eQTL in any window —
at V:10.76; `tub-2` (−3.35), `egl-30` and `scm-1` at I:1.88.

</div>

<div class="caveat">

<span class="ch">What this comparison cannot say</span>

**One expression value per strain**, so there is no within-strain
replication and no p-value on a parental difference. The z-score
describes how unusual the gap is against the 207-strain spread, nothing
more.

**The expression study is whole-animal, one stage, one condition.** A
gene whose parental difference exists only in the relevant tissue, or
only on RNAi food, would appear null here. The absence of a difference
is therefore weak evidence against a regulatory mechanism, not proof of
its absence.

**An eQTL in the 207 panel is a statement about wild variation in
general**, not about these two strains. A gene can differ between the
parents without having a mappable eQTL, and can have a strong eQTL while
these two parents happen to share an allele.

</div>

## Figure S15 — where T96’s pocket sits in the charge distribution

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_sid2_local_charge.R`<br> **Supports**
Figure 4C, the quantitative half<br> **Scope** residues 21–188 · 12 Å
neighbourhoods · pH 4.4

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_local_charge.png" alt="Histogram of local net charge across the ectodomain, with the three uptake histidines marked and vertical lines at T96 and T96K." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_local_charge
</p>

</div>

Figure 4C shows *where* the positive pocket is. This shows *how unusual*
it is, which is the part that can be argued with.

<div class="aside">

<span class="ch">The claim, and its size</span>

T96’s neighbourhood is `+1.24` e, the **82nd percentile** of the 168
ectodomain residues (median `0.00`); T96K takes it to `+2.24` e, the
**98th**. The pocket is made by two lysines, K93 at `6.6` Å and K132 at
`6.8` Å (`4.5` and `4.4` Å nearest heavy atom).

Two things stop this being over-read. `81` of the 168 residues have
positive local charge at this pH, so a positive pocket by itself is
unremarkable — the percentile is the claim, not the sign. And the domain
is net **acidic** overall, so this is a basic pocket in an acidic
domain, not a polybasic surface of the kind SID-1 uses.

The three uptake histidines are marked for context: H32 at `−0.15`, H168
at `+0.85`, H175 at `+1.46` e. H175 — the one `37.8` Å from T96 and so
absent from Figure 4C’s field — sits in the most positive environment of
the three.

</div>

## Figure S16 — model confidence, and the proximity null

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_sid2_model_confidence.R`<br>
**Supports** Figure 4C, defensively<br> **Scope** residues 21–188 ·
AlphaFold3 pLDDT

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_sid2_model_confidence.png" alt="Ectodomain cartoon coloured by pLDDT, the T96 zoom in the same colouring, and the distribution of Ca distances from T96 with the annotated residues marked." width="100%" />
<p class="filecap">
SUPP_FIG_XX_sid2_model_confidence
</p>

</div>

Figure 4C is coloured by charge and carries no confidence encoding, so
this is where the model’s quality where it matters can be checked, and
where the proximity null is shown rather than argued.

<div class="aside">

<span class="ch">Two defensive points</span>

**The model is sound where the claim is made.** `117` of the 168
ectodomain residues reach pLDDT ≥ 70, and T96 itself is at `79.3` —
confident, though at the edge of the well-modelled core rather than deep
inside it. The low-confidence region is the lumenal cap above T96, and
nothing rests on the cap.

**The proximity argument is shown failing, not omitted.** `41.9%` of the
ectodomain lies within 20 Å of T96, and two of the three uptake
histidines fall inside that radius — binomial `p = 0.38`. The zoom keeps
the two evidence classes apart, the histidines against the *qt13* allele
D34, which is the distinction the released panel collapsed.
`METHODS.txt` states that no proximity claim is made; this figure is why
that can be stated rather than asserted.

</div>

## Figure S17 — Figure 2 without the cross QTL

<div class="meta">

**Script** `scripts/Figure2.R`<br> **Supports** Figure 2, as its
control<br> **Scope** the same two pooled GWAS scans · cross QTL tracks
removed

</div>

<div class="plate">

<img src="plots/Figure2_no_cross_qtl.png" alt="The same mirrored Manhattan as Figure 2 - mig-6 above the axis, pos-1 below - with the cross QTL arrow tracks, their labels and the cross legend removed." width="100%" />
<p class="filecap">
Figure2_no_cross_qtl
</p>

</div>

Figure 2 draws nine cross QTL as arrows over a mirrored Manhattan, which
asks the reader to judge two sources of evidence at once. This is the
same panel with the arrows, their labels and the cross legend removed,
written by the same script from the same data in the same call —
`build(LOD_MIN, "Figure2_no_cross_qtl", arrows = FALSE)` — so nothing
but the QTL layer differs.

<div class="aside">

<span class="ch">What it is for</span>

**It separates the claim from the annotation.** The arrows carry
position only, and the association peaks beneath them are the
independent measurement. Seeing the Manhattan alone is how a reader
checks that a peak is present where an arrow points, rather than reading
the arrow as the evidence.

**It is a control, not an alternative.** Figure 2 is the figure the
argument uses. This one exists so the overlay can be removed and the
underlying scan inspected, which is a question a reviewer is entitled to
ask and which no amount of caption prose answers as directly.

</div>

## Figure S18 — what the 37 kb interval contains

<div class="meta">

**Script** `scripts/SUPP_FIG_XX_nil_interval_genes.R`<br> **Supports**
Figure 3, by saying what is inside the interval it resolves<br>
**Scope** chromosome III 13.6577–13.6950 Mb · 12 genes · 27 parental
differences

</div>

<div class="plate">

<img src="plots/SUPP_FIG_XX_nil_interval_genes.png" alt="Gene models across the 37 kb interval on chromosome III, with sites where JU1793 and JU2466 differ drawn above them: two lollipops for the protein-altering missense variants in sid-2, and a rug of the remaining differences." width="100%" />
<p class="filecap">
SUPP_FIG_XX_nil_interval_genes
</p>

</div>

Figure 3 resolves the QTL to 13.6577–13.6950 Mb but does not say what is
in there. This does. The interval holds **six protein-coding genes** —
*flp-15* (4.4 kb), *aqp-11* (1.5 kb), *dyf-2* (11.6 kb), *sid-2* (3.0
kb), *cul-2* (5.4 kb) and *cyn-2* (1.1 kb) — plus six non-coding genes
of 52–260 bp. *sid-2* lies **inside** *dyf-2* on the opposite strand,
which is why it takes its own row.

<div class="derived">

The two parents of the cross differ at **27 sites** across the 37 kb.
**Two of them alter a protein, and both are missense in *sid-2***:
`5V→5L` at 13,679,460, where JU1793 carries the alternate allele, and
**`96T→96K`** at 13,680,248, where JU2466 does. The rest are 3
synonymous, and 22 intronic, UTR or unannotated. **No difference in the
interval is HIGH-impact** — no stop gained or lost, no frameshift, no
splice-site or start-lost change in any of the twelve genes.

The census is close to complete rather than merely suggestive: of the
959 variant sites in the interval, exactly **one** has a no-call in
either parent, and **neither parent carries a divergent region**
overlapping it, so the absence of other coding differences is not an
alignment artefact. What this call set cannot see is structural
variation — it is short-read SNVs and small indels, so an inversion, a
large deletion or anything in a region that failed to align would not
appear.

</div>

<div class="aside">

<span class="ch">Why the rug is drawn</span>

**The denominator is the argument.** Two lollipops on an otherwise empty
axis would read as a search that found two hits. The rug shows the 25
other places the parents differ, so the panel reads as what it is: a
census of the interval in which two differences out of 27 are
protein-altering, and both fall in the same gene.

**Impact classes are bcftools csq terms mapped to the snpEff ladder**,
and the mapping is written out in `scripts/make_nil_interval_tables.R`
rather than assumed — HIGH is stop gained/lost, start lost, frameshift
and splice acceptor/donor; MODERATE is missense and inframe indels.

</div>

## The panel split at two loci

**Script** `scripts/diagnostic_genotype_splits.R`<br> **Figure**
`plots/diagnostics/genotype_splits.png`

Two loci, asked separately and then together: `IV:15,323,414`, the
strongest marker in the pooled scan, and `III:13,680,248`, *sid-2* T96K
— the allele the cross and the NILs implicate. Higher VST is more
resistant.

<img src="plots/diagnostics/genotype_splits.png" alt="Boxplots of the pooled pos-1 response split by genotype at the chromosome IV peak, at sid-2 T96K, and by both jointly." width="100%" />

<table>
<caption>
The marginal association each marker carries in the scan itself.
</caption>
<thead>
<tr>
<th style="text-align:left;">
Locus
</th>
<th style="text-align:right;">
AF
</th>
<th style="text-align:right;">
β
</th>
<th style="text-align:right;">
p (Wald)
</th>
<th style="text-align:right;">
−log₁₀p
</th>
<th style="text-align:right;">
Rank of 464,045
</th>
<th style="text-align:right;">
Single-locus R²
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
IV:15.32 Mb (GWAS peak)
</td>
<td style="text-align:right;">
0.208
</td>
<td style="text-align:right;">
+0.0258
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
8.84
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.159
</td>
</tr>
<tr>
<td style="text-align:left;">
sid-2 T96K
</td>
<td style="text-align:right;">
0.359
</td>
<td style="text-align:right;">
-0.0050
</td>
<td style="text-align:right;">
0.243
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
101,105
</td>
<td style="text-align:right;">
0.007
</td>
</tr>
</tbody>
</table>

<div class="derived">

Derived from plots/diagnostics/TABLE_genotype_splits.tsv and
TABLE_genotype_splits_marginal.tsv

</div>

<div class="aside">

<span class="ch">Chromosome IV is what this panel can see</span>

Splitting at the chromosome IV peak separates `48` strains carrying the
resistant allele (mean `+0.023`) from `182` carrying the other
(`−0.029`): Δ = `−0.0514`, Wilcoxon `p = 1.2e-07`, which is **12% of the
phenotypic range**. Its marginal statistic is the strongest in the scan
— rank `1` of 464,045, `−log₁₀p = 8.84`, single-locus R² `0.159` — and
in an additive two-locus fit it is the only term that carries anything:
β `−0.052`, `p = 8e-10`, adjusted R² `0.152`.

So there is more than *sid-2* segregating for this trait in the wild
population, and chromosome IV is the part of it the pooled GWAS is
actually powered to find.

</div>

<div class="caveat">

<span class="ch">Why *sid-2* is invisible here, and why that is not a
contradiction</span>

Splitting at T96K separates `83` strains carrying 96K (mean `−0.024`)
from `147` carrying 96T (`−0.015`). The direction is right — 96T more
resistant, matching the allele swaps — but Δ = `+0.0089` with Wilcoxon
`p = 0.99`, single-locus R² `0.007`, and a marginal rank of `101,105` of
464,045. In the additive fit, β = `−0.002`, `p = 0.74`.

**The two loci are not independent in the panel.** Fisher `p = 8.1e-05`,
odds ratio `0.20`: only `6` of 230 strains carry 96K on the chromosome
IV resistant background. The allele is largely confined to one
background, so a marginal test at T96K is asking a question this panel
cannot cleanly answer — which is a different statement from the allele
having no effect.

The evidence for *sid-2* is the cross, the NIL series and the allele
swaps, not the mapping. T96K moves JU1793 from `95%` to `53%` hatching,
roughly two-thirds of the `64`-point parental gap, and the direction
holds in all three backgrounds tested. What this diagnostic adds is the
honest account of the GWAS side: the scan does not support *sid-2*, it
is not powered to, and the reason is visible in the genotype table
rather than a matter of assertion.

One further caution on interval arithmetic. The r² ≥ 0.7 interval around
the chromosome III peak reaches *sid-2* only through a **single** marker
— the chromosome’s terminal marker at 13.7835 Mb, r² = `0.7014`, a
thousandth above the cutoff. Every other marker beyond 13.6 Mb,
*sid-2*’s own included, sits at r² = `0.5275`. At r² ≥ 0.71 the interval
collapses to 14 kb. So “the LD interval contains *sid-2*” is a
convex-hull effect at 0.7 and only a genuine statement at 0.5, and the
causal argument does not need it either way.

</div>

# Open before submission

Everything above is generated and verified. These are the items that
still need a decision or a number from outside the repository.

**Caption missing.** Figure 4 panel D has no entry in
`FIGURE_CAPTIONS.txt`. The wild-variation panel was added after that
file was generated, so the panel D caption above was written from the
script header and the variant table rather than copied from the caption
file. Fold it back so the two agree.

**Naming.** The caption file still describes Figures 1, 3 and 4 by their
superseded filenames and carries captions for the layout variants now in
`plots/legacy/`. The curated set is `Figure1_pos1`, `Figure2`,
`Figure3_quad`, `Figure4_sid2` — worth renumbering to 1–4 and S1–S12 in
one pass before submission.

**Verify.** Published *sid-2* allele positions (residues 32, 34, 168,
175, 199) come from a UniProt annotation for `G5EEV9` recorded in
earlier structure-modelling work. Check the residue numbers against
UniProt directly — they carry the proximity argument in Figure 4C.

**To fill.** `METHODS.txt` has 17 `[TO FILL]` markers — husbandry,
library prep, sequencing platform, editing protocol, the VST transform
definition, what the bootstrap resampled, and GEMMA / GATK / DeepTMHMM /
AlphaFold versions. The Dryad DOI is a placeholder in
`DATA_AVAILABILITY.md`.

**Repository.** Two figure scripts in `scripts/legacy/` still resolve
their working directory through the RStudio API and cannot run headless.
`.git` remains \~380 MB because the untracked bulk data is still in
history; shrinking it needs a history rewrite and a force-push.

# Figure manifest

<table>
<thead>
<tr>
<th style="text-align:left;">
Number
</th>
<th style="text-align:left;">
File
</th>
<th style="text-align:right;">
Size (KB)
</th>
<th style="text-align:right;">
Modified
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Figure S1
</td>
<td style="text-align:left;">
SUPP_FIG_XX_simulation_depth
</td>
<td style="text-align:right;">
448
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S2
</td>
<td style="text-align:left;">
SUPP_FIG_XX_dilution_validation
</td>
<td style="text-align:right;">
463
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure 1
</td>
<td style="text-align:left;">
Figure1_pos1
</td>
<td style="text-align:right;">
811
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S3
</td>
<td style="text-align:left;">
SUPP_FIG_XX_baugh_per_sample_frequencies
</td>
<td style="text-align:right;">
396
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S4
</td>
<td style="text-align:left;">
SUPP_FIG_XX_bootstrap_propagation_checks
</td>
<td style="text-align:right;">
448
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S5
</td>
<td style="text-align:left;">
SUPP_FIG_XX_downsample_per_sample
</td>
<td style="text-align:right;">
269
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S6
</td>
<td style="text-align:left;">
SUPP_FIG_XX_original_pos1_dfreq_rep_correlation
</td>
<td style="text-align:right;">
341
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S7
</td>
<td style="text-align:left;">
SUPP_FIG_plate_vs_paaby_vs_pos1original
</td>
<td style="text-align:right;">
199
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure 2
</td>
<td style="text-align:left;">
Figure2
</td>
<td style="text-align:right;">
1335
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S8
</td>
<td style="text-align:left;">
SUPP_FIG_XX_pooled_phenotype_ranks
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S9
</td>
<td style="text-align:left;">
SUPP_FIG_XX_cross_contrast_panels
</td>
<td style="text-align:right;">
1011
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure 3
</td>
<td style="text-align:left;">
Figure3_quad
</td>
<td style="text-align:right;">
157
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S10
</td>
<td style="text-align:left;">
SUPP_FIG_XX_nil_hatching_full
</td>
<td style="text-align:right;">
163
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure 4
</td>
<td style="text-align:left;">
Figure4_sid2
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S11
</td>
<td style="text-align:left;">
SUPP_FIG_XX_n2_swap_dose
</td>
<td style="text-align:right;">
240
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S12
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_allele_swaps_full
</td>
<td style="text-align:right;">
562
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S13
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_allele_in_panel
</td>
<td style="text-align:right;">
486
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S14
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_electrostatics
</td>
<td style="text-align:right;">
794
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S15
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_local_charge
</td>
<td style="text-align:right;">
96
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S16
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_model_confidence
</td>
<td style="text-align:right;">
441
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S17
</td>
<td style="text-align:left;">
Figure2_no_cross_qtl
</td>
<td style="text-align:right;">
1515
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S18
</td>
<td style="text-align:left;">
SUPP_FIG_XX_nil_interval_genes
</td>
<td style="text-align:right;">
112
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
<tr>
<td style="text-align:left;">
Figure S19
</td>
<td style="text-align:left;">
SUPP_FIG_XX_sid2_briggsae_alignment
</td>
<td style="text-align:right;">
202
</td>
<td style="text-align:right;">
2026-09-09 21:44
</td>
</tr>
</tbody>
</table>

<div class="tnote">

All twenty-one figures rebuild from `supplemental_data/` with `data/`
absent, and are pixel-identical across repeated runs. Captions
transcribed from `FIGURE_CAPTIONS.txt`; every number in the caption
prose was taken from the generating scripts’ console output, and every
table marked *derived* is recomputed from the deposit each time this
file knits.

</div>
