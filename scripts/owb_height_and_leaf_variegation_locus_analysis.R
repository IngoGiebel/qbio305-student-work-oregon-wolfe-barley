# Load libraries ----------------------------------------------------------

# For determining the project root
library(here)
# For transforming and representing data
library(tidyverse)
# For plotting
library(ggpubr)
# For Levene's test
library(car)
# For linkage analysis
library(qtl)


# File paths --------------------------------------------------------------

dta_owb_genetic_map_file_path <- here::here(
  "data",
  "dta_owb_genetic_map.csv"
)

dta_owb_height_and_spike_length_by_genotype_file_path <- here::here(
  "data",
  "dta_owb_height_and_spike_length_by_genotype.csv"
)

dta_owb_linkage_map_qualt_phenotypes_file_path <- here::here(
  "data",
  "dta_owb_linkage_map_qualt_phenotypes.csv"
)


# Set R options -----------------------------------------------------------

options(scipen = 1)


# Load the raw data -------------------------------------------------------

# Read in the raw quantitative data of the different OWB genotypes
dta_quant_raw <- car::strings2factors(readr::read_csv(
  dta_owb_height_and_spike_length_by_genotype_file_path
))

# Read in the genetic map of OWB
dta_owb_genetic_map <- readr::read_csv(
  dta_owb_genetic_map_file_path
)

# Read in linkage map and qualitative phenotypes as denoted in marker annotation
dta_linkage_map_qualt_phenotypes <- qtl::read.cross(
  format = "csv",
  file = dta_owb_linkage_map_qualt_phenotypes_file_path,
  genotypes = c("a", "b"),
  alleles = c("a", "b"),
  crosstype = "dh"
)


# ANOVA height vs genotype (raw data) -------------------------------------

# Densitiy plot comparing the height of the different genotypes
ggpubr::ggdensity(
  data = dta_quant_raw,
  x = "Height",
  add = "mean",
  color = "Genotype",
  fill = "Genotype",
  xlab = "Height in cm",
  ylab = "Density",
) +
  ggplot2::theme(legend.justification = c("left", "top")) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      title.position = "top",
      nrow = 6,
      byrow = TRUE
    )
  )

# Summary statistics per genotype
dta_quant_raw |>
  dplyr::select(Genotype, Height) |>
  dplyr::group_by(Genotype) |>
  dplyr::summarise(
    Count = n(),
    Mean = mean(Height),
    Sd = sd(Height)
  )

# Fit an analysis of variance model for the height vs genotype
height_raw_aov <- aov(Height ~ Genotype, data = dta_quant_raw)

# Summary of the ANOVA analysis
summary(height_raw_aov)

# Conclusion: The genotype has a significant effect on the plant height
# (p-value < 2e-16).


# ANOVA - Check the homogeneity of variance assumption --------------------

# The ANOVA test assumes that, the data are normally distributed and the
# variance across groups are homogeneous. This is assessed with some
# diagnostic plots and tests.

# Plot the residuals vs fit in order to check the homogeneity of variance
# assumption
plot(height_raw_aov, which = 1)

# Levene's test for homogeneity of variance across groups
car::leveneTest(Height ~ Genotype, data = dta_quant_raw)

# Conclusion: the variances are roughly homogeneous. The p-value for the
# Levene's test of the raw data is 0.7938. This means that the variances
# within the groups are not significantly different from each other, i.e.,
# the homogeneity assumption of the variance is met.

# Outliers: the rows 14, 40, and 182 are detected as outliers, which can
# severely affect normality and homogeneity of the variance. Therefore, these
# outliers (in the genotypes OWB33, OWB76, and OWB54) will be removed to meet
# the test assumptions.


# ANOVA - Check the normality assumption ----------------------------------

# Plot a histogram of the residuals
ggplot2::ggplot(
  data = data.frame(
    1:length(height_raw_aov$residuals),
    height_raw_aov$residuals
  ),
  mapping = ggplot2::aes(x = height_raw_aov$residuals),
) +
  ggplot2::geom_histogram(
    mapping = ggplot2::aes(y = ..density..),
    color = "darkgrey",
    bins = 50
  ) +
  ggplot2::stat_function(
    fun = dnorm,
    args = list(
      mean = mean(height_raw_aov$residuals),
      sd = sd(height_raw_aov$residuals)
    ),
    color = "darkred"
  ) +
  ggplot2::labs(
    x = "Residuals of the ANOVA height ~ genotype",
    y = "Density"
  )

# Q-Q plot of the ANOVA residuals
plot(height_raw_aov, which = 2)

# Cook's distance plot
plot(height_raw_aov, which = 4)

# Shapiro-Wilk test for normality
shapiro.test(height_raw_aov$residuals)

# Conclusion: According to the Shapiro-Wilk test for normality, it is unlikely
# that the raw data is normally distributed (p-values < 2.2e-16). Nevertheless,
# the plots reveal that the data conform at least roughly to a normal
# distribution, in particular without the three outliers.
# Therefore, the further analysis will be done using ANOVA with cleansed data
# (without the outliers).


# ANOVA height vs genotype (cleansed data w/out outliers) -----------------

# Remove outliers
dta_quant_clean <- dta_quant_raw[-c(14, 40, 182), ]

# Fit an analysis of variance model for the height (cleansed data)
height_clean_aov <- aov(Height ~ Genotype, data = dta_quant_clean)

# Summary of the ANOVA analysis
summary(height_clean_aov)

# Plot a histogram of the residuals
ggplot2::ggplot(
  data = data.frame(
    1:length(height_clean_aov$residuals),
    height_clean_aov$residuals
  ),
  mapping = ggplot2::aes(x = height_clean_aov$residuals),
) +
  ggplot2::geom_histogram(
    mapping = ggplot2::aes(y = ..density..),
    color = "darkgrey",
    bins = 50
  ) +
  ggplot2::stat_function(
    fun = dnorm,
    args = list(
      mean = mean(height_clean_aov$residuals),
      sd = sd(height_clean_aov$residuals)
    ),
    color = "darkred"
  ) +
  ggplot2::labs(
    x = "Residuals of the ANOVA height ~ genotype",
    y = "Density"
  )

# Q-Q plot of the ANOVA residuals
plot(height_clean_aov, which = 2)

# Cook's distance plot
plot(height_clean_aov, which = 4)

# Shapiro-Wilk test for normality
shapiro.test(height_clean_aov$residuals)

# p-value = 7.823e-14
# Still very unlikely, but much more likely than for the raw data.

# Conclusion: According to the Shapiro-Wilk test for normality, it is still
# unlikely that the cleansed data is normally distributed (p-value = 7.823e-14).
# Nevertheless, the plots reveal that the data conform at least roughly to a
# normal distribution, in particular without the removed three outliers.


# ANOVA of the plant height for all markers -------------------------------

# Tibble for the markers, theirs positions and significance values for the
# plant height
marker_pval_df <- tibble::tibble(
  marker = character(),
  chr = character(),
  pos = numeric(),
  pval = numeric())

# Loop through all markers and conduct ANOVAs for each one of them
for (m in 5:ncol(dta_quant_clean)) {
  marker_name <- colnames(dta_quant_clean)[m]
  # Extract the marker's genetic position from the genetic map
  marker_chr_pos <- dplyr::filter(dta_owb_genetic_map, marker == marker_name)
  # Conduct an one-way ANOVA of the plant height vs the marker
  marker_aov <- aov(
    as.formula(paste("Height", "~", marker_name)),
    data = dta_quant_clean)
  # Extract the marker's p-value from the ANOVA table
  pval <- summary(marker_aov)[[1]][1, 5]
  # Add a new row for the marker, its position, and its p-value
  marker_pval_df <- marker_pval_df |> tibble::add_row(marker_chr_pos, pval)
}

# Create a data frame that follows the requirements of the
# Manhattan command of the qqman library
marker_pval_manhattan_plot_df <- marker_pval_df |>
  dplyr::mutate(CHR = as.numeric(gsub("H", "", chr))) |>
  dplyr::rename(BP = pos) |>
  dplyr::mutate(SNP = 1:n(), .before = pval) |>
  dplyr::select(CHR, BP, SNP, pval)

# Calculate the Bonferroni corrected significance threshold. The Bonferroni
# correction adjusts the chosen significance level (here 5 %) by dividing it
# by the number of tests.
sig_thld_bonfcorr <- 0.05 / nrow(marker_pval_df)

# Manhattan plot of the genome-wide p-values for the markers
qqman::manhattan(
  marker_pval_manhattan_plot_df,
  genomewideline = -log(sig_thld_bonfcorr),
  suggestiveline = FALSE,
  logp = TRUE,
  p = "pval",
  lwd = 3,
  ylab = "p-value (neg. log)",
  main = "Marker effect on the plant height"
)

# Determine all gene positions which have p-value less than the Bonferroni
# corrected p-value
marker_pval_sig_df <- marker_pval_df |>
  dplyr::filter(pval < sig_thld_bonfcorr) |>
  dplyr::rename(Chromosome = chr) |>
  dplyr::group_by(Chromosome) |>
  dplyr::summarise(
    Count = n(),
    "Min. p-value" = min(pval),
    "Median p-value" = median(pval),
    "Min. chr position" = min(pos),
    "Max. chr position" = max(pos)
  )

# Result:

# Chr. | Count | Min. p-val | Median p-val | Min. chr pos | Max. chr pos

# 1H      37     1.74e- 7     4.24e- 6       52.7           157.
# 2H     180     1.10e-42     1.82e-14       1.26           187.
# 3H      68     5.43e- 8     2.57e- 6       58.6           146.
# 5H      38     1.48e- 8     4.45e- 7       137.           238.
# 6H       1     6.37e- 6     6.37e- 6       129.           129.

# Conclusion: The plant height is mainly significantly effected by
# a genomic region within gene 2. But there are also markers within
# the genes 1H, 3H, 5H, and 6H which have a significant impact on
# the plant height.


# Map the qualitative traits ----------------------------------------------

# Plot the markers on a linkage map
qtl::plotMap(
  dta_linkage_map_qualt_phenotypes,
  main = "",
  show.marker.names = TRUE
)

# Mapping the locus responsible for the leaf variegation
# To map the locus, calculate all pairwise recombination frequencies
# and LOD scores
leaf_variegation <- qtl::tryallpositions(
  dta_linkage_map_qualt_phenotypes,
  marker = "leaf_variegation",
  error.prob = 0
)
# Show the best linkage to a marker from each chromosome
summary(leaf_variegation)

# Move the locus to the best marker position
dta_linkage_map_qualt_phenotypes <- qtl::movemarker(
  dta_linkage_map_qualt_phenotypes,
  marker = "leaf_variegation",
  newchr = "2H",
  newpos = 192.8
)

# Update the linkage map plot
qtl::plotMap(
  dta_linkage_map_qualt_phenotypes,
  main = "",
  show.marker.names = TRUE
)

# Conclusion: The leaf variegation of the Oregon Wolfe barley collection
# is significantly effected by the chromosome 2H at position 192.8.

# Remark: According to
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8066651/
#
# the leaf variegation is controlled by the gene Wst. Therefore, our findings
# are in accordance with the findings of that paper.
# See figure
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8066651/figure/plants-10-00694-f003/
# for an illustration of the Wst locus in gene 2H.
