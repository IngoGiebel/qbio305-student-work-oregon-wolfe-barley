# Load libraries ----------------------------------------------------------

# For determining the project root
library(here)
# For plotting
library(ggplot2)
# For qq-plotting
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

# Remove scientific notation of very low integers
options(scipen = 999)


# QTL analysis - load the raw data ----------------------------------------

dta_quant_raw <- read.csv(dta_owb_height_and_spike_length_by_genotype_file_path)

dta_owb_genetic_map <- read.csv(dta_owb_genetic_map_file_path)

# Read in linkage map and qualitative phenotypes as denoted in marker annotation
dta_linkage_map_qualt_phenotypes <- qtl::read.cross(
  format = "csv",
  file = dta_owb_linkage_map_qualt_phenotypes_file_path,
  genotypes = c("a", "b"),
  alleles = c("a", "b"),
  crosstype = "dh"
)

# ANOVA height vs genotype (raw data) -------------------------------------

# Fit an analysis of variance model for the height
height_raw_aov <- aov(Height ~ Genotype, data = dta_quant_raw)

# Plot a histogram of the raw residuals
ggplot2::ggplot(
  data = data.frame(
    1:length(height_raw_aov$residuals),
    height_raw_aov$residuals
  ),
  mapping = ggplot2::aes(x = height_raw_aov$residuals)
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
  )

# Plot qq-plot
car::qqPlot(height_raw_aov$residuals, envelope = FALSE)

# Conclusion: The residuals are not normally distributed. The first visual
# inspection of the plant heights (disregarding the genotype) does not uncover
# any outliers.

# Check whether or not the dependent height variable is normally distributed
hist(dta_quant_raw$Height)

# Check for homoscedasticity and outliers
plot(height_raw_aov$residuals, which = 1)

# Conclusion: There are two greater outliers (in the genotypes OWB33 and OWB76).

# Create a linear model of the raw data and print a summary
summary(lm(Height ~ Genotype, data = dta_quant_raw))


# ANOVA height vs genotype (clean data w/out outliers) --------------------

# Remove outliers
dta_quant_clean <- dta_quant_raw[-c(14, 40), ]

# Fit an analysis of variance model for the height (cleaned data)
height_clean_aov <- aov(Height ~ Genotype, data = dta_quant_clean)

# Plot a histogram of the raw residuals
ggplot2::ggplot(
  data = data.frame(
    1:length(height_clean_aov$residuals),
    height_clean_aov$residuals
  ),
  mapping = ggplot2::aes(x = height_clean_aov$residuals)
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
  )

# Plot qq-plot
car::qqPlot(height_clean_aov$residuals, envelope = FALSE)

# Check normality and homoscedasticity again
shapiro.test(height_clean_aov$residuals)
plot(height_clean_aov$residuals, which = 1)

# Check the ANOVA table
summary(height_clean_aov)

# Conclusion: The two outliers are removed. The cleaned data is roughly
# normally distributed (according to the visual inspection). The genotype
# has a significant effect on the plant height.


# ANOVA of the plant height for all markers -------------------------------

# Loop through all markers and conduct ANOVAs for each one of them
lis_out <- list()
for (m in 5:ncol(dta_quant_clean)) {
  # Extract the plant height and the m-th column (=marker)
  loop_dta <- dta_quant_clean[, c(3, m)]
  # Save the name of the current marker
  loop_marker_name <- colnames(loop_dta)[2]
  # In the model in the aov command, we will have to define the independent
  # variable name. Since the marker name changes in each loop, we need to
  # to change the column names here to have the same marker name in each loop.
  colnames(loop_dta) <- c("Height", "Allele")
  # Conduct one-way ANOVA
  loop_aov <- aov(Height ~ Allele, loop_dta)
  # Extract the allele's p-value
  loop_pval_allele <- summary(loop_aov)[[1]][1, 5]
  # Extract the marker's genetic position from the genetic map
  loop_map <- dta_owb_genetic_map[
    which(dta_owb_genetic_map$marker == loop_marker_name),
  ]
  # Create the output data frame
  loop_df_out <- data.frame(loop_map, pval_allele = loop_pval_allele)
  # Save the output data frame in the list
  lis_out[[m]] <- loop_df_out
}
# Combine the loop's output data frames into a single data frame
res_aov <- do.call("rbind", lis_out)


# Manhattan plot ----------------------------------------------------------

# Create a data frame that follows the requirements of the
# Manhattan command of the qqman library (run as is)
df_manhattan_plot <- data.frame(
  CHR = as.numeric(gsub("H", "", res_aov$chr)),
  BP = res_aov$pos,
  SNP = 1:nrow(res_aov),
  pval_marker = res_aov$pval_allele
)

# Calculate the negative logarithm of the Bonferroni corrected significance
# threshold
sig_threshold_BonfCorrected <- -log(0.05 / nrow(res_aov))

# Plot genome-wide p-values for markers
qqman::manhattan(
  df_manhattan_plot,
  genomewideline = sig_threshold_BonfCorrected,
  suggestiveline = FALSE,
  logp = TRUE,
  p = "pval_marker",
  lwd = 3,
  ylab = "-log10(p-values)",
  main = "Marker effect on the plant height"
)

# Conclusion: The plant height is significantly effected by a genomic region
# within gene 2.


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
