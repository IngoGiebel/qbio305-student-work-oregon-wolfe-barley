library(ggplot2)    # ggplot2 library for plotting                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
library(car)        # car library to make a QQplot
options(scipen=999) # remove scientific notation of very low integers

if (!require("ggplot2",character.only = TRUE)){
  install.packages("ggplot2")
} 
library(ggplot2)    

if (!require("qqman",character.only = TRUE)){
  install.packages("qqman")
}
library(qqman)


# load the data
dta.raw <- read.csv("dta---ANOVA_QTLanalysis.csv")

# view data
head(dta.raw, n = 15)

#H0: means of all groups are equal 
#HA: at least one group mean differs significantly

# conduct the raw ANOVA
model.raw <- aov(Height ~ Genotype, 
                 data = dta.raw)

# extract the residuals
resid.raw <- model.raw$residuals
resid.raw

# plot the histogram
ggplot(data = data.frame(1:length(resid.raw), resid.raw), 
       aes(x = resid.raw)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.raw), 
                            sd = sd(resid.raw)),
                color="darkred")

# plot qq-plot
qqPlot(resid.raw, envelope = F)

# check for homoscedasticity and outliers 
plot(model.raw, which = 1)

# remove outliers
dta.clean <- dta.raw[-c(14,40,182),]

# conduct the ANOVA with the cleaned data
model.clean <- aov(Height ~ Genotype, data= dta.clean) 

# extract the residuals of the new ANOVA
resid.clean <- model.clean$residuals

# check normality and homoscedasticity again
shapiro.test(resid.clean)

ggplot(data = data.frame(1:length(resid.clean), resid.clean), 
       aes(x = resid.clean)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.clean), 
                            sd = sd(resid.clean)),
                color="darkred")
plot(model.clean, which = 1)
plot(model.clean, which = 2)

# create the final model
model.aov <- aov(Height ~ Genotype, 
                 data = dta.clean)
# check the ANOVA table
summary(model.aov)


# p-value is less than a chosen significance level (often 0.05),
# we can reject the null hypothesis and conclude that there is 
# a statistically significant difference in plant heights among the groups


dta <- read.csv("dta---ANOVA_QTLanalysis.csv")
map <- read.csv("input_map.csv")

# view data
dim(dta)
dta[1:10,1:10]

# extract the columns group, plant_height and the column containing the 6th marker ??
dta.singleM <- dta[,c(2,3,4,10)]
dta.singleM

# rename the plant_height and the marker column
colnames(dta.singleM) <- c("group","plant_height", "spike_length", "M_bPb6318")
dta.singleM

# plot a boxplot of allele vs plant height
ggplot(dta.singleM, aes(x=M_bPb6318, y=plant_height, fill=M_bPb6318)) + 
  geom_boxplot()

# line in the middle = median -> skewed distribution, not normal
# mean = middle of the two boxes
# QTA analysis now it is done with linear regression, not with lots of anovas

############################## question ##############################
# Create the same boxplot, but for groups instead of alleles.
# Do you expect to find a significant effect of this marker? Why?
# Answer:

ggplot(dta.singleM, aes(x=group, y=plant_height, fill=group)) + 
  geom_boxplot()

######################################################################

# conduct ANOVA
res.aov.singleM <- aov(plant_height ~ M_bPb6318, dta.singleM)
summary(res.aov.singleM)

############################## question ##############################
# What are the p-values? Do alleles or groups have a significant effect and
# What does that mean?
# Answer:0.0989 p value is higher that 0.05. Therefore plant_height and marker don't have 
#a significant effect on groups

######################################################################

# loop through all markers and conduct ANOVAs for each one of them
lis.out <- list()
for(m in 5:ncol(dta)){
  # extract the first two columns (group + spike length) and the m-th column (=marker)
  loop.dta <- dta[,c(2,3,m)]
  # save the name of the current marker
  loop.marker.name <- colnames(loop.dta)[3]
  # in the model in the aov command, we will have to define the independent 
  # variable name. since the marker name changes in each loop, we need to 
  # to change the column names here to have the same marker name in each loop
  colnames(loop.dta) <- c("group","trait", "allele")
  # conduct one-way ANOVA
  loop.aov <- aov(trait ~ allele + group, loop.dta)
  # extract the allele's p-value
  loop.pval.allele <- summary(loop.aov)[[1]][1,5]
  # extract the marker's genetic position from the genetic map
  loop.map <- map[which(map$marker == loop.marker.name),]
  # create the output data frame
  loop.df.out <- data.frame(loop.map,
                            pval.allele=loop.pval.allele)
  # save the output data frame in the list
  lis.out[[m]] <- loop.df.out
}

# combine the loop's output data frames into a single data frame
res.aov <- do.call("rbind", lis.out)

# have a look at the output data frame
head(res.aov)
dim(res.aov)

# calculate the negative logarithm of the Bonferroni corrected significance threshold
sig.threshold.BonfCorrected <- -log(0.05/nrow(res.aov))

# create a data frame that follows the requirements of the manhattan command of the qqman library
dta.plot <- data.frame(CHR = as.numeric(gsub("H","",res.aov$chr)),
                       BP = res.aov$pos,
                       SNP = 1:nrow(res.aov),
                       pval.allele = res.aov$pval.allele)

# plot genome-wide p-values for markers
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.allele", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Marker effect")

# Do you have any regions in the genome with significant marker effects?
# chromosome 2 

                           
# GENETIC LINKAGE MAP

#load the qtl library for linkage analysis
library(qtl)

# read filtered marker data of chromosomes 1H and 2H and create linkage map
chr1_2 <- read.cross("csv", ".", "barley_marker_data_chr1_2.csv", genotypes=c("a","b"), 
                     alleles=c("a", "b"), crosstype="dh", estimate.map = FALSE) #crosstype:double haplotype

#get overview of data
summary(chr1_2) 

#horizontal "line" is individual with a lot of missing genotyping data
plotMissing(chr1_2) 

#estimate pairwise recombination fractions
chr1_2 <- est.rf(chr1_2) 
head(chr1_2[["rf"]][,1:6])

#plot pairwise recombination frequency of markers
plotRF(chr1_2, alternate.chrid=TRUE) 

#partition markers into linkage groups
chr1_2 <- formLinkageGroups(chr1_2, reorgMarkers=TRUE) 

#check linkage groups
plotRF(chr1_2, alternate.chrid=TRUE) 

#plot the markers on a linkage map
plotMap(chr1_2, main="", show.marker.names=T) 


# try with chromosome 3 and 4 ##################################################

chr3_4 <- read.cross("csv", ".", "barley_marker_data_chr3_4.csv", genotypes=c("a","b"), 
                     alleles=c("a", "b"), crosstype="dh", estimate.map = FALSE)

# Get overview of the data
summary(chr3_4)
plotMissing(chr3_4)

# Calculate and plot recombination fractions
chr3_4 <- est.rf(chr3_4) 
head(chr3_4[["rf"]][,1:6])

plotRF(chr3_4, alternate.chrid=TRUE)

# Form linkage groups and plot linkage map
chr3_4 <- formLinkageGroups(chr3_4, reorgMarkers=TRUE) #partition markers into linkage groups

plotRF(chr3_4, alternate.chrid=TRUE) #check linkage groups, do you see 2 chromosomes?

plotMap(chr3_4, main="", show.marker.names=T) #plot the markers on a linkage map

################################################################################

# Mapping of leaf variegation trait = linkage map

#read in linkage map and qualitative phenotypes as denoted in marker annotation
owb <- read.cross("csv", ".", "owb_linkage_map_qualt_phenotypes_leafvar.csv", genotypes=c("a","b"), 
                  alleles=c("a", "b"), crosstype="dh")

#have a look at the complete linkage map
plotMap(owb, main="", show.marker.names=T)

#Example: mapping the locus responsible for leaf variegation 
#To map the locus, calculate all pairwise recombination frequencies and LOD scores
leaf_var <- tryallpositions(owb, "leaf_var", error.prob=0)

#show the best linkage to a marker from each chromosome
summary(leaf_var)
#View(leaf_var)

#move the locus to the best marker position
owb <- movemarker(owb, "leaf_var", "2H", 192.7990)

# chromosome 2H , position 192.7990 -> linked to highest lod score = lod score:3.8330018846980835612
# scssr08447-qter (10-*) closest marker name 

#update map
plotMap(owb, main="", show.marker.names=T)


