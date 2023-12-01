############################### load libraries and data (run as is)
if (!require("ggplot2",character.only = TRUE)){
  install.packages("ggplot2")
} 
library(ggplot2)    

if (!require("qqman",character.only = TRUE)){
  install.packages("qqman")
}
library(qqman)

options(scipen=999) # remove scientific notation of very low integers

dta <- read.csv("student_data_with_gt.csv")
map <- read.csv("input_map.csv")
###############################

# always have a look at the data first
dim(dta)
dta[1:10,1:10]

# extract the columns group, plant_height and the column containing the 6th marker
dta.singleM <- dta[???]
# to have more generic code, rename the plant_height and the marker column
colnames(dta.singleM) <- c("group", "trait", "allele")

############################## question ##############################
# What is the name of the marker that you have extracted?
# Does plant height has the correct data structure? Use the str command to 
# check and convert the data if necessary (use as.numeric).
# Answer:

######################################################################

# plot a boxplot of allele vs plant height
ggplot(dta.singleM, aes(x=???, y=plant_height, fill=allele)) + 
  geom_boxplot()

############################## question ##############################
# Create the same boxplot, but for groups instead of alleles.
# Do you expect to find a significant effect of this marker? Why?
# Answer:

######################################################################

# conduct ANOVA
res.aov.singleM <- aov(???)
summary(res.aov.singleM)

############################## question ##############################
# What are the p-values? Do alleles or groups have a significant effect and
# What does that mean?
# Answer:

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
head(???)
dim(???)

# calculate the negative logarithm of the Bonferroni corrected significance threshold
sig.threshold.BonfCorrected <- -log(???)

# create a data frame that follows the requirements of the 
# manhattan command of the qqman library (run as is)
dta.plot <- data.frame(CHR = as.numeric(gsub("H","",res.aov$chr)),
                       BP = res.aov$pos,
                       SNP = 1:nrow(res.aov),
                       pval.marker = res.aov$pval.allele,
                       pval.group = res.aov$pval.group)

# plot genome-wide p-values
# for markers
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.marker", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Marker effect")
# for group
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.group", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Group effect")


############################## question ##############################
# Do you have any regions in the genome with significant marker effects?
# If yes where?
# Answer:

######################################################################


############################## task ##################################
# repeat the analysis with spike length 
# (tip: you need to change only two numbers in the entire script)
######################################################################



############################## question ##############################
# Compare your marker effect plots of plant height and spike length
# Why are they so similar? (you cannot know that, just speculate)
# Answer:

######################################################################

