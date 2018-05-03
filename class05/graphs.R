#' ---
#' title: "Bioinformatics Class 5"
#' author: "Jason Patrick Bennett"
#' output:
#'   html_document:
#'     code_folding: hide
#' ---

# Class 5 Graphs

# --------------------- Section 1 ------------------------

# boxplot
boxplot( rnorm(1000,0), horizontal = TRUE )

# histogram
hist( rnorm(1000,0) )

# summary shows basic statistics about the data being plotted
summary( rnorm(1000,0) )

# Read first data file
weight_df <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)

# Plot the dataframe
plot(weight_df, type="b", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Infant Age v. Weight")

# 1B -------------------------

# Read second data file
feature_df <- read.table("bimm143_05_rstats/feature_counts.txt", sep="\t", header=TRUE)

#par( mar=c(5,11,4,2) )

barplot(feature_df[,2], horiz=TRUE, names.arg=feature_df[,1], main="Total Number of Features", las=1)






# ----------------------------- Section 2 ------------------------------

# 2A ------------------------------

# Read file
mfFile <- "bimm143_05_rstats/male_female_counts.txt"
mf_count <- read.table(mfFile, header=TRUE, sep="\t")
color_vector=rainbow(nrow(mf_count))

# Barplot of file: rainbow colors
barplot(mf_count$Count, col=color_vector, names.arg=mf_count$Sample)

# Barplot of file: red and blue
barplot(mf_count$Count, col=c("red2", "blue2"), names.arg=mf_count$Sample)

# 2B ------------------------------

# Read file
upDownFile <- "bimm143_05_rstats/up_down_expression.txt"
up_down <- read.table(upDownFile, header=TRUE, sep="\t")
palette(c("Blue", "Grey", "Red"))

# Scatterplot of Condition1 v. Condition2
plot(up_down$Condition1, up_down$Condition2, col=up_down$State)

# 2C ------------------------------

# Run r file and read text file
source("bimm143_05_rstats/color_to_value_map.r")
exprMethFile <- "bimm143_05_rstats/expression_methylation.txt"

# Read Expression Methylation File
exp_meth <- read.delim(exprMethFile)
high <- max(exp_meth$expression)
low <- min(exp_meth$expression)
highLow <- c(high, low)

# Scatterplot of Promoter v. Gene
setOfColors <- (colorRampPalette(c("grey", "red")) (100))
methColorVector <- map.colors(exp_meth$expression, highLow, setOfColors)
plot(exp_meth$promoter.meth, exp_meth$gene.meth, col=methColorVector)