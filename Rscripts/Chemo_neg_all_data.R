library(tibble)
library(factoextra)
library(ggplot2)

#Importing the data
Chemo_neg <- read.csv("~/Documents/4th year at Mac/Spring 2018/Knights Lab/Chemo_neg_all_samples_compounds.csv")

#Import the csv file with labels that can be used to correctly label the data
Chemo_neg_labels <- read.csv("~/Documents/4th year at Mac/Spring 2018/Knights Lab/Chemo_neg_all_samples_labels.csv")

#Changing the first columns from data to row names
Chemo_neg <- column_to_rownames(df = Chemo_neg, var = "Compound")
Chemo_neg_labels <- column_to_rownames(df = Chemo_neg_labels, var = "file_name")

#Changing the data set from a frame to a matrix:
Chemo_neg <- as.matrix(Chemo_neg)

#Get the samples to overlap between Chemo_neg and Chemo_neg_labels
common.ids <- intersect(colnames(Chemo_neg), rownames(Chemo_neg_labels))
common.ids <- sort(common.ids)
Chemo_neg <- Chemo_neg[, common.ids]
Chemo_neg_labels <- Chemo_neg_labels[common.ids,]

#Running a PCA analysis using prcomp (source: https://cdn.intechopen.com/pdfs-wm/52527.pdf)
Chemo_neg <- Chemo_neg[rowSums(Chemo_neg) != 0, ]
prcomp(Chemo_neg, retx = T, center = T, scale = F) 

newzero <- min(Chemo_neg[Chemo_neg > 0]) / 2
Chemo_neg[Chemo_neg == 0] <- newzero

log_Chemo_neg <- log(Chemo_neg, 2)

paretoscale <- function(z) {
  rowmean <- apply(z, 1, mean) 
  rowsd <- apply(z, 1, sd) 
  rowsqrtsd <- sqrt(rowsd) 
  rv <- sweep(z, 1, rowmean, "-") #subtracts the row mean from each value
  rv <- sweep(rv, 1, rowsqrtsd, "/") #divides by the square root of the stdv for each value
  return(rv)
  }

pareto.log_Chemo_neg <- paretoscale(log_Chemo_neg)
pca <- prcomp(t(pareto.log_Chemo_neg), center=F, scale=F)
pcaresults <- summary(pca)

score.Chemo_neg <- as.data.frame(pcaresults$x)

data_full <- score.Chemo_neg[, c(1:3)] # subset columns 1-3

#Merge the data from 'data_full' and 'Chemo_neg_labels' to include all pertinent info. in a single table
#This also makes it easier to plot using ggplot2
data_merged <- merge(data_full, Chemo_neg_labels, by=0)

ggplot(data_merged, aes(PC1, PC2, group = patient_id)) + geom_point(aes(color=CHEMO_PRE_POST)) + geom_line(size=0.1) + scale_color_manual(values = c("PRE_CHEMO" = "coral", "POST_CHEMO" = "cyan3", "POOLED" = "gray"))

ggplot(data_merged, aes(PC1, PC2, group = patient_id)) + geom_point(aes(color=HCT_PRE_POST)) + geom_line(size=0.1) + scale_color_manual(values = c("PRE_HCT" = "coral", "POST_HCT" = "cyan3", "POOLED" = "gray" ))

ggplot(data_merged, aes(PC1, PC2, group = patient_id)) + geom_point(aes(color=Patient_status)) + geom_line(size=0.1) + scale_color_manual(values = c("Pre_chemo" = "darkorchid1", "Post_chemo_pre_HCT" = "deepskyblue", "Post_HCT" = "olivedrab3", "POOLED" = "gray"))

#Look up how to make dots bigger (so they're more visible)
#I tried to use geom_segment initially for the line segments ~ didn't work 
#Remove NA's from patient ID's --> not necessary; R just excluded them from the analysis
