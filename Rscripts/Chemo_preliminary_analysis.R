#Importing the data
Chemo_neg_pre_post_data_for_r_studio <- read.csv("~/Documents/4th year at Mac/Spring 2018/Knights Lab/Chemo_neg_pre_post_data_for_r_studio.csv")
Chemo_neg_pre_post_2 <- Chemo_neg_pre_post_data_for_r_studio
library(tibble)
library(factoextra)
library(ggplot2)

#Changing the first column from data to row names
Chemo_neg_pre_post_2 <- column_to_rownames(df = Chemo_neg_pre_post_2, var = "Compound")

#Checking the data type/class
class(Chemo_neg_pre_post_2)

#Changing the data set from a frame to a matrix:
Chemo_neg_pre_post_matrix <- Chemo_neg_pre_post_2
Chemo_neg_pre_post_matrix <- as.matrix(Chemo_neg_pre_post_matrix)

#Running the Wilcoxon signed rank test to compare same patient, pre- & post-chemo samples:
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_53, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_53, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_54, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_54, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_75, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_75, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_78, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_78, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_81, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_81, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)
wilcox.test(Chemo_neg_pre_post_data_for_r_studio$Pre_chemo_112, Chemo_neg_pre_post_data_for_r_studio$Post_chemo_112, mu = 0, alt = "two.sided", paired = T, conf.int = T, conf.level = 0.95, exact = F)

#--------------------------------------------------------------------------------------------------------------
#Running a PCA analysis using prcomp: Second attempt (source: http://www.sthda.com/english/wiki/print.php?id=207#prcomp-and-princomp-functions)

#head(Chemo_neg_pre_post_matrix[, 1:6])
#res.pca <- prcomp(Chemo_neg_pre_post_matrix, scale = T)
#names(res.pca)
#head(res.pca$sdev)

#Eigenvalues
#eig <- (res.pca$sdev)^2
#variance <- eig*100/sum(eig)
#cumvar <- cumsum(variance)
#eig.Chemo_neg_pre_post_matrix <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)
#head(eig.Chemo_neg_pre_post_matrix)
#summary(res.pca)

#Visualizing principal component importance w/ a scree plot
#barplot(eig.Chemo_neg_pre_post_matrix[, 2], names.arg=1:nrow(eig.Chemo_neg_pre_post_matrix), main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col = "steelblue")

#Add connected line segments to the plot
#lines(x = 1:nrow(eig.Chemo_neg_pre_post_matrix), eig.Chemo_neg_pre_post_matrix[, 2], type="b", pch=19, col="red")

#Screeplot with eigenvalues
#fviz_screeplot(res.pca, ncp=10, choice="eigenvalue")

#--------------------------------------------------------------------------------------------------------------
#Running a PCA analysis using prcomp (source: https://cdn.intechopen.com/pdfs-wm/52527.pdf)

Chemo <- Chemo_neg_pre_post_matrix #Shortening the name
Chemo <- Chemo[rowSums(Chemo) != 0, ]
prcomp(Chemo, retx = T, center = T, scale = F) #initial prcomp 

# replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2) #Defining a function to replace zeros in the matrix
# Chemo <- apply(Chemo, 1, replacezero) #Replacing zeros in the data

newzero <- min(Chemo[Chemo > 0]) / 2
Chemo[Chemo == 0] <- newzero

logdata <- log(Chemo, 2) #Creating a log version of the data
paretoscale <- function(z) {
  rowmean <- apply(z, 1, mean) 
  rowsd <- apply(z, 1, sd) 
  rowsqrtsd <- sqrt(rowsd) 
  rv <- sweep(z, 1, rowmean, "-") #subtracts the row mean from each value
  rv <- sweep(rv, 1, rowsqrtsd, "/") #dividing by the square root of the stdv for each value
  return(rv)
  }

pareto.logdata <- paretoscale(logdata)
pca <- prcomp(t(pareto.logdata), center=F, scale=F)
pcaresults <- summary(pca)
head(pcaresults)

# scree.data <- as.data.frame(pcaresults$importance)
score.data <- as.data.frame(pcaresults$x)
# loadings.data <- as.data.frame(pcaresults$rotation)
# write.csv(scree.data, "pca_scree.csv")
# write.csv(score.data, "pca_scores.csv")
# write.csv(loadings.data, "pca_loadings.csv")

# data <- read.csv("pca_scores.csv", header = T)
data <- score.data[, c(1:3)] # subset columns 1-3

data$Group <- c("Pre-chemo", "Post-chemo", "Pre-chemo", "Post-chemo", "Pre-chemo", "Post-chemo", "Pre-chemo", "Post-chemo", "Pre-chemo", "Post-chemo", "Pre-chemo", "Post-chemo")
data$Patient <- c("Patient_53", "Patient_53", "Patient_54", "Patient_54", "Patient_75", "Patient_75", "Patient_78", "Patient_78", "Patient_81", "Patient_81", "Patient_112", "Patient_112")

ggplot(data, aes(PC1, PC2)) + geom_point(aes(shape=Patient, color=Group))

# stat_ellipse(aes(fill=Group)
# geom_text(aes(label=rownames(data))
# c("Pre-chemo", "Post-chemo")
