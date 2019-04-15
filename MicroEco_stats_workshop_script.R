#Microbial ecology statistics workshop
#06/08/17
#Presented by Damien Finn; damien.finn@uqconnect.edu.au

#Libraries:
library("reshape2")
library("plyr")
library("vegan")
library("ggplot2")
library("RColorBrewer")
library("gplots")

df <- read.table(".../otu_table.txt", sep = "\t", header = T)
#Optional: change the OTU names to something more human friendly
colnames(df) <- as.character(gsub(".*f__", "", colnames(df)))
#Note dimensions of data; check if it's imported correctly
View(df)
nrow(df)
ncol(df)

#****Basic diversity metrics: H', Beta, Richness / Evenness****
Spe.df <- df[, 2:4791]
rows <- df[,1]
rownames(Spe.df) <- rows
#Alpha diversity (Shannon) of individual Wetland sites
H <- diversity(Spe.df)
#And evenness
J <- H/log(specnumber(Spe.df))
#And richness
Rich <- rarefy(Spe.df, min(rowSums(Spe.df)))
#And beta diversity
beta <- vegdist(Spe.df, binary=TRUE)
mean(beta)


#****Univariate statistics applied to 'big' data****
#Transform data to relative abundance as opposed to raw sequence counts
Spe.df <- sweep(Spe.df, 1, rowSums(Spe.df), '/')
#Scale data with Hellinger to account for zero species paradox
Spe.df <- decostand(Spe.df, 'hel')
#Create empty vector, Z, to store P values
Z <- data.frame("P.values" = rep(NA, 4790))
#Create categorical variables to be tested: here I am testing against three
#Types of wetlands, Quistococha (QU), San Jorge (SJ) and Buena Vista (BV)
df[,1]
Wetlands <- c("QU", "SJ", "SJ", "SJ", "SJ", "QU", "QU", "QU", "BV", "BV", "BV", "BV")
Wetlands <- factor(Wetlands)
#Loop to test which species differ between wetlands
for(i in 1:4790){
  res <- aov(Spe.df[, i] ~ Wetlands)
  Z$P.values[i] <- summary(res)[[1]][["Pr(>F)"]]
  Z$Bonferonni = p.adjust(Z$P.values, method = "bonferroni")
  write.table(Z, "/Users/Damien/Documents/R/Stats_workshop_Juin17/otus/AOV_Pvals_scale.txt", sep = "\t", row.names = F)  
} 
#Output is as follows: col 1 = aov P values; col 2 = P vals adjusted for
#False discovery rate with the more stringent Bonferroni method.
#Remember! At P = 0.05, that's still 1/20 chance of false discovery!
#With 4790 OTUs, we can expect ~ 240 'significant differences' to be false!
#Take subset of Z
sub_Z <- as.data.frame(Z[,2])
#Transpose OTU table and set names
Spe.df_t <- t(Spe.df)
Spe.df_Pvals <- cbind(Spe.df_t, sub_Z)
colnames(Spe.df_Pvals) <- Wetlands
#Sort species with significant differences based on Bon adjusted P values
#And generate heatmap on hellinger transformed rel abundance species data
B.df <- subset(Spe.df_Pvals, Spe.df_Pvals[,13] < 0.05)
write.table(B.df, ".../AOV_otu_table.txt", sep = "\t")
B.df <- read.table(".../AOV_otu_table.txt", sep = "\t", header = T)
B.df <- head(arrange(B.df, desc(Pvals)), n = 50)
rows <- B.df[,1]
B.df <- B.df[, 2:13]
rownames(B.df) <- rows
B.df <- data.matrix(B.df)
my_palette <- colorRampPalette(c("red", "orange", "yellow")) (n = 299)
heatmap.2(B.df, col = my_palette, margins = c(5, 22), dendrogram = "column", density.info = "none", trace = "none")
#Loop over lines 52-57 to increase or reduce stringency as necessary
#Primary different OTUs suggest BV more unique than SJ and QU, which is
#A result of specific Acidobacteria, Nitrospirae and Chloroflexi


#****Multivariate statistics applied to 'big' data****

#Using ranked order of (dis)similarity methods: NMDS and ANOSIM
#NMDS
Spe.df <- df[, 2:4791]
rows <- df[,1]
rownames(Spe.df) <- rows
Spe.df <- sweep(Spe.df, 1, rowSums(Spe.df), '/')
Spe.df.mds <- metaMDS(Spe.df, dist = "bray")
ordiplot(Spe.df.mds, display = "sites", type = "text")
ordiellipse(Spe.df.mds, Wetlands, conf = 0.95, label = TRUE, lty = 2)

#ANOSIM
dist.Spe.df <- dist(Spe.df, method = "euclidean")
anosim.res <- anosim(dist.Spe.df, Wetlands, permutations = 1000, distance = "bray")
summary(anosim.res)
#Summary: Note dissimilarity ranks between and within classes:
#See that QU and SJ display greater within class dissimilarity than BV?
#P = 0.002, suggesting Wetlands are different. 

#The Fisherian approach: Eigen coupled with pseudo F-statistic
#PCA
Spe.hel <- decostand(Spe.df, 'hel')
pca.res <- rda(Spe.hel, scale = TRUE)
head(summary(pca.res), n = 10)
colvec <- c("blue", "red", "red", "red", "red", "blue", "blue", "blue", "green", "green", "green", "green")
biplot(pca.res, type = "n", ylab = "PC2 (13.25%)", xlab = "PC1 (34.14%)")
with(df, points(pca.res, display = "sites", col = colvec, pch = 16))
#text(pca.res, display = "species", cex = 1, font = 2, col = "black")
with(df, ordiellipse(pca.res, Wetlands, kind = "sd", conf = 0.95, label = FALSE, col = c("gray"), lty = c(2)))
#xStart <- c(0,0)
#yStart <- c(0,0)
#P1 <- scores(pca.res, choices=c(1), display=c("species"))
#P2 <- scores(pca.res, choices=c(2), display=c("species"))
#arrows(xStart, yStart, x1 = P1, y1 = P2, length = 0.25, angle = 10, code = 2, col = "gray", lty = 1)
legend(16, 10, pch=c(16), col = c("blue", "red", "green"), legend =c("QU", "SJ", "BV"), bty = "n")

#perMANOVA
Spe.dm <- data.matrix(Spe.hel, rownames.force=NA)
permanova.res <- adonis(Spe.dm ~ Wetlands, permutations = 1000, method = "bray")
print(permanova.res$aov.tab)
#Note that P value is lower (P = 0.001) but null hypothesis still rejected
#Why is the P value lower here? perMANOVA is more sensitive to data
#Transformations than ANOSIM, so if we're following solid theory (i.e.
#To use Hellinger) then we can improve our inference.

#So which approach is better?
#Pros for rank order:
#a) less sensitive to data transformations because compares ranked order
#of dissimilarity within/between groups to the mean dissimilarity
#Cons:
#a) less sensitive to data transformations! 
#b) more sensitive to heteroskedasticity than perMANOVA 

#Pros and cons for PCA/perMANOVA approach
#Pros:
#a) Using a weighted ratio of dissimilarity of within vs between groups
#which makes it more robust to certain types of experimental design
#b) if using correct data transformations can give better results
#c) As eigenvectors retain relationship of separation in space (now 
#measured as eigenvalues) we can actually go on and do a lot more with
#the PCA, such as Principle Components Regression Analysis for 
#modelling variability explained by certain explanatory variables
#Cons:
#a) whereas ANOSIM results can easily be directly compared to other studies
#That have the same degrees of freedom (remember in ANOSIM the dissim-
#ilarity of the points are all just relative to each other), perMANOVA
#can't be because it's dependent on distances derived from the spe counts

#For further multivariate analyses, if sound environmental chemistry
#Data in addition to species data, consider redundancy analysis which
#Looks for a linear relationship between environment ~ biology. As the
#sum of squares can be calculated from this relationship, we can even
#get an R2 value. Good if you want to demonstrate sound relationship
#Between community composition and environmental conditions


#****Machine learning approaches to ecology data****

#Boosted regression to explain % variability of any OTU to an ecosystem
#Process
library("gbm")
library("grid")

df <- read.table(".../otu_table_phylum.txt", sep = "\t", header = T)

#Make sure df is good
Spe.df <- df[, 2:65]
rows <- df[,1]
rownames(Spe.df) <- rows
Spe.df <- sweep(Spe.df, 1, rowSums(Spe.df), '/')
#Create vector of biological function: in this case, I am going to use
#Dissolved CH4 (%) from three wetland sites, up to 1 m deep, every 10 cm:
CH4 <- c(2.91, 0.89, 0.79, 2.03, 2.71, 0.01, 1.35, 0.4, 1.99, 0.34, 0.11, 0.53, 0.11, 0.12, 0.32, 0.27, 0.43, 0.33, 2.42, 3.86, 2.49, 1.51, 1.28, 0.12, 2.48, 0.3, 1.2, 0.6)
Spe.df.ch4 <- cbind(Spe.df, CH4)
#Take a subsample of dataset to split into training vs test
smp <- floor(0.5 * nrow(Spe.df.ch4))
#Create boosting function
tmp1 <- function(Spe.df.ch4){
  n.trees = seq(from = 1, to = 5000)
  train_ind <- sample(seq_len(nrow(Spe.df.ch4)), size = smp)
  train <- Spe.df.ch4[train_ind, ]
  test <- Spe.df.ch4[-train_ind, ]
  boost = gbm(CH4 ~ .-CH4, data = train, distribution = "gaussian", n.trees = 5000, shrinkage = 0.01, interaction.depth = 2, n.minobsinnode = 2)
  predmat = predict(boost, newdata = test, distribution = "gaussian", n.trees = n.trees, shrinkage = 0.01, interaction.depth = 2, n.minobsinnode = 2)
  berr = with(Spe.df.ch4[-train_ind,], apply((predmat - CH4)^2, 2, mean))
}
#Run boosting with ID = 2, nminobs = 2, trees = 5000
tmp1.results <- replicate(1000, tmp1(Spe.df.ch4))
#Plot how the cross validation looks
x1 <- c(1:5000)
y1 <- rowMeans(tmp1.results[ , 1:1000])
z1 <- apply(tmp1.results[ , 1:1000], 2, sd)
plot(x1, y1, ylab = "Mean Squared Error", xlab = "No. of Trees", cex.lab = 1.2, las = 1, main = "ID 2")
segments(x1, (y1 - (z1*0.25)), x1, (y1 + (z1*0.25)))
#Note how the tail of the MSE drops, and then increases with number of
#Trees due to over-fitting of the gbm. For this example, 500 trees is
#A good number. Could try further optimising with ID of 1, 4, 8.

#With cross-validated model, run accordingly
mod1 <- function(Spe.df.ch4){
  boost = gbm(CH4 ~ . -CH4, data = Spe.df.ch4, distribution = "gaussian", n.trees = 500, shrinkage = 0.01, interaction.depth = 2, n.minobsinnode = 2)
  summary(boost)
}
mod1.results <- replicate(1000, mod1(Spe.df.ch4))
#This generates 1000 replicates of the gbm
#Now must do some formatting magic to make sense of the gbm output
for(i in 1:1000){ ifelse(i==1, mod1.res <- tapply(mod1.results[[2*i]], mod1.results[[2*(i-1)+1]], I), mod1.res <- rbind(mod1.res, tapply(mod1.results[[2*i]], mod1.results[[2*(i-1)+1]], I)))}
df <- melt(mod1.res)
df.m.s <- ddply(df, .(Var2), summarise, mean = round(mean(value), 2), sd = round(sd(value), 2))
df.m.s.2 <- data.frame(Var2 = rownames(df.m.s), df.m.s, row.names = NULL)
df.m.s.2$Var2 <- factor(df.m.s$Var2, levels = df.m.s[order(df.m.s$mean), "Var2"])

bp <- ggplot(df.m.s.2, aes(y = Var2, x = mean)) +
  geom_point(stat = "identity") +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd, width = 0.1), position = position_dodge(0))
bp + theme(axis.title.x = element_text(size= 16, vjust = -0.15, face = "bold"), axis.title.y = element_text(size = 16, vjust = 0.35, face = "bold"), axis.text.x = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 10.5, face = "bold"), plot.margin = unit(c(0.8, 0.95, 0.95, 0.95), "cm")) +
  scale_y_discrete(name = "Variable") +
  scale_x_continuous(name = "Relative Influence (%)")

#Thus, we see that the phyla associated with CH4 (%) are
#Archaea within the Parvarchaeota and Euryarchaeota.



