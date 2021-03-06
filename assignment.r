###############################
######## Assignment 8 #########
###############################

# Measuring pairwise dissimilarities

rm(list = ls())
#setwd("/Users/emanuelastruffolino/Desktop/SequenceCourse/Lezione 8_22102012")
#getwd()
library (TraMineRextras)

#1. Re-create the state sequence object biofam.seq and the matrix 
   #dOM of pairwise OM dissimilarities based on the properties matrix 
   #considered in the previous assignment.

data(biofam)
mycol<-brewer.pal(,"RdBu")
biofam.lab<-c("Parent", "Left", "Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced")
biofam.shortlab<-c("P", "L", "M", "LM", "C","LC","LMC", "D")
biofam.seq <- seqdef(biofam[,10:25],states=biofam.shortlab,labels=biofam.lab)
summary (biofamseq)

properties <- matrix(c(# left, married, child, divorced
  0, 0, 0, 0,  # parent
  1, 0, 0, 0,  # left
  0, 1, .5, 0, # marr
  1, 1, 0, 0,  # left+marr
  0, 0, 1, 0,  # child
  1, 0, 1, 0,  # left+child
  1, 1, 1, 0,  # left+marr+child
  .5, 1, .5, 1 # divorced
), 8, 4, byrow=TRUE)
scm <- as.matrix(dist(properties))
indel <- .5*max(scm)
distOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=scm, full.matrix = FALSE)

#2. Create the hierarchical cluster tree object once with the Ward method 
    #and once with WPGMA (McQuitty) method, and using the sequence weights 
    #in each case (Tip: retrieve the weights with attr(biofam.seq, "weight").)

weight<-attr(biofam.seq, "weight")
clust.ward<-hclust(as.dist(distOM), method="ward", members=weight)
clust.WPGMA<-hclust(as.dist(distOM), method="mcquitty", members=weight)

#3. Display both hierarchical trees side by side.
par(mfrow=c(1,2))
plot(clust.ward, labels=FALSE, main="Ward")
plot(clust.WPGMA, labels=FALSE, main="mcquitty")

#4. Select the three-cluster solution from the Ward analysis, and label the 
    #clusters by looking at the I-plots by cluster.
biofam.cl3h <- cutree(clust.ward, k = 3)
biofam.cl3h[1:10]
seqIplot(biofam.seq, group = group.p(biofam.cl3h), border = NA,  sortv = "from.start")

cl3h.labels <- c("Parent-LeftMarried-Child", "EarlyLeft",
                 "Belated left-Married")
#biofam.cl3h.factor <- factor(biofam.cl3h, levels = c(1, 2, 3), labels = cl3h.labels)
seqIplot(biofam.seq, group = biofam.cl3h, border = NA,  sortv = "from.start")

#5. Make and comment the silhouette plot of the retained solution. 
    #(Tip: Use the silhouette function of the cluster package.)
library (cluster)
plot(silhouette(biofam.cl3h, dmatrix=as.matrix(distOM)))

#6. Look at the partition quality measures returned for the Ward-three-cluster 
    #solution by the wcClusterQuality function (library WeightedCluster) when 
    #specifying the sequence weights. Compare the average silhouette with the 
    #value obtained with summary(silhouette(...)) which does not account for 
    #weights.
library (WeightedCluster)
wcClusterQuality(distOM, biofam.cl3h, weight = weight)
summary(silhouette(biofam.cl3h, dmatrix=as.matrix(distOM)))

#Warning in install.packages: package ‘WeightedCluster’ is not available (for R version 2.15.1)

#7. Study with logistic regressions how the cluster membership is related to 
    #the sex, birthyear, and the language of the questionnaire.

cluster1 <- as.numeric(biofam.cl3h==1)
reg1 <- glm(cluster1 ~ sex + birthyr + plingu02, data=biofam, family=binomial)
summary(reg1)

cluster2 <- as.numeric(biofam.cl3h==2)
reg2 <- glm(cluster2 ~ sex + birthyr + plingu02, data=biofam, family=binomial)
summary(reg2)

cluster3 <- as.numeric(biofam.cl3h==3)
reg3 <- glm(cluster3 ~ sex + birthyr + plingu02, data=biofam, family=binomial)
summary(reg3)

#8. In order to make a PAM clustering, first examine the evolution (values  
    #and plot) of the quality indicators for the solutions for 
    #k = 2, . . . , 20. What value k do you retain?
quality <- wcKMedRange(distOM, 1:10, weight = weight)
plot (quality)
#I can't check because I can't install the WeightedCluster package

#9. Label the clusters of the PAM solution for the chosen k and plot the 
    #mean time spent in the states by cluster.
