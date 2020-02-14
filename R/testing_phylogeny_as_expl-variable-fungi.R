library(vegan)
library(ape)
library(phangorn)



#read in Mapping file with the Details for each Sample
map <- read.csv("/home/robert/Projects/microbiome-a.incompertus/data/Metadata_with_Host.tsv",sep="\t")

names(map)[1]<-"SampleID"

#read in OTU table which i transposed in Python before

OTU <- read.csv("/home/robert/Projects/microbiome-a.incompertus/data/OTU_table_97-fungi.csv", sep=",", header=TRUE)

#OTU table has a weird column name, fix that
names(OTU)[1]<-"SampleID"

#rearrange rows in mapping file to match relative abundance matrix

tab <- map[match(OTU$SampleID, as.character(map$SampleID)),]
summary(OTU$SampleID==tab$SampleID)

#merge the two Dataframes

fungi <- merge(tab,OTU, by.x="SampleID", by.y = "SampleID")
scaledOTU <- scale(fungi[1:nrow(fungi),24:ncol(fungi)], scale=T)
summary(scaledOTU, display=NULL)
mean(scaledOTU)

#delete the Mock-Community and negatives
fungi_new <- fungi[-c(63,64,65,66),]
fungi_new$Location<-factor(fungi_new$Location)
fungi_new$Gender<-factor(fungi_new$Gender)
fungi_new$Raffaelea<-factor(fungi_new$Raffaelea)


#Beetle Phylogeny might be importend, because each tree can be a different beetle haplotype. Lets test that.

#import the phylogeny for tree

##################################################################  Maximum-Liklihood tree ##################################################################################

tree.ML <- read.nexus('/home/robert/Projects/microbiome-a.incompertus/scripts/phylogenetic_analysis_beetle/without_outgroup/Ai-COI-137_samples_24Jan2019_aligned_without_outgroup.nxs')


#plot tree
plot(tree.ML)


#The tree and alignment have less tips/entries than the dataframe has rows. How many?

differ.ML<- setdiff(x=fungi_new$Sample_ID, y=tree.ML$tip.label)

#now select only the samples i have in the tree from the fungi dataframe

fungi.beetle.phylogeny <- fungi_new[fungi_new$Sample_ID %in% tree.ML$tip.label,]

#Match the dataframe to the tree$tip.labels

matched.fungi <- fungi.beetle.phylogeny[match(as.character(tree.ML$tip.label), fungi.beetle.phylogeny$Sample_ID),]

#Lest check if it worked
all(tree.ML$tip.label == matched.fungi$Sample_ID)


#OTUs for fungi.beetle.phylogeny

OTU.beetle.phylogeny <- matched.fungi[,24:ncol(matched.fungi)]


normi.phylogeny <- decostand(OTU.beetle.phylogeny, method="hellinger")

#Get the environmental variable I thnk are a key factor

ambrosia.phylogeny.env <- matched.fungi[, c('Location', 'Gender', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

#ambrosia.phylogeny.env$Tree <- factor(with(ambrosia.phylogeny.env, paste(Location, Tree, sep ='_')))
ambrosia.phylogeny.env <- droplevels(ambrosia.phylogeny.env)

summary(ambrosia.phylogeny.env)



# The column collected should be a factor
ambrosia.phylogeny.env$Collected <- factor(ambrosia.phylogeny.env$Collected)


# set up full and null models for ordistep
fungi.phylogeny.cap1 <- capscale(normi.phylogeny ~ ., data=ambrosia.phylogeny.env, dist="bray")
fungi.phylogeny.cap0 <- capscale(normi.phylogeny ~ 1, data=ambrosia.phylogeny.env, dist="bray")


#perform forward and backward selection of explanatory variables 

step.phylogeny.env <- ordistep(fungi.phylogeny.cap0, scope=formula(fungi.phylogeny.cap1))

step.phylogeny.env$anova


#Pairwise Distance from Phylogentic Tree

beetle.distance <- as.dist(cophenetic.phylo(tree.ML))

#compute distance on alignment

#beetle.distance2 <- as.dist(distance.fasta.Ai)



beetle.distance.pco0<-capscale(beetle.distance~1)  # principle coordinates matrix is in scores(tree.pco,1:8,'sites')
#beetle.distance2.pco0<-capscale(beetle.distance2~1


#Check the structure of the distance and find the coordinate matrix. 


#summary(beetle.distance2.pco0)

#Get the scores

beetle.distance.scores<-as.data.frame(scores(beetle.distance.pco0, 1:9,'sites'))
#beetle.distance.scores2<-as.data.frame(scores(beetle.distance2.pco0,1:6,'sites'))


# set up full and null models for ordistep
fungi.beetle.cap1 <- capscale(normi.phylogeny ~ ., data=beetle.distance.scores, dist="bray")
fungi.beetle.cap0 <- capscale(normi.phylogeny ~ 1, data=beetle.distance.scores, dist="bray")


step.phylogeny <- ordistep(fungi.beetle.cap0, scope = formula(fungi.beetle.cap1))

step.phylogeny$anova

#transform the phylogenetic distance

beetle.pcnm <- as.data.frame(scores(pcnm(dist(beetle.distance))))

dim(beetle.pcnm)


#Select the OTUs from the matched.fungi dataframe

ambrosia.var.phylogeny <- varpart(OTU.beetle.phylogeny, beetle.distance.scores[, paste('MDS', 1:7, sep='')], ~Location, ~Host, data=matched.fungi[,1:21])

plot(ambrosia.var.phylogeny, bg=1:3, Xnames=c('Phylogeny', 'Location', 'Host'))

############################################################## coalescent tree ##########################################################################################

#Import the Alignment
Ai.phy <- read.phyDat('/home/robert/Projects/microbiome-a.incompertus/scripts/phylogenetic_analysis_beetle/without_outgroup/Ai-COI-137_samples_24Jan2019_aligned_without_outgroup.phy', format="phylip")

Ai.fasta <- read.FASTA('/home/robert/Projects/microbiome-a.incompertus/scripts/phylogenetic_analysis_beetle/without_outgroup/Ai-COI-137_samples_24Jan2019_aligned_without_outgroup.fasta')


#out <- coalescentMCMC(Ai.fasta, ntrees = 100000, burnin = 1000)

#plot(out)

#dim(out)

#colnames(out)

#R <- getMCMCtrees(1)

#getMCMCstats()


#plot(R)

#Compute the distance matrix for the alignment
distance.fasta.Ai <- dist.ml(Ai.phy)


#Generate the UPGMA tree
treeUPGMA <- upgma(distance.fasta.Ai)

#Plot the tree

plot(treeUPGMA, main="UPGMA")

#The tree and alignment have less tips/entries than the dataframe has rows. How many?

differ.UPGMA<- setdiff(x=fungi_new$Sample_ID, y=treeUPGMA$tip.label)

#now select only the samples i have in the tree from the fungi dataframe

fungi.beetle.phylogeny.UPGMA <- fungi_new[fungi_new$Sample_ID %in% treeUPGMA$tip.label,]

#Match the dataframe to the tree$tip.labels

matched.fungi.UPGMA <- fungi.beetle.phylogeny.UPGMA[match(as.character(treeUPGMA$tip.label), fungi.beetle.phylogeny.UPGMA$Sample_ID),]

#Lest check if it worked
all(treeUPGMA$tip.label == matched.fungi.UPGMA$Sample_ID)


#OTUs for fungi.beetle.phylogeny

OTU.beetle.phylogeny.UPGMA <- matched.fungi.UPGMA[,24:ncol(matched.fungi.UPGMA)]


normi.phylogeny.UPGMA <- decostand(OTU.beetle.phylogeny.UPGMA, method="hellinger")

#Get the environmental variable I thnk are a key factor

ambrosia.phylogeny.env.UPGMA <- matched.fungi.UPGMA[, c('Location', 'Gender', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

#ambrosia.phylogeny.env$Tree <- factor(with(ambrosia.phylogeny.env, paste(Location, Tree, sep ='_')))
ambrosia.phylogeny.env.UPGMA <- droplevels(ambrosia.phylogeny.env.UPGMA)

summary(ambrosia.phylogeny.env.UPGMA)



# The column collected should be a factor
ambrosia.phylogeny.env.UPGMA$Collected <- factor(ambrosia.phylogeny.env.UPGMA$Collected)


# set up full and null models for ordistep
fungi.phylogeny.cap1.UPGMA <- capscale(normi.phylogeny.UPGMA ~ ., data=ambrosia.phylogeny.env.UPGMA, dist="bray")
fungi.phylogeny.cap0.UPGMA <- capscale(normi.phylogeny.UPGMA ~ 1, data=ambrosia.phylogeny.env.UPGMA, dist="bray")


#perform forward and backward selection of explanatory variables 

step.phylogeny.env.UPGMA <- ordistep(fungi.phylogeny.cap0.UPGMA, scope=formula(fungi.phylogeny.cap1.UPGMA))

step.phylogeny.env.UPGMA$anova


#Pairwise Distance from Phylogentic Tree

beetle.distance.UPGMA <- as.dist(cophenetic.phylo(treeUPGMA))

#compute distance on alignment

#beetle.distance2 <- as.dist(distance.fasta.Ai)



beetle.distance.pco0.UPGMA<-capscale(beetle.distance.UPGMA~1)  # principle coordinates matrix is in scores(tree.pco,1:8,'sites')
#beetle.distance2.pco0<-capscale(beetle.distance2~1)



#Check the structure of the distance and find the coordinate matrix. 


summary(beetle.distance.pco0.UPGMA)

#Get the scores

beetle.distance.scores.UPGMA<-as.data.frame(scores(beetle.distance.pco0.UPGMA,1:6,'sites'))
#beetle.distance.scores2<-as.data.frame(scores(beetle.distance2.pco0,1:6,'sites'))


# set up full and null models for ordistep
fungi.beetle.cap1.UPGMA <- capscale(normi.phylogeny.UPGMA ~ ., data=beetle.distance.scores.UPGMA, dist="bray")
fungi.beetle.cap0.UPGMA <- capscale(normi.phylogeny.UPGMA ~ 1, data=beetle.distance.scores.UPGMA, dist="bray")


step.phylogeny.UPGMA <- ordistep(fungi.beetle.cap0.UPGMA, scope = formula(fungi.beetle.cap1.UPGMA))

step.phylogeny.UPGMA$anova


#transform the phylogenetic distance

beetle.pcnm.UPGMA <- as.data.frame(scores(pcnm(dist(beetle.distance.UPGMA))))

dim(beetle.pcnm.UPGMA)


#Select the OTUs from the matched.fungi dataframe

ambrosia.var.phylogeny.UPGMA <- varpart(OTU.beetle.phylogeny.UPGMA, beetle.distance.scores.UPGMA[, paste('MDS', 1:6, sep='')], ~Location, ~Host, data=matched.fungi.UPGMA[,1:21])

plot(ambrosia.var.phylogeny.UPGMA, bg=1:3, Xnames=c('Phylogeny', 'Location', 'Host'))

#####################################################################################################################################################################################

#test for the significance of phylogeny 

sig.rda.phylogeny <- rda(OTU.beetle.phylogeny, beetle.distance.scores[, paste('MDS', 1:6, sep='')], cbind(ambrosia.var.phylogeny$Location, ambrosia.var.phylogeny$Host))

anova(sig.rda.phylogeny)

plot(sig.rda.phylogeny)


k <- scores(sig.rda.phylogeny, display = "species")


l <- scores(sig.rda.phylogeny, display = "sites")

head(k[order(k[, 1], decreasing=T), ])

head(k[order(k[, 1], decreasing=T), ], 10)
head(k[order(k[, 1], decreasing=F), ], 10)


head(l[order(l[, 1], decreasing=T), ])


layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
pdf('test.pdf', height=15, width=10)
for(i in 1:6) {
  plot(treeUPGMA,label.offset=0.001)
  tiplabels(pch=21,cex=3*(beetle.distance.scores.UPGMA[,i]-min(beetle.distance.scores.UPGMA[,i])))
  text(0.005,5,paste('phyloMDS',i,sep=''))
}
dev.off()

############################################################################################################################################################################################
#remove sample 136 so sample RM79

#Import the Alignment
Ai.phy <- read.phyDat('/home/robert/Projects/microbiome-a.incompertus/scripts/phylogenetic_analysis_beetle/withoutC21andOutgroup/Ai-COI-137_samples_24Jan2019_aligned_without_outgroup_c21.phy', format="phylip")

#plot(R)

#Compute the distance matrix for the alignment
distance.fasta.Ai <- dist.ml(Ai.phy)


#Generate the UPGMA tree
treeUPGMA <- upgma(distance.fasta.Ai)
#treeML <- read.nexus('/home/robert/Projects/fungal_microbiom_Ai/scripts/97_otu/Beta_Diversity_97/Beetle_Phylogeny/ML/withoutC21/Ai.aligned_June2018withoutC21')
#Plot the tree

plot(treeUPGMA, main="UPGMA")

#The tree and alignment have less tips/entries than the dataframe has rows. How many?

differ.UPGMA<- setdiff(x=fungi_new$Sample_ID, y=treeUPGMA$tip.label)

#now select only the samples i have in the tree from the fungi dataframe

fungi.beetle.phylogeny.UPGMA <- fungi_new[fungi_new$Sample_ID %in% treeUPGMA$tip.label,]

#Match the dataframe to the tree$tip.labels

matched.fungi.UPGMA <- fungi.beetle.phylogeny.UPGMA[match(as.character(treeUPGMA$tip.label), fungi.beetle.phylogeny.UPGMA$Sample_ID),]

#Lest check if it worked
all(treeUPGMA$tip.label == matched.fungi.UPGMA$Sample_ID)


#remove the sample 136
new.matched.UPGMA <- matched.fungi.UPGMA[!matched.fungi.UPGMA$SampleID=="RM79", ]


#OTUs for fungi.beetle.phylogeny

OTU.beetle.phylogeny.UPGMA.new <- new.matched.UPGMA[,24:ncol(new.matched.UPGMA)]


normi.phylogeny.UPGMA.new <- decostand(OTU.beetle.phylogeny.UPGMA.new, method="hellinger")

#Get the environmental variable I thnk are a key factor

ambrosia.phylogeny.env.UPGMA.new <- new.matched.UPGMA[, c('Location', 'Gender', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

#ambrosia.phylogeny.env$Tree <- factor(with(ambrosia.phylogeny.env, paste(Location, Tree, sep ='_')))
ambrosia.phylogeny.env.UPGMA.new <- droplevels(ambrosia.phylogeny.env.UPGMA.new)

summary(ambrosia.phylogeny.env.UPGMA.new)



# The column collected should be a factor
ambrosia.phylogeny.env.UPGMA.new$Collected <- factor(ambrosia.phylogeny.env.UPGMA.new$Collected)


# set up full and null models for ordistep
fungi.phylogeny.cap1.UPGMA.new <- capscale(normi.phylogeny.UPGMA.new ~ ., data=ambrosia.phylogeny.env.UPGMA.new, dist="bray")
fungi.phylogeny.cap0.UPGMA.new <- capscale(normi.phylogeny.UPGMA.new ~ 1, data=ambrosia.phylogeny.env.UPGMA.new, dist="bray")


#perform forward and backward selection of explanatory variables 

step.phylogeny.env.UPGMA.new <- ordistep(fungi.phylogeny.cap0.UPGMA.new, scope=formula(fungi.phylogeny.cap1.UPGMA.new))

step.phylogeny.env.UPGMA.new$anova


#Pairwise Distance from Phylogentic Tree

beetle.distance.UPGMA.new <- as.dist(cophenetic.phylo(treeUPGMA))

#compute distance on alignment

#beetle.distance2 <- as.dist(distance.fasta.Ai)



beetle.distance.pco0.UPGMA.new<-capscale(beetle.distance.UPGMA.new~1)  # principle coordinates matrix is in scores(tree.pco,1:8,'sites')
#beetle.distance2.pco0<-capscale(beetle.distance2~1)



#Check the structure of the distance and find the coordinate matrix. 


summary(beetle.distance.pco0.UPGMA.new)

#Get the scores

beetle.distance.scores.UPGMA.new<-as.data.frame(scores(beetle.distance.pco0.UPGMA.new,1:6,'sites'))
#beetle.distance.scores2<-as.data.frame(scores(beetle.distance2.pco0,1:6,'sites'))


# set up full and null models for ordistep
fungi.beetle.cap1.UPGMA.new <- capscale(normi.phylogeny.UPGMA.new ~ ., data=beetle.distance.scores.UPGMA.new, dist="bray")
fungi.beetle.cap0.UPGMA.new <- capscale(normi.phylogeny.UPGMA.new ~ 1, data=beetle.distance.scores.UPGMA.new, dist="bray")


step.phylogeny.UPGMA.new <- ordistep(fungi.beetle.cap0.UPGMA.new, scope = formula(fungi.beetle.cap1.UPGMA.new))

step.phylogeny.UPGMA.new$anova


#transform the phylogenetic distance

beetle.pcnm.UPGMA.new <- as.data.frame(scores(pcnm(dist(beetle.distance.UPGMA.new))))

dim(beetle.pcnm.UPGMA.new)


#Select the OTUs from the matched.fungi dataframe

ambrosia.var.phylogeny.UPGMA.new <- varpart(OTU.beetle.phylogeny.UPGMA.new, beetle.distance.scores.UPGMA.new[, paste('MDS', 1:6, sep='')], ~Location, ~Host, data=matched.fungi.UPGMA[,1:21])

plot(ambrosia.var.phylogeny.UPGMA.new, bg=1:3, Xnames=c('Phylogeny', 'Location', 'Host'))

ambrosia.var.phylogeny.UPGMA.new
#test for the significance of phylogeny 

sig.rda.phylogeny.new <- rda(OTU.beetle.phylogeny.UPGMA.new, beetle.distance.scores.UPGMA.new[, paste('MDS', 1:6, sep='')], cbind(ambrosia.var.phylogeny.UPGMA.new$Location, ambrosia.var.phylogeny.UPGMA.new$Host))

anova(sig.rda.phylogeny.new)

plot(sig.rda.phylogeny.new)


k <- scores(sig.rda.phylogeny.new, display = "species")


l <- scores(sig.rda.phylogeny.new, display = "sites")

head(k[order(k[, 1], decreasing=T), ])

head(k[order(k[, 1], decreasing=T), ], 10)
head(k[order(k[, 1], decreasing=F), ], 10)


head(l[order(l[, 1], decreasing=T), ])


layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
pdf('test.pdf', height=15, width=10)
for(i in 1:6) {
  plot(treeUPGMA,label.offset=0.001)
  tiplabels(pch=21,cex=3*(beetle.distance.scores.UPGMA[,i]-min(beetle.distance.scores.UPGMA[,i])))
  text(0.005,5,paste('phyloMDS',i,sep=''))
}
dev.off()


