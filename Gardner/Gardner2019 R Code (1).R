#8/26/21
library(phyloseq)
library(vegan)
library(ggplot2)
library(Biostrings)
library(DECIPHER)
library(dplyr)
library(speedyseq)
library(MASS)

getwd()
setwd("/Users/jchung800/Desktop/R/Gardner")
      
Gardner2019.seqtab <- read.table("Gardner2019.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Gardner2019.taxonomy <- read.table ("Gardner2019 taxonomy.tabular",header=TRUE,row.names=1,sep="\t")

colnames(Gardner2019.taxonomy)
colnames(Gardner2019.seqtab) #SRR numbers

Gardner2019.metadata<-read.csv("Gardner2019.metadata.csv",header=TRUE,row.names=1)
rownames(Gardner2019.metadata) #SRR numbers
length(colnames(Gardner2019.seqtab)) #38
length(rownames(Gardner2019.metadata)) #38

#make phyloseq
OTU <- otu_table(Gardner2019.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Gardner2019.metadata)
TAX <- tax_table(as.matrix(Gardner2019.taxonomy))
Gardner2019.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Gardner2019.ps 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 2555 taxa and 38 samples ]:
#   sample_data() Sample Data:        [ 38 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 2555 taxa by 7 taxonomic ranks ]:
#   taxa are rows
saveRDS(Gardner2019.ps, "Gardner2019.ps")
hist(taxa_sums(Gardner2019.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Gardner2019.ps)
Gardner.filter.ps<-prune_samples(sample_sums(Gardner2019.ps) > 999, Gardner2019.ps) # prune with low sample number
Gardner.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 2555 taxa and 34 samples ]:
#   sample_data() Sample Data:        [ 34 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 2555 taxa by 7 taxonomic ranks ]:
#   taxa are rows
Gardner.filter.ps<-prune_taxa(taxa_sums(Gardner.filter.ps) > 49, Gardner.filter.ps)
Gardner.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 1156 taxa and 34 samples ]:
#   sample_data() Sample Data:        [ 34 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 1156 taxa by 7 taxonomic ranks ]:
#   taxa are rows
saveRDS(Gardner.filter.ps, "Gardner.filter.ps")
# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}
#convert to ASV number
library(Biostrings)
Gardner.dna <- Biostrings::DNAStringSet(taxa_names(Gardner.filter.ps))
names(Gardner.dna) <- taxa_names(Gardner.filter.ps)
Gardner2019.filter.ps <- merge_phyloseq(Gardner.filter.ps, Gardner.dna)
taxa_names(Gardner2019.filter.ps) <- paste0("Gardner2019.ASV", seq(ntaxa(Gardner2019.filter.ps)))
Gardner2019.filter.ps
rownames (Gardner2019.filter.ps@sam_data)
rownames (Gardner2019.filter.ps@tax_table)

#transform to relative abundance
Gardner.rel.ps = transform_sample_counts(Gardner2019.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Gardner.otu<-psotu2veg(Gardner.rel.ps)
Gardner2019.metadata<-pssd2veg(Gardner.rel.ps)
#species glom 

#Bray-Curtis Dissimilarity
library(vegan)
Gardner.comm.bc.dist <- vegdist(Gardner.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")
saveRDS(Gardner.comm.bc.dist, "Gardner Bacterial Bray- Curtis")

Gardner.rename.otu <- Gardner.otu
rownames(Gardner.rename.otu)
rownames(x) <- newfile
adonis(Gardner.comm.bc.dist ~ lat_lon, data = Gardner2019.metadata, permutations = 99)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# lat_lon   18    6.5119 0.36177  1.1875 0.58764   0.11
# Residuals 15    4.5696 0.30464         0.41236       
# Total     33   11.0815                 1.00000       
adonis(Gardner.comm.bc.dist ~ Host + lat_lon, data = Gardner2019.metadata, permutations = 99)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Host       3    4.0241 1.34138  7.3611 0.36314   0.01 **
# lat_lon   18    4.8707 0.27059  1.4849 0.43953   0.01 **
# Residuals 12    2.1867 0.18223         0.19733          
# Total     33   11.0815                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Gardner.simperOutput<-simper(Gardner.otu,Gardner2019.metadata$lat_lon,perm=99)
Gardner.simperOutput


#subset the Alphaproteobacteria
Gardner.alpha.ps= subset_taxa(Gardner.rel.ps, class=="Alphaproteobacteria") 
Gardner.alpha.ps #731 taxa (28.6% of the taxa)
title = "plot_bar; Gardner Alphaproteobacteria Abundance"
plot_bar(Gardner.alpha.ps, "lat_lon", "Abundance", title=title) #highest abundance in 4.28 S, 55.7 E

#subset the gammaproteobacteria
Gardner.gamma.ps= subset_taxa(Gardner.rel.ps, class=="Gammaproteobacteria") 
Gardner.gamma.ps #665 taxa
title = "plot_bar; Gardner Gammaprotebacteria Abundance"
plot_bar(Gardner.gamma.ps, "lat_lon", "Abundance", title=title) #highest abundance in 4.28 S, 55.7 E very high 

Gardner.actin.ps= subset_taxa(Gardner.rel.ps, class=="Actinobacteria") 
Gardner.actin.ps #57taxa
title = "plot_bar; Gardner Actinobacteria Abundance"
plot_bar(Gardner.actin.ps, "lat_lon", "Abundance", title=title) #highest abundance in 4.28 S, 55.7 E very high 

Gardner.endoz.ps= subset_taxa(Gardner.rel.ps, family=="Endozoicomonadaceae") 
Gardner.endoz.ps #113 taxa
title = "plot_bar; Gardner Endozoicomonadaceae Abundance"
plot_bar(Gardner.endoz.ps, "lat_lon", "Abundance", title=title) #highest abundance in 4.28 S, 55.7 E very high 

Gardner.rhizo.ps= subset_taxa(Gardner.rel.ps, order=="Rhizobiales") 
Gardner.rhizo.ps #111 taxa
title = "plot_bar; Gardner Rhizobiales Abundance"
plot_bar(Gardner.rhizo.ps, "lat_lon", "Abundance", title=title) #highest abundance in 4.307 S, 55.730212 E (very high, aprox 0.37

library (DECIPHER)
library (speedyseq)
Gardner.refseq<- Gardner.rel.ps@refseq
Gardner.refseq
# align
aln<-AlignSeqs(Gardner.refseq,normPower=0, processors=10) #<5 minutes
# distance matrix
dm<-DistanceMatrix(aln,includeTerminalGaps = FALSE,processors = 10) #quick
# find clusters
Gardner.clusters<-IdClusters(dm,method="UPGMA",type="clusters",cutoff=0.01) #quick
# note, no identity errors

# save intermediate objects
saveRDS(aln,"Gardner.alignment.RDS")
saveRDS(dm,"Gardner.dm.RDS")
saveRDS(Gardner.clusters,"Gardner.clusters.RDS")

# how to handle the non-finite values of the distance matrix before UPGMA fails?
Gardner.dm.dataframe<-as.data.frame(dm)
is.na(Gardner.dm.dataframe)<-sapply(Gardner.dm.dataframe, is.infinite)
Gardner.dm.dataframe[is.na(Gardner.dm.dataframe)]<-1
Gardner.dm2<-as.matrix(Gardner.dm.dataframe)
saveRDS(Gardner.dm2,"Gardner.dm2.RDS")

#Making a tree
library(phangorn)
Gardner.upgma<-upgma(Gardner.dm2)
plot(Gardner.upgma,type="fan",cex=0.8,show.tip.label=FALSE)
saveRDS(Gardner.upgma,"Gardner.upgma.tree.RDS")


#9/1
#Diversity Tests
Gardner.alphadiv <- diversity(Gardner.otu,MARGIN=2, index="invsimpson")
hist(Gardner.alphadiv)
library(ggplot2)
ggplot(Gardner2019.metadata,aes(x=Host, y=alpha)) + geom_boxplot()
Gardner.anova <- aov(alpha~Host,Gardner2019.metadata)
summary(Gardner.anova)
library(agricolae)
tukey_result <- HSD.test(Gardner.anova, "Host", group=TRUE)
print(tukey_result)

library (phyloseq)
plot_richness(Gardner2019.metadata, color="Type",x="Host", taxa_are_rows=TRUE)
plot_richness(OTU,x="Coral", color="Host", measures=c("Chao1","Shannon"))
plot_richness(Gardner2019.metadata)


#9/12
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
meta <- data.frame (phyloseq::Yu2020.metadata(Yu2020.rel.ps))
colorCode <- c()

#9/13
plot_tree(Gardner2019.filter.ps, color="SampleType", label.tips="Phylum", ladderize= "left", justify = "left", size = "Abundance")

phyloseq::plot_bar(Gardner.rel.ps, fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Host, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H

#mantel test
Gardner.Bray.distmat <- as.matrix(Gardner.comm.bc.dist)
GardnerITS.Bray.distmat <- as.matrix(GardnerITS.comm.bc.dist)
mantel(Gardner.Bray.distmat,GardnerITS.Bray.distmat)
bac.matrix <- Gardner.Bray.distmat
its.matrix <- GardnerITS.Bray.distmat
dim(bac.matrix)
dim(its.matrix)
rownames(bac.matrix)
rownames(its.matrix)

#different number of samples, so need to make all these the same sample names for both the Bacteria.OTU file and the ITS.OTU file

rownames(Gardner2019.metadata)
rownames(Gardner.otu)
newrownames<-Gardner2019.metadata$Sample.Name
newrows2<-gsub("16S","",newrownames)
newrows2
newrows2<-gsub("-","",newrows2)
newrows2
rownames(Gardner.otu)<-newrows2

newrows3<-Gardner2019.ITS.metadata$Sample.Name
newrows3
newrows4<-gsub("ITS2","",newrows3)
newrows4<-gsub("-","",newrows4)
newrows4

rownames(Gardner2019.ITS.metadata)
rownames(Gardner2019.ITS.otu)<-newrows4

subset.Pollock.metadata<-Pollock.metadata[which(row.names(Pollock.metadata) %in% newNames),]
subset.Gardner.ITS.otu<-Gardner2019.ITS.otu[which(row.names(Gardner2019.ITS.otu) %in% row.names(Gardner.otu)),]

subset.Gardner.ITS.otu<-Gardner2019.ITS.otu[which(row.names(Gardner2019.ITS.otu) %in% row.names(Gardner.otu)),]
dim(subset.Gardner.ITS.otu) #33x264
subset.Gardner.bac.otu<-Gardner.otu[which(row.names(Gardner.otu) %in% row.names(subset.Gardner.ITS.otu)),]
dim(subset.Gardner.bac.otu) #33x1156

rownames(subset.Gardner.ITS.otu)
rownames(subset.Gardner.bac.otu)
subset.Gardner.ITS.otu<-subset.Gardner.ITS.otu[rownames(subset.Gardner.bac.otu),]

new.bac.dist<-vegdist(subset.Gardner.bac.otu, method = "bray")
new.ITS.dist<-vegdist(subset.Gardner.ITS.otu, method = "bray")

mantel.test(as.matrix(new.bac.dist),as.matrix(new.ITS.dist), graph=TRUE)
bac.mds<-metaMDS(new.bac.dist)
plot(bac.mds,type="t")
ef<-envfit(bac.mds,subset.Gardner.ITS.otu,permu=999)
ef
write.csv(ef, "Gardner.ef.csv")

plot(bac.mds, display = "sites")
plot(ef, p.max = 0.01)

plot(bac.mds, display = "sites")
plot(ef, p.max = 0.003)

plot(bac.mds, display = "sites")
plot(ef, p.max = 0.0015)

bac.mds$points #put this in ggplot
plot(bac.mds$points)

bac.mds2<-metaMDS(subset.Gardner.bac.otu)
plot(bac.mds2,display="species")
ef2 <- envfit(bac.mds2, subset.Gardner.ITS.otu,permu=999)
plot(ef2, p.max = 0.0015, cex=0.6)

spnames<-colnames(subset.Gardner.bac.otu)

plot1<-plot(bac.mds2,display="species")
identify(plot1, "species", labels=spnames)
## click on the plot where you want labels to be ##
## then hit ESC and wait ###

bac.df<-as.data.frame(subset.Gardner.bac.otu)
its.df<-as.data.frame(subset.Gardner.ITS.otu)

plot(its.df$Gardner2019.ITS.ASV1,bac.df$Gardner2019.ASV1) #no intermediate,if this ITS is present, then this ASV is abundant 
Gardner.rel.ps@tax_table #this is an Actinobacteria (Brevibacterium aurantiacum)

plot(its.df$Gardner2019.ITS.ASV2,bac.df$Gardner2019.ASV2) #no intermediate,if this ITS is present, then this ASV is abundant 
Gardner.rel.ps@tax_table #this is an Actinobacteria (Brevibacterium aurantiacum)

plot(its.df$Gardner2019.ITS.ASV20,bac.df$Gardner2019.ASV10) #no intermediate,if this ITS is present, then this ASV is abundant 


testDF<-as.data.frame(cbind(subset.Gardner.bac.otu[,1:5],subset.Gardner.ITS.otu[,1:5]))
pairs(testDF)
Gardner.rel.ps@tax_table[1:3,]

#to do 
#redo this with glom and do with yu and camp 

Gardner.species.glom <- tax_glom(Gardner.rel.ps, taxrank="species", NArm = FALSE)
Gardner.species.glom

# isolate the sequences
G19ITS.refseq <- Gardner2019.ITS.filter.ps@refseq
# align
aln<-AlignSeqs(G19ITS.refseq,normPower=0, processors=10) #<5 minutes
# distance matrix
dm<-DistanceMatrix(aln,includeTerminalGaps = FALSE,processors = 10) #quick
# find clusters
G19ITS.clusters<-IdClusters(dm,method="UPGMA",type="clusters",cutoff=0.01) #quick
# note, no identity errors

# save intermediate objects
saveRDS(Yu.seq.ps,"Yu.seq.ps.RDS")
saveRDS(aln,"Yu.alignment.RDS")
saveRDS(dm,"Yu.dm.RDS")
saveRDS(Yu.clusters,"Yu.clusters.RDS")

# how to handle the non-finite values of the distance matrix before UPGMA fails?
G19ITS.dm.dataframe<-as.data.frame(dm)
is.na(G19ITS.dm.dataframe)<-sapply(G19ITS.dm.dataframe, is.infinite)
G19ITS.dm.dataframe[is.na(G19ITS.dm.dataframe)]<-1
G19ITS.dm2<-as.matrix(G19ITS.dm.dataframe)
saveRDS(G19ITS.dm2,"G19ITS.dm2.RDS")

#how do we get that into a tree?
install.packages("phangorn")
library(phangorn)
G19ITS.upgma<-upgma(G19ITS.dm2)
plot(G19ITS.upgma,type="fan",cex=0.8,show.tip.label=FALSE)

saveRDS(Yu.upgma,"Yu.upgma.tree.RDS")

#add tree to phyloseq object
G19ITS.tree.ps<-G19ITS.rel.ps
phy_tree(G19ITS.tree.ps)<-G19ITS.upgma
G19ITS.tree.ps



# glom
G19ITS.glom.005<-tree_glom(G19ITS.tree.ps,resolution=0.005)
G19ITS.glom.005
plot(G19ITS.glom.005@phy_tree)

G19ITS.glom.005@tax_table #ITS are not classified by dada2

writeXStringSet(G19ITS.glom.005@refseq, "G19ITS.glom.005.fasta", format="fasta")

library(phytools)
fastDist(Yu2020ITS.new.glom.005@phy_tree,"Yu2020ITS.ASV3","Yu2020ITS.ASV2")
# only 3%



  # need to plot tree of sequences to see if there is a reverse-complement problem
  
  # isolate the sequences
  G19ITS.refseq<-Gardner2019.ITS.filter.ps@refseq
G19ITS.refseq
# align
aln<-AlignSeqs(G19ITS.refseq,normPower=0, processors=10) #<5 minutes
# distance matrix
dm<-DistanceMatrix(aln,includeTerminalGaps = FALSE,processors = 10) #quick
# find clusters
G19ITS.clusters<-IdClusters(dm,method="UPGMA",type="clusters",cutoff=0.01) #quick
# note, no identity errors

# save intermediate objects
saveRDS(aln,"G19ITS.alignment.RDS")
saveRDS(dm,"G19ITS.dm.RDS")
saveRDS(G19ITS.clusters,"G19ITS.clusters.RDS")

# how to handle the non-finite values of the distance matrix before UPGMA fails?
G19ITS.dm.dataframe<-as.data.frame(dm)
is.na(G19ITS.dm.dataframe)<-sapply(G19ITS.dm.dataframe, is.infinite)
G19ITS.dm.dataframe[is.na(G19ITS.dm.dataframe)]<-1
G19ITS.dm2<-as.matrix(G19ITS.dm.dataframe)
saveRDS(G19ITS.dm2,"G19ITS.ITS.dm2.RDS")

#how do we get that into a tree?
#install.packages("phangorn")
library(phangorn)
G19ITS.upgma<-upgma(G19ITS.dm2)
plot(G19ITS.upgma,type="fan",cex=0.8,show.tip.label=FALSE)
# cool!
saveRDS(G19ITS.upgma,"G19ITS.upgma.tree.RDS")

#add tree to phyloseq object
G19ITS.tree.ps<-Gardner2019.ITS.filter.ps
phy_tree(G19ITS.tree.ps)<-G19ITS.upgma
G19ITS.tree.ps
saveRDS(G19ITS.tree.ps,"G19ITS.tree.ps.RDS")

# based on this tree, it looks like we need a 99% (0.005) glom

# we also need an R function to automatically check this issue!

# glom
G19ITS.glom.005<-tree_glom(G19ITS.tree.ps,resolution=0.005)
G19ITS.glom.005
plot(G19ITS.glom.005@phy_tree)

G19ITS.glom.005@tax_table #ITS are not classified by dada2

writeXStringSet(G19ITS.glom.005@refseq, "G19ITS.glom.005.fasta", format="fasta")

G19ITS.glom.005
taxa_names(G19ITS.glom.005)
G19ITS.2merge <- G19ITS.glom.005
mytaxaNames <- as.character (G19ITS.2merge@refseq)
mytaxaNames
taxa_names(G19ITS.2merg) <- mytaxaNames
head(taxa_names(G19ITS.2merge))


# converts the ref_seq back into the taxa names (from ASVs)
G19ITS.glom.005
taxa_names(G19ITS.glom.005)
G19ITS.2merge<-G19ITS.glom.005
mytaxaNames<-as.character(G19ITS.2merge@refseq)
mytaxaNames
mytaxaNames<-unname(mytaxaNames)
taxa_names(G19ITS.2merge)<-mytaxaNames
head(taxa_names(G19ITS.2merge))

# break the previous object and recombine into a new one
G19ITS.2merge
OTU<-otu_table(G19ITS.2merge,taxa_are_rows=FALSE)
SAMPLEDATA<-sample_data(G19ITS.2merge)
TAX<-tax_table(G19ITS.2merge)
G19ITS.2merge.new<-phyloseq(OTU, SAMPLEDATA, TAX)

# only the columns that are the same name are going to merge, so need to check
colnames(G19ITS.2merge.new@sam_data)

#making heatmap
Gardner.comm<-Gardner.otu
Gardner.relAbundComm <- decostand(Gardner.comm, method = "total") #transforms to relative abundance

print(apply(Gardner.relAbundComm, 1, sum)) #yay
Gardner.relAbundComm[,1:5]
head(colSums(Gardner.relAbundComm)) #note, need to resort

Gardner.sortedRelAbundComm<-Gardner.relAbundComm[,order(-colSums(relAbundComm))]
Gardner.top50comm<-Gardner.sortedRelAbundComm[,1:50]
dim(Gardner.sortedRelAbundComm)
dim(Gardner.top50comm)
Gardner.top50comm[1:5,1:5]

head(colSums(Gardner.top50comm))

class(Gardner.top50comm)
Gardner.top50comm<-as.data.frame(Gardner.top50comm)
mydata<-Gardner.top50comm
mydata$names<-rownames(Gardner.top50comm)
head(mydata)

comm.melted<-melt(mydata)
class(comm.melted)
head(comm.melted)

# rough view #use reshape2 package
ggplot(comm.melted, aes(x=names, y=variable, fill=value)) + geom_tile()

tax_table(Gardner.rel.ps)[10:20]
comm<-Gardner.otu
metadata<-Gardner2019.metadata
colnames(metadata)
# compare species Simpson diversity between host species
invsimpson.comm <- diversity(comm, index = "invsimpson", MARGIN= 1)
boxplot(invsimpson.comm ~ metadata$Host, xlab = "Host", ylab = "InvSimpson Index")
