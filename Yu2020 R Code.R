#8/16/21
getwd()
setwd("/home/jc857653/R/Documents/Yu") 

#### Yu data files

Yu2020.seqtab<-read.table("Yu2020.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Yu2020.taxonomy<-read.table("Yu2020.taxonomy.tabular",header=TRUE,row.names=1,sep="\t")

colnames(Yu2020.taxonomy)
colnames(Yu2020.seqtab) #SRR numbers

Yu2020.metadata<-read.csv("Yu2020.metadata.correct.csv",header=TRUE,row.names=1)
rownames(Yu2020.metadata) #SRR numbers

length(colnames(Yu2020.seqtab)) #15
length(rownames(Yu2020.metadata)) #15

#sort the metadata into the same sequence of SRR as the seqtab
Yu2020.metadata.sorted<-Yu2020.metadata[colnames(Yu2020.seqtab),]

all.equal(colnames(Yu2020.seqtab),rownames(Yu2020.metadata.sorted)) #TRUE

#make phyloseq
OTU <- otu_table(Yu2020.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Yu2020.metadata.sorted)
TAX <- tax_table(as.matrix(Yu2020.taxonomy))
Yu2020.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Yu2020.ps

Yu2020.noChloroplast.ps<-subset_taxa(Yu2020.ps,order != "Chloroplast")
Yu2020.noChloroplast.ps

head(Yu2020.noChloroplast.ps@tax_table)
myname<-grep("Mitochondria",Yu2020.noChloroplast.ps@tax_table,value=TRUE)
myname

Yu2020.noChloroMit.ps<-subset_taxa(Yu2020.noChloroplast.ps,family != "Mitochondria")
Yu2020.noChloroMit.ps

hist(taxa_sums(Yu2020.ps),breaks=10000,xlim=c(20,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Yu2020.noChloroMit.ps)
Yu.filter.ps<-prune_samples(sample_sums(Yu2020.noChloroMit.ps) > 999, Yu2020.noChloroMit.ps) # prune with low sample number
Yu.filter.ps
Yu.filter.ps<-prune_taxa(taxa_sums(Yu.filter.ps) > 49, Yu.filter.ps)
Yu.filter.ps
min(sample_sums(Yu.filter.ps)) #2724
50/2724 #1.8%

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

library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(Yu.filter.ps))
names(dna) <- taxa_names(Yu.filter.ps)
Yu2020.filter.ps <- merge_phyloseq(Yu.filter.ps, dna)
taxa_names(Yu2020.filter.ps) <- paste0("Yu2020.ASV", seq(ntaxa(Yu2020.filter.ps)))
Yu2020.filter.ps
head(taxa_names(Yu2020.filter.ps))

#transform to relative abundance
Yu.rel.ps = transform_sample_counts(Yu2020.filter.ps, function(x) x/sum(x)) #changed to Yu2020.ps without filter 
#extract vegan format otu table and metadata
Yu.otu<-psotu2veg(Yu.rel.ps)
Yu.metadata<-pssd2veg(Yu.rel.ps)

#Bray-Curtis Dissimilarity
library(vegan)
comm.bc.dist <- vegdist(Yu.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis dissimilarity")
#look at most abundant rows
Yu2020.sorted.otu<-Yu.otu[,order(-colSums(Yu.otu))]
Yu2020.sorted.otu[1:15,1:10]
mycolumns<-colnames(Yu2020.sorted.otu)[1:20]
mycolumns
Yu.rel.ps@tax_table[mycolumns]


adonis(comm.bc.dist ~ , data = Yu.metadata, permutations = 99)
simperOutput<-simper(Yu.otu,Yu.metadata$Treatment,perm=99)

#ordination 
Yu.mds <- isoMDS(comm.bc.dist)
stressplot(Yu.mds, comm.bc.dist)
ordiplot (Yu.mds, type = "t")
#metaMDS
Yu2.mds <- metaMDS (Yu.otu, trace= FALSE)
plot (Yu.mds, type = "t")
#rank index 
rankindex (scale(Yu.otu), Yu.otu, c("euc", "man", "bray", "jac", "kul"))
#euclidean distance
Yu.dis <-  vegdist(decostand(Yu.otu, "norm"), "euclid")
#Hellinger distance
Yu.dis2 <- vegdist(decostand(Yu.otu, "hell"), "euclidean")
#procrustes rotation
Yu.pro <- procrustes (Yu.mds, Yu.mds) 
Yu.pro #proceustes sum of squares 267.7 
plot (Yu.pro) 
plot (Yu.pro, kind=2) #what does kind=2 mean?
#PCA analysis
Yu.pca <- rda (Yu.otu, scale= TRUE)
Yu.pca
plot (Yu.pca, scaling=3)
#correspondence analysis
Yu.ca <- cca(Yu.otu)
Yu.ca
plot (Yu.ca)
chisq.test(Yu.otu/sum(Yu.otu))

BiocManager::install("DECIPHER")
library(DECIPHER)

#8/18/21
# how to remove Chloroplasts
Yu2020.noChloroplast.ps<-subset_taxa(Yu2020.ps,order != "Chloroplast")
Yu2020.noChloroplast.ps

head(Yu2020.noChloroplast.ps@tax_table)
myname<-grep("Mitochondria",Yu2020.noChloroplast.ps@tax_table,value=TRUE)
myname

Yu2020.noChloroMit.ps<-subset_taxa(Yu2020.noChloroplast.ps,family != "Mitochondria")
Yu2020.noChloroMit.ps

Yu.filter.ps<-filter_taxa(Yu2020.noChloroMit.ps, function(x){sum(x > 0) > 5}, prune = TRUE) 
Yu.filter.ps
Yu.filter.ps<-prune_samples(sample_sums(Yu.filter.ps) > 999, Yu.filter.ps) # prune with low sample number
Yu.filter.ps
Yu.filter.ps<-prune_taxa(taxa_sums(Yu.filter.ps) > 49, Yu.filter.ps)

sample_sums(Yu2020.noChloroMit.ps)
Yu.filter.ps<-prune_samples(sample_sums(Yu2020.noChloroMit.ps) > 999, Yu2020.noChloroMit.ps) # prune with low sample number
Yu.filter.ps
Yu.filter.ps<-prune_taxa(taxa_sums(Yu.filter.ps) > 49, Yu.filter.ps)
Yu.filter.ps
min(sample_sums(Yu.filter.ps)) #2724
50/2724 #1.8%

# after filtering, create the ASVs

library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(Yu.filter.ps))
names(dna) <- taxa_names(Yu.filter.ps)
Yu2020.filter.ps <- merge_phyloseq(Yu.filter.ps, dna)
taxa_names(Yu2020.filter.ps) <- paste0("Yu2020.ASV", seq(ntaxa(Yu2020.filter.ps)))
Yu2020.filter.ps
head(taxa_names(Yu2020.filter.ps))

saveRDS(Yu2020.filter.ps,"Yu2020.filter.ps.RDS")



#transform to relative abundance
Yu.rel.ps = transform_sample_counts(Yu2020.filter.ps, function(x) x/sum(x)) #changed to Yu2020.ps without filter 
#extract vegan format otu table and metadata
Yu.otu<-psotu2veg(Yu.rel.ps)
Yu.metadata<-pssd2veg(Yu.rel.ps)

#transform filtered data to relative abundance
sample_sums(Yu.filter.ps)
Yu2020.rel.ps = transform_sample_counts(Yu2020.filter.ps, function(x) x/sum(x))
sample_sums(Yu2020.rel.ps)

library(vegan)
comm.bc.dist <- vegdist(Yu.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis dissimilarity")

#look at most abundant rows
Yu2020.sorted.otu<-Yu2020.otu[,order(-colSums(Yu2020.otu))]
Yu2020.sorted.otu[1:15,1:10]

mycolumns<-colnames(Yu2020.sorted.otu)[1:20]
mycolumns
Yu2020.rel.ps@tax_table[mycolumns]

BiocManager::install("DECIPHER")
library(DECIPHER)

# isolate the sequences
Yu.refseq<-Yu2020.rel.ps@refseq
Yu.refseq
# align
aln<-AlignSeqs(Yu.refseq,normPower=0, processors=10) #<5 minutes
# distance matrix
dm<-DistanceMatrix(aln,includeTerminalGaps = FALSE,processors = 10) #quick
# find clusters
Yu.clusters<-IdClusters(dm,method="UPGMA",type="clusters",cutoff=0.01) #quick
# note, no identity errors

# save intermediate objects
saveRDS(Yu.seq.ps,"Yu.seq.ps.RDS")
saveRDS(aln,"Yu.alignment.RDS")
saveRDS(dm,"Yu.dm.RDS")
saveRDS(Yu.clusters,"Yu.clusters.RDS")

# how to handle the non-finite values of the distance matrix before UPGMA fails?
Yu.dm.dataframe<-as.data.frame(dm)
is.na(Yu.dm.dataframe)<-sapply(Yu.dm.dataframe, is.infinite)
Yu.dm.dataframe[is.na(Yu.dm.dataframe)]<-1
Yu.dm2<-as.matrix(Yu.dm.dataframe)
saveRDS(Yu.dm2,"Yu.dm2.RDS")

#how do we get that into a tree?
install.packages("phangorn")
library(phangorn)
Yu.upgma<-upgma(Yu.dm2)
plot(Yu.upgma,type="fan",cex=0.8,show.tip.label=FALSE)
# cool!
saveRDS(Yu.upgma,"Yu.upgma.tree.RDS")

#add tree to phyloseq object
Yu2020.tree.ps<-Yu2020.rel.ps
phy_tree(Yu2020.tree.ps)<-Yu.upgma
Yu2020.tree.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2789 taxa and 33 samples ]
# sample_data() Sample Data:       [ 33 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 2789 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2789 tips and 2788 internal nodes ]
# refseq()      DNAStringSet:      [ 2789 reference sequences ]
saveRDS(Yu.2020tree.ps,"Yu2020.tree.ps.RDS")

# get speedyseq
install.packages("remotes")
remotes::install_github("mikemc/speedyseq")

library(speedyseq)
Yu2020.glom.005<-tree_glom(Yu2020.tree.ps,resolution=0.005)
Yu2020.glom.005
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 391 taxa and 15 samples ]:
#   sample_data() Sample Data:        [ 15 samples by 31 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 391 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 391 tips and 390 internal nodes ]:
#   refseq()      DNAStringSet:       [ 391 reference sequences ]
# taxa are rows


# yay, now we have an object with 99% OTUs!
colnames(Yu2020.glom.005@otu_table)
head(rownames(Yu2020.glom.005@otu_table))
# it is not worth changing the names, except ASV to OTU at publication stage
saveRDS(Yu2020.glom.005,"Yu2020.glom.005.RDS")


# how to combine into genera with tax_glom
# find out capitalization
colnames(Yu2020.glom.005@tax_table) #genus
Yu2020.genus.glom<-tax_glom(Yu2020.glom.005, taxrank="genus")

Yu2020.otu<-psotu2veg(Yu2020.genus.glom)
Yu2020.metadata<-pssd2veg(Yu2020.genus.glom)
comm.bc.dist <- vegdist(Yu2020.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis dissimilarity")

colnames(Yu2020.metadata)
adonis(comm.bc.dist ~ Treatment, data = Yu2020.metadata, permutations = 999)

simperOutput<-simper(Yu2020.otu,Yu2020.metadata$Treatment,perm=999)
simperOutput

# sample comparison of means
smallASVs<-simperOutput$`bacteria-HS_bacteria-R`$species[1:2]
small.otu<-Yu2020.otu[,smallASVs]
small.otu
small.otu.means<-aggregate(small.otu,by=list(Yu2020.metadata$Treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Yu2020.ASV1)

Yu2020.genus.glom@tax_table[smallASVs]

#9/3/21
comm<-Yu2020.otu
metadata<-Yu2020.metadata
colnames(metadata)
# compare species Simpson diversity between host species
invsimpson.comm <- diversity(comm, index = "invsimpson")
boxplot(invsimpson.comm ~ metadata$Treatment, xlab = "Treatment", ylab = "InvSimpson Index")

plotData<-as.data.frame(invsimpson.comm)
plotData

plotData$treatment<-metadata$Treatment
plotData
metadata$Treatment
rownames(metadata)

colnames(plotData)[1]<-"invSimpson"
plotData
ggplot(plotData, aes(x="treatment", y = "invSimpson")) + geom_boxplot()

install.packages("devtools")
devtools::install_github('aliceyiwang/mvabund')
library(mvabund)

#9/12 creating graphs 
phyloseq::plot_bar(Yu2020.rel.ps, fill="Phylum") +
  geom_bar(aes(color= Phylum, fill= Phylum), stat= "identity", position="stack") + 
  labs= (x= "", y= "Relative Abundance\n") +
  facet_wrap(~Status, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ward <- as.dendrogram(hclust(comm.bc.dist, method = "ward.D2"))
meta <- data.frame (phyloseq::Yu2020.metadata(Yu2020.rel.ps))
colorCode <- c('HS' ="red", 'RAHS' = "blue",'C'= "tomato",'A'= "yellow", 'RA'="springgreen" )
labels_colors(ward) <- colorCode[Yu2020.metadata$Treatment][order.dendrogram(ward)]
plot (ward)
  