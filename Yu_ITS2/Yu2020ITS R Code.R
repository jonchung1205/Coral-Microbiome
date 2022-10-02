#9/11/21
getwd()
setwd("C:/Users/jonch/OneDrive/Desktop/R/Yu_ITS2")

library (phyloseq)
library (vegan) 

Yu2020ITS.seqtab <- read.table ("Yu.ITS.seqtab.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Yu2020ITS.taxonomy <- read.table("Yu.ITS.taxonomy.tabular",header=TRUE,row.names=1,sep="\t" )
Yu2020ITS.metadata<-read.csv("YuITS2020.metadata.csv",header=TRUE,row.names=1)

colnames(Yu2020ITS.taxonomy)
colnames(Yu2020ITS.seqtab) 

#make phyloseq
OTU <- otu_table(Yu2020ITS.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Yu2020ITS.metadata)
TAX <- tax_table(as.matrix(Yu2020ITS.taxonomy))
Yu2020ITS.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Yu2020ITS.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1145 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 1145 taxa by 7 taxonomic ranks ]

hist(taxa_sums(Yu2020ITS.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Yu2020ITS.ps)
Yu2020ITS.filter.ps<-prune_samples(sample_sums(Yu2020ITS.ps) > 999, Yu2020ITS.ps) # prune with low sample number
Yu2020ITS.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1145 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 1145 taxa by 7 taxonomic ranks ]   #is it supposed to be the same?
Yu2020ITS.filter.ps<-prune_taxa(taxa_sums(Yu2020ITS.filter.ps) > 49, Yu2020ITS.filter.ps)
Yu2020ITS.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1001 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 1001 taxa by 7 taxonomic ranks ]

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
Yu2020ITS.dna <- Biostrings::DNAStringSet(taxa_names(Yu2020ITS.filter.ps))
names(Yu2020ITS.dna) <- taxa_names(Yu2020ITS.filter.ps)
Yu2020ITS.filter.ps <- merge_phyloseq(Yu2020ITS.filter.ps, Yu2020ITS.dna)
taxa_names(Yu2020ITS.filter.ps) <- paste0("Yu2020ITS.ASV", seq(ntaxa(Yu2020ITS.filter.ps)))
Yu2020ITS.filter.ps
rownames (Yu2020ITS.filter.ps@sam_data)
rownames (Yu2020ITS.filter.ps@tax_table)

#transform to relative abundance
Yu2020ITS.rel.ps = transform_sample_counts(Yu2020ITS.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Yu2020ITS.otu<-psotu2veg(Yu2020ITS.rel.ps)
Yu2020ITS.metadata<-pssd2veg(Yu2020ITS.rel.ps)

#Bray-Curtis Dissimilarity
library(vegan)
comm.bc.dist <- vegdist(Yu2020ITS.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")

adonis(comm.bc.dist ~ Treatment, data = Yu2020ITS.metadata, permutations = 99)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Treatment  4         2     0.5       1 0.28571      1
# Residuals 10         5     0.5         0.71429       
# Total     14         7                 1.00000      
Yu2020ITS.simperOutput<-simper(Yu2020ITS.otu,Yu2020ITS.metadata$Treatment,perm=99)
Yu2020ITS.simperOutput
