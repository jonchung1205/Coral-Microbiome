#9/11/21
getwd()
setwd("/home/jc857653/R/Documents/Gardner ITS2")

Gardner2019.ITS.seqtab <- read.table ("Gardner.ITS.seqtab.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Gardner2019.ITS.taxonomy <- read.table("Gardner.ITS.taxonomy.tabular",header=TRUE,row.names=1,sep="\t" )
Gardner2019.ITS.metadata<-read.csv("GardnerITS.metadata.csv",header=TRUE,row.names=1)

colnames(Gardner2019.ITS.taxonomy)
colnames(Gardner2019.ITS.seqtab) 

#make phyloseq
OTU <- otu_table(Gardner2019.ITS.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Gardner2019.ITS.metadata)
TAX <- tax_table(as.matrix(Gardner2019.ITS.taxonomy))
Gardner2019.ITS.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Gardner2019.ITS.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 360 taxa and 48 samples ]:
#   sample_data() Sample Data:        [ 48 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 360 taxa by 7 taxonomic ranks ]:
#   taxa are rows

hist(taxa_sums(Gardner2019.ITS.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Gardner2019.ITS.ps)
Gardner2019.ITS.filter.ps<-prune_samples(sample_sums(Gardner2019.ITS.ps) > 999, Gardner2019.ITS.ps) # prune with low sample number
Gardner2019.ITS.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 360 taxa and 46 samples ]:
#   sample_data() Sample Data:        [ 46 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 360 taxa by 7 taxonomic ranks ]:
Gardner2019.ITS.filter.ps<-prune_taxa(taxa_sums(Gardner2019.ITS.filter.ps) > 49, Gardner2019.ITS.filter.ps)
Gardner2019.ITS.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 264 taxa and 46 samples ]:
#   sample_data() Sample Data:        [ 46 samples by 33 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 264 taxa by 7 taxonomic ranks ]:
#   taxa are rows

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
Gardner2019.ITS.dna <- Biostrings::DNAStringSet(taxa_names(Gardner2019.ITS.filter.ps))
names(Gardner2019.ITS.dna) <- taxa_names(Gardner2019.ITS.filter.ps)
Gardner2019.ITS.filter.ps <- merge_phyloseq(Gardner2019.ITS.filter.ps, Gardner2019.ITS.dna)
taxa_names(Gardner2019.ITS.filter.ps) <- paste0("Gardner2019.ITS.ASV", seq(ntaxa(Gardner2019.ITS.filter.ps)))
Gardner2019.ITS.filter.ps
rownames (Gardner2019.ITS.filter.ps@sam_data)
rownames (Gardner2019.ITS.filter.ps@tax_table)

#transform to relative abundance
Gardner2019.ITS.rel.ps = transform_sample_counts(Gardner2019.ITS.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Gardner2019.ITS.otu<-psotu2veg(Gardner2019.ITS.rel.ps)
Gardner2019.ITS.metadata<-pssd2veg(Gardner2019.ITS.rel.ps)

#Bray-Curtis Dissimilarity
library(vegan)
comm.bc.dist <- vegdist(Gardner2019.ITS.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")
