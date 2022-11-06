getwd()
setwd("/Users/jchung800/Desktop/R/Chen ITS2")

ChenITS2020.seqtab <- read.table("ChenITS2.seqtab",header=TRUE,row.names=1,sep="\t")
ChenITS2020.taxonomy <- read.table ("ChenITS2.taxonomy",header=TRUE,row.names=1,sep="\t")

colnames(ChenITS2020.taxonomy)
colnames(ChenITS2020.seqtab) #SRR numbers

ChenITS2020.metadata<-read.csv("ChenITS2.metadata.csv",header=TRUE,row.names=1)
rownames(ChenITS2020.metadata) #SRR numbers
length(colnames(ChenITS2020.seqtab)) 
length(rownames(ChenITS2020.metadata)) 

#make phyloseq
OTU <- otu_table(ChenITS2020.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(ChenITS2020.metadata)
TAX <- tax_table(as.matrix(ChenITS2020.taxonomy))
ChenITS2020.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
ChenITS2020.ps 


hist(taxa_sums(ChenITS2020.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(ChenITS2020.ps)
Chen.filter.ps<-prune_samples(sample_sums(ChenITS2020.ps) > 19, ChenITS2020.ps) # prune with low sample number
Chen.filter.ps

Chen.filter.ps<-prune_taxa(taxa_sums(Chen.filter.ps) > 19, Chen.filter.ps)
Chen.filter.ps

saveRDS(Chen.filter.ps, "Chen.filter.ps")
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
Chen.dna <- Biostrings::DNAStringSet(taxa_names(Chen.filter.ps))
names(Chen.dna) <- taxa_names(Chen.filter.ps)
ChenITS2020.filter.ps <- merge_phyloseq(Chen.filter.ps, Chen.dna)
taxa_names(ChenITS2020.filter.ps) <- paste0("ChenITS2020.ASV", seq(ntaxa(ChenITS2020.filter.ps)))
ChenITS2020.filter.ps
rownames (ChenITS2020.filter.ps@sam_data)
rownames (ChenITS2020.filter.ps@tax_table)

#transform to relative abundance
Chen.rel.ps = transform_sample_counts(ChenITS2020.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Chen.otu<-psotu2veg(Chen.rel.ps)
ChenITS2020.metadata<-pssd2veg(Chen.rel.ps)
#species glom 

#Bray-Curtis Dissimilarity
library(vegan)
Chen.comm.bc.dist <- vegdist(Chen.otu, method = "bray", na.rm = FALSE)
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(Chen.comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")
saveRDS(Chen.comm.bc.dist, "Chen Bacterial Bray- Curtis")

adonis(Chen.comm.bc.dist ~ tmp, data = ChenITS2020.metadata, permutations = 999)

adonis(Chen.comm.bc.dist ~ HOST, data = ChenITS2020.metadata, permutations = 999)


Chen.simperOutput<-simper(Chen.otu,ChenITS2020.metadata$HOST,perm=99)
Chen.simperOutput

Chen.ruegeria.ps= subset_taxa(Chen.rel.ps, genus=="Ruegeria") 
Chen.ruegeria.ps #7 taxa
title = "plot_bar; Chen Ruegeria Abundance"
plot_bar(Chen.ruegeria.ps, "outcome", "Abundance", title=title) 

Chen.vibrio.ps= subset_taxa(Chen.rel.ps, family=="Vibrionaceae") 
Chen.vibrio.ps #63 taxa
title = "plot_bar; Chen Vibrionaceae Abundance"
plot_bar(Chen.vibrio.ps, "outcome", "Abundance", title=title) 

test= "Diseased_C24e"
gsub("_.+","", test)
gsub(".+_","",test)
gsub("[a-z]","","C24e")

test="Diseased_C24e"
gsub("_.+","",test)
ChenITS2020.metadata$treatment<-gsub("_.+","",ChenITS2020.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(ChenITS2020.metadata)
ChenITS2020.metadata$genotype<-gsub(".+_","",ChenITS2020.metadata$outcome)
ChenITS2020.metadata$genotype<-gsub("[a-z]","",ChenITS2020.metadata$genotype)

# need to pool some of the genotypes into a new category
View(ChenITS2020.metadata)

adonis(Chen.comm.bc.dist ~ genotype, data = ChenITS2020.metadata, permutations = 999)


test="Diseased_C24e"
gsub("_.+","",test)
ChenITS2020.metadata$treatment<-gsub("_.+","",ChenITS2020.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(ChenITS2020.metadata)
ChenITS2020.metadata$genotype<-gsub(".+_","",ChenITS2020.metadata$outcome)
ChenITS2020.metadata$genotype<-gsub("[a-z]","",ChenITS2020.metadata$genotype)

adonis(Chen.comm.bc.dist ~ Host+ genotype + treatment, data = ChenITS2020.metadata, permutations = 999)

# need to pool some of the genotypes into a new category for comparison at end of results

View(ChenITS2020.taxonomy)

Chen.vibrio.ps= subset_taxa(Chen.rel.ps, family=="Fokiniaceae") 
Chen.vibrio.ps #24 taxa
title = "plot_bar; Chen Fokiniaceae Abundance"
plot_bar(Chen.vibrio.ps, "outcome", "Abundance", title=title) 

Chen.vibrio.ps@tax_table
writeXStringSet(Chen.vibrio.ps@refseq, "ChenRickett.fasta", format="fasta")
