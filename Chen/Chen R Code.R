getwd()
setwd("/Users/jchung800/Desktop/R/Chen")

Chen2020.seqtab <- read.table("Chen.16s.seqtab",header=TRUE,row.names=1,sep="\t")
Chen2020.taxonomy <- read.table ("Chen.16s.taxonomy",header=TRUE,row.names=1,sep="\t")
head(rownames(Chen2020.seqtab))
colnames(Chen2020.taxonomy)
colnames(Chen2020.seqtab) #SRR numbers

Chen2020.metadata<-read.csv("Chen.16s.metadata.csv",header=TRUE,row.names=1)
rownames(Chen2020.metadata) #SRR numbers
length(colnames(Chen2020.seqtab)) #95
length(rownames(Chen2020.metadata)) #95

#make phyloseq
OTU <- otu_table(Chen2020.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Chen2020.metadata)
TAX <- tax_table(as.matrix(Chen2020.taxonomy))
Chen2020.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Chen2020.ps 


hist(taxa_sums(Chen.filter.ps),breaks=10000,xlim=c(0,200))

# filter taxa out that do not occur in at least 5 samples:
sample_sums(Chen2020.ps)
Chen.filter.ps<-prune_samples(sample_sums(Chen2020.ps) > 19, Chen2020.ps) # prune with low sample number
Chen.filter.ps
Chen.filter.ps<-prune_taxa(taxa_sums(Chen.filter.ps) > 19, Chen.filter.ps)
Chen.filter.psChen.filter.ps<-filter_taxa(Chen.filter.ps, function(x) sum(x>0) > 4, TRUE)
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
Chen2020.filter.ps <- merge_phyloseq(Chen.filter.ps, Chen.dna)
taxa_names(Chen2020.filter.ps) <- paste0("Chen2020.ASV", seq(ntaxa(Chen2020.filter.ps)))
Chen2020.filter.ps
rownames (Chen2020.filter.ps@sam_data)
rownames (Chen2020.filter.ps@tax_table)

#transform to relative abundance
Chen.rel.ps = transform_sample_counts(Chen2020.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Chen.otu<-psotu2veg(Chen.rel.ps)
Chen2020.metadata<-pssd2veg(Chen.rel.ps)
#species glom 

#Bray-Curtis Dissimilarity
library(vegan)
Chen.comm.bc.dist <- vegdist(Chen.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(Chen.comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")
saveRDS(Chen.comm.bc.dist, "Chen Bacterial Bray- Curtis")

tail(min(sample_sums(Chen.filter.ps)))
genusChen <-tax_glom(Chen.filter.ps, taxrank="genus", NArm=FALSE)

adonis(Chen.comm.bc.dist ~ tmp, data = Chen2020.metadata, permutations = 999)
            #Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#tmp         1      0.46 0.45971 0.93118 0.00427   0.75
#Residuals 217    107.13 0.49369         0.99573       
#Total     218    107.59                 1.00000        

adonis(Chen.comm.bc.dist ~ HOST, data = Chen2020.metadata, permutations = 999)
            #Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#HOST        1     0.692 0.69240  1.4056 0.00644  0.007 **
#Residuals 217   106.898 0.49262         0.99356          
#Total     218   107.590                 1.00000      
#.03 when permutations= 99
#HOST + tmp-> HOST r^2= .006, tmp r^2= .668

Chen.simperOutput<-simper(Chen.otu,Chen2020.metadata$HOST,perm=99)
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
Chen2020.metadata$treatment<-gsub("_.+","",Chen2020.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(Chen2020.metadata)
Chen2020.metadata$genotype<-gsub(".+_","",Chen2020.metadata$outcome)
Chen2020.metadata$genotype<-gsub("[a-z]","",Chen2020.metadata$genotype)

# need to pool some of the genotypes into a new category
View(Chen2020.metadata)

adonis(Chen.comm.bc.dist ~ genotype, data = Chen2020.metadata, permutations = 999)


test="Diseased_C24e"
gsub("_.+","",test)
Chen2020.metadata$treatment<-gsub("_.+","",Chen2020.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(Chen2020.metadata)
Chen2020.metadata$genotype<-gsub(".+_","",Chen2020.metadata$outcome)
Chen2020.metadata$genotype<-gsub("[a-z]","",Chen2020.metadata$genotype)

adonis(Chen.comm.bc.dist ~ Host+ genotype + treatment, data = Chen2020.metadata, permutations = 999)

# need to pool some of the genotypes into a new category for comparison at end of results

View(Chen2020.taxonomy)

Chen.vibrio.ps= subset_taxa(Chen.rel.ps, family=="Fokiniaceae") 
Chen.vibrio.ps #24 taxa
title = "plot_bar; Chen Fokiniaceae Abundance"
plot_bar(Chen.vibrio.ps, "outcome", "Abundance", title=title) 

Chen.vibrio.ps@tax_table
writeXStringSet(Chen.vibrio.ps@refseq, "ChenRickett.fasta", format="fasta")
