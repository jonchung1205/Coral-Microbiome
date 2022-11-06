getwd()
setwd("/Users/jchung800/Desktop/R/Rosales")

Rosales2019.seqtab <- read.table("Rosales.seqtab",header=TRUE,row.names=1,sep="\t")
Rosales2019.taxonomy <- read.table ("Rosales.taxonomy",header=TRUE,row.names=1,sep="\t")

colnames(Rosales2019.taxonomy)
colnames(Rosales2019.seqtab) #SRR numbers

Rosales2019.metadata<-read.csv("Rosales.metadata.csv",header=TRUE,row.names=1)
rownames(Rosales2019.metadata) #SRR numbers
length(colnames(Rosales2019.seqtab)) #95
length(rownames(Rosales2019.metadata)) #95
Rosales2019.metadata$outcome<-gsub("_.+","",Rosales2019.metadata$outcome)

#make phyloseq
OTU <- otu_table(Rosales2019.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Rosales2019.metadata)
TAX <- tax_table(as.matrix(Rosales2019.taxonomy))
Rosales2019.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Rosales2019.ps 
# phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 2783 taxa and 95 samples ]:
#sample_data() Sample Data:        [ 95 samples by 31 sample variables ]:
#tax_table()   Taxonomy Table:     [ 2783 taxa by 7 taxonomic ranks ]:
#   taxa are rows

hist(taxa_sums(Rosales2019.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Rosales2019.ps)
Rosales.filter.ps<-prune_samples(sample_sums(Rosales2019.ps) > 999, Rosales2019.ps) # prune with low sample number
Rosales.filter.ps
# phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 2783 taxa and 94 samples ]:
#sample_data() Sample Data:        [ 94 samples by 31 sample variables ]:
#tax_table()   Taxonomy Table:     [ 2783 taxa by 7 taxonomic ranks ]:
#   taxa are rows
Rosales.filter.ps<-prune_taxa(taxa_sums(Rosales.filter.ps) > 49, Rosales.filter.ps)
Rosales.filter.ps
# phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 1721 taxa and 94 samples ]:
#sample_data() Sample Data:        [ 94 samples by 31 sample variables ]:
#tax_table()   Taxonomy Table:     [ 1721 taxa by 7 taxonomic ranks ]:
#   taxa are rows
saveRDS(Rosales.filter.ps, "Rosales.filter.ps")
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
Rosales.dna <- Biostrings::DNAStringSet(taxa_names(Rosales.filter.ps))
names(Rosales.dna) <- taxa_names(Rosales.filter.ps)
Rosales2019.filter.ps <- merge_phyloseq(Rosales.filter.ps, Rosales.dna)
taxa_names(Rosales2019.filter.ps) <- paste0("Rosales2019.ASV", seq(ntaxa(Rosales2019.filter.ps)))
Rosales2019.filter.ps
rownames (Rosales2019.filter.ps@sam_data)
rownames (Rosales2019.filter.ps@tax_table)

#transform to relative abundance
Rosales.rel.ps = transform_sample_counts(Rosales2019.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Rosales.otu<-psotu2veg(Rosales.rel.ps)
Rosales2019.metadata<-pssd2veg(Rosales.rel.ps)
#species glom 

#Bray-Curtis Dissimilarity
library(vegan)
Rosales.comm.bc.dist <- vegdist(Rosales.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(Rosales.comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")
saveRDS(Rosales.comm.bc.dist, "Rosales Bacterial Bray- Curtis")

adonis(Rosales.comm.bc.dist ~ outcome, data = Rosales2019.metadata, permutations = 99)
#           Df SumsOfSqs MeanSqs F.Model R2 Pr(>F)
#outcome   93    26.099       0       0  1      1
#Residuals  0     0.000     Inf          0       
#Total     93    26.099                  1       

Rosales.simperOutput<-simper(Rosales.otu,Rosales2019.metadata$outcome,perm=99)
Rosales.simperOutput

Rosales.ruegeria.ps= subset_taxa(Rosales.rel.ps, genus=="Ruegeria") 
Rosales.ruegeria.ps #7 taxa
title = "plot_bar; Rosales Ruegeria Abundance"
plot_bar(Rosales.ruegeria.ps, "outcome", "Abundance", title=title) 

Rosales.vibrio.ps= subset_taxa(Rosales.rel.ps, family=="Vibrionaceae") 
Rosales.vibrio.ps #63 taxa
title = "plot_bar; Rosales Vibrionaceae Abundance"
plot_bar(Rosales.vibrio.ps, "outcome", "Abundance", title=title) 

library(viridis)

phyloseq::plot_bar(Rosales.rel.ps, fill = "genus") + 
  geom_bar(aes(color = genus, fill = genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ outcome, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H
test= "Diseased_C24e"
gsub("_.+","", test)
gsub(".+_","",test)
gsub("[a-z]","","C24e")

test="Diseased_C24e"
gsub("_.+","",test)
Rosales2019.metadata$treatment<-gsub("_.+","",Rosales2019.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(Rosales2019.metadata)
Rosales2019.metadata$genotype<-gsub(".+_","",Rosales2019.metadata$outcome)
Rosales2019.metadata$genotype<-gsub("[a-z]","",Rosales2019.metadata$genotype)

# need to pool some of the genotypes into a new category
View(Rosales2019.metadata)

adonis(Rosales.comm.bc.dist ~ genotype, data = Rosales2019.metadata, permutations = 999)


test="Diseased_C24e"
gsub("_.+","",test)
Rosales2019.metadata$treatment<-gsub("_.+","",Rosales2019.metadata$outcome)

gsub(".+_","",test)
gsub("[a-z]","","C24e")

colnames(Rosales2019.metadata)
Rosales2019.metadata$genotype<-gsub(".+_","",Rosales2019.metadata$outcome)
Rosales2019.metadata$genotype<-gsub("[a-z]","",Rosales2019.metadata$genotype)

adonis(Rosales.comm.bc.dist ~ Host+ genotype + treatment, data = Rosales2019.metadata, permutations = 999)

# need to pool some of the genotypes into a new category for comparison at end of results

View(Rosales2019.taxonomy)

Rosales.vibrio.ps= subset_taxa(Rosales.rel.ps, family=="Fokiniaceae") 
Rosales.vibrio.ps #24 taxa
title = "plot_bar; Rosales Fokiniaceae Abundance"
plot_bar(Rosales.vibrio.ps, "outcome", "Abundance", title=title) 

Rosales.vibrio.ps@tax_table
writeXStringSet(Rosales.vibrio.ps@refseq, "RosalesRickett.fasta", format="fasta")


# be careful because the subset_taxa also removes "NA"
Rosales2019.bt<-subset_taxa(Rosales2019.ps, (family!="Mitochondria" | is.na(family)))
Rosales2019.bt # 2701 taxa
Rosales2019.bt<-subset_taxa(Rosales2019.bt, (order!="Chloroplast" | is.na(order)))
Rosales2019.bt # 2481 taxa
colnames(Rosales2019.bt@tax_table)
nokingdom<-subset_taxa(Rosales2019.bt, is.na(domain)) #none exist
bt.taxa<-as.data.frame(Rosales2019.bt@tax_table)
unique(bt.taxa$domain) # Bacteria, Eukaryota, Archaea
Rosales2019.bt<-subset_taxa(Rosales2019.bt, (domain!="Eukaryota"))
Rosales2019.bt # 2474 taxa

#ASVs were removed if they did not appear in at least 4 samples
Rosales2019.bt2<-filter_taxa(Rosales2019.bt, function(x) sum(x>0) > 3, TRUE)
Rosales2019.bt2 #494 taxa!


