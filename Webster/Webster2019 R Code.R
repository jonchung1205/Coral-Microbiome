#8/31/21

getwd()
setwd("/Users/jchung800/Desktop/R/Webster")

Webster2019.seqtab <- read.table ("Webster2019.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Webster2019.taxonomy <- read.table("Webster2019.taxonomy.tabular",header=TRUE,row.names=1,sep="\t" )
Webster2019.metadata<-read.csv("Webster2019.metadata.csv",header=TRUE,row.names=1)

colnames(Webster2019.taxonomy)
colnames(Webster2019.seqtab) #SRR numbers

rownames(Webster2019.metadata) #SRR numbers
length(colnames(Webster2019.seqtab)) #109
length(rownames(Webster2019.metadata)) #109

#make phyloseq
OTU <- otu_table(Webster2019.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Webster2019.metadata)
TAX <- tax_table(as.matrix(Webster2019.taxonomy))
Webster2019.ps<-phyloseq(OTU, SAMPLEDATA, TAX)
Webster2019.ps
saveRDS(Webster2019.ps, "Webster2019.ps.RDS")
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 2538 taxa and 109 samples ]:
#   sample_data() Sample Data:        [ 109 samples by 37 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 2538 taxa by 7 taxonomic ranks ]:

hist(taxa_sums(Webster2019.ps),breaks=10000,xlim=c(0,200))
# filter taxa out that do not occur in at least 5 samples:
sample_sums(Webster2019.ps)
Webster.filter.ps<-prune_samples(sample_sums(Webster2019.ps) > 999, Webster2019.ps) # prune with low sample number
Webster.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 2538 taxa and 98 samples ]:   trimmed down to 98 samples
#   sample_data() Sample Data:        [ 98 samples by 37 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 2538 taxa by 7 taxonomic ranks ]:
Webster.filter.ps<-prune_taxa(taxa_sums(Webster.filter.ps) > 49, Webster.filter.ps)
Webster.filter.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 1335 taxa and 98 samples ]: trimmed down to 1335 taxa
#   sample_data() Sample Data:        [ 98 samples by 37 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 1335 taxa by 7 taxonomic ranks ]:
#   taxa are rows

#removing the chloroplast
Webster2019.noChloroplast.ps<-subset_taxa(Webster2019.ps,order != "Chloroplast")
Webster2019.noChloroplast.ps

head(Webster2019.noChloroplast.ps@tax_table)
myname<-grep("Mitochondria",Webster2019.noChloroplast.ps@tax_table,value=TRUE)
myname

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
Webster.dna <- Biostrings::DNAStringSet(taxa_names(Webster.filter.ps))
names(Webster.dna) <- taxa_names(Webster.filter.ps)
Webster2019.filter.ps <- merge_phyloseq(Webster.filter.ps, Webster.dna)
taxa_names(Webster.filter.ps) <- paste0("Webster2019.ASV", seq(ntaxa(Webster2019.filter.ps)))
Webster2019.filter.ps
rownames (Webster2019.filter.ps@sam_data)
rownames (Webster2019.filter.ps@tax_table)

#transform to relative abundance
Webster.rel.ps = transform_sample_counts(Webster2019.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Webster.otu<-psotu2veg(Webster.rel.ps)
Webster2019.metadata<-pssd2veg(Webster.rel.ps)
Webster.rel.ps_genus <- tax_glom(Webster.rel.ps,taxrank = "class")
Webster.top20 <- names(sort(taxa_sums(Webster.rel.ps_genus),TRUE)[1:20])
plot_heatmap(Webster.top20)
#Bray-Curtis Dissimilarity
library(vegan)
comm.bc.dist <- vegdist(Webster.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")

adonis(comm.bc.dist ~ treatment, data = Webster2019.metadata, permutations = 99)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  3     2.033 0.67783  1.6887 0.05114   0.01 **
# Residuals 94    37.731 0.40139         0.94886          
# Total     97    39.764                 1.00000          
# ---

Webster.simperOutput<-simper(Webster.otu,Webster2019.metadata$treatment,perm=99)
Webster.simperOutput  #ASVs of interest- ASV1- family Rhodobacteraceae, genus Ruegeria, ASV2- Candidatus Crytpoplasma califor
Webster.rel.ps@tax_table

#subset Alphaproteobacteria class
Webster.alpha.ps= subset_taxa(Webster.rel.ps, class=="Alphaproteobacteria") 
Webster.alpha.ps #968 taxa (72.5% of the taxa)
title = "plot_bar; Webster Alphaproteobacteria Abundance"
plot_bar(Webster.alpha.ps, "treatment", "Abundance", title=title) #highest abundance in stress_2100 (abundance approx 28) 

#subset Gammaproteobacteria class
Webster.gamma.ps= subset_taxa(Webster.rel.ps, class=="Gammaproteobacteria") 
Webster.gamma.ps #17 taxa
title = "plot_bar; Webster Gammaproteobacteria Abundance"
plot_bar(Webster.gamma.ps, "treatment", "Abundance", title=title) #highest abundance in stress_ambient (0.19 abundance)

#subset Planctomycetes class
Webster.planc.ps= subset_taxa(Webster.rel.ps, class=="Planctomycetes") 
Webster.planc.ps #94 taxa
title = "plot_bar; Webster Planctomycetes Abundance"
plot_bar(Webster.planc.ps, "treatment", "Abundance", title=title) #highest abundance in control(0.65) and slightly less in stress_2100 (0.64)

#subset endoz family
Webster.endoz.ps= subset_taxa(Webster.rel.ps, family=="Endozoicomonadaceae") 
Webster.endoz.ps #9 taxa
title = "plot_bar; Webster Endozoicomonadaceae Abundance"
plot_bar(Webster.endoz.ps, "treatment", "Abundance", title=title) #highest abundance in stress_2100 (approx 0.11)

#9/1/21
#subset rickettsiales order
Webster.ricket.ps= subset_taxa(Webster.rel.ps, order=="Rickettsiales") 
Webster.ricket.ps #92 taxa
title = "plot_bar; Webster Rickettsiales Abundance"
plot_bar(Webster.ricket.ps, "treatment", "Abundance", title=title) #highest abundance in stress_2100 (approx 5.8)

#subset Rhizobiales order
Webster.rhizo.ps= subset_taxa(Webster.rel.ps, order=="Rhizobiales") 
Webster.rhizo.ps #293 taxa
title = "plot_bar; Webster Rhizobiales Abundance"
plot_bar(Webster.rhizo.ps, "treatment", "Abundance", title=title) #highest abundance in stress_2100 (approx 8.8)

#finding species with greatest contributions
#control_field_control
smallASVs<-Webster.simperOutput$`control_field_control`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV4, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV10, xlab= "Treatment", ylab= "Relative Abundance")

#control_stress_ambient
smallASVs<-Webster.simperOutput$`control_field_control`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV5, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV14, xlab= "Treatment", ylab= "Relative Abundance")

#control_stress_2100
smallASVs<-Webster.simperOutput$`control_stress_2100`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV5, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV14, xlab= "Treatment", ylab= "Relative Abundance")
#field_control_stress_ambient
smallASVs<-Webster.simperOutput$`control_stress_ambient`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV5, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV10, xlab= "Treatment", ylab= "Relative Abundance")
#field_control_stress_2100
smallASVs<-Webster.simperOutput$`control_stress_ambient`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV5, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV10, xlab= "Treatment", ylab= "Relative Abundance")
#stress_ambient_stress_2100
smallASVs<-Webster.simperOutput$`control_stress_ambient`$species[1:5]
smallASVs <- Webster.simperOutput$control_field_control$species
small.otu<-Webster.otu[,smallASVs]
small.otu.means<-aggregate(small.otu,by=list(Webster2019.metadata$treatment),FUN=mean)
small.otu.means
barplot(small.otu.means$Webster2019.ASV2, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV1, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV5, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV3, xlab= "Treatment", ylab= "Relative Abundance")
barplot(small.otu.means$Webster2019.ASV19, xlab= "Treatment", ylab= "Relative Abundance")

#rarefaction
filter.otu <- psotu2veg(Webster2019.filter.ps)
rarefy(filter.otu, 1000)
rarecurve(filter.otu, 1000, xlab="sample size", ylab= "species", label=TRUE )


phyloseq::plot_bar(Webster.rel.ps, fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H
#diversity tests
Webster.comm<-Webster.otu
Webster.metadata<-Webster2019.metadata
colnames(Webster.metadata)
# compare species Simpson diversity between host species
Webster.invsimpson.comm <- diversity(Webster.comm, index = "invsimpson")
Webster.shannon <- diversity(Webster.comm, index = "shannon")
boxplot(Webster.invsimpson.comm ~ Webster.metadata$treatment, xlab = "Treatment", ylab = "InvSimpson Index")
boxplot(Webster.shannon ~ Webster.metadata$treatment, xlab = "Treatment", ylab = "Shannon Index")


#making heatmap
Webster.comm<-Webster.otu
Webster.relAbundComm <- decostand(Webster.comm, method = "total") #transforms to relative abundance

print(apply(Webster.relAbundComm, 1, sum)) #yay
Webster.relAbundComm[,1:5]
head(colSums(Webster.relAbundComm)) #note, need to resort

Webster.sortedRelAbundComm<-Webster.relAbundComm[,order(-colSums(Webster.relAbundComm))]
Webster.top50comm<-Webster.sortedRelAbundComm[,1:50]
dim(Webster.sortedRelAbundComm)
dim(Webster.top50comm)
Webster.top50comm[1:5,1:5]

head(colSums(Webster.top50comm))

class(Webster.top50comm)
Webster.top50comm<-as.data.frame(Webster.top50comm)
mydata<-Webster.top50comm
mydata$names<-rownames(Webster.top50comm)
head(mydata)

comm.melted<-melt(mydata)
class(comm.melted)
head(comm.melted)

# rough view #use reshape2 package
ggplot(comm.melted, aes(x=names, y=variable, fill=value)) + geom_tile()

tax_table(Webster.rel.ps)[10:20]

#ordination
library(MASS)
Webster.dis <- Webster.comm.bc.dist
Webster.mds <- isoMDS(Webster.dis)
stressplot(Webster.mds, Webster.dis)
ordiplot (Webster.mds, type = "t")
#metaMDS
Webster2.mds <- metaMDS (Webster.otu, trace= FALSE)
plot (Webster2.mds, type = "t")
mds_data <- as.data.frame(Webster2.mds$points)
mds_data2 <- rownames(mds_data)
mds_data3 <- dplyr::left_join(mds_data, invsimpson.comm)

biplot(comm.bc.dist,Y=comm.bc.clust, plot.axes=c(1,2), dir.axis1=1, dir.axis2=1)
mds.stuff <- cmdscale(comm.bc.dist,eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),X=mds.values[,1],Y=mds.values[,2])
mds.data
#making a graph 
ggplot(data=mds.data,aes(x=X, y=Y, label=Sample))+
  geom_text()+ 
  theme_bw()+ 
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) + 
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) + 
  ggtitle("MDS plot using Euclidean distance")
  