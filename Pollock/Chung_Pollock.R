#8/9/21
getwd()
setwd("/Users/jchung800/Desktop/R/Pollock")
#Pollock Dataset
Pollock123.seqtab<-read.table("Pollock.Set123.Galaxy2150.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Pollock123.taxonomy<-read.table("Pollock.Set123.Galaxy2151.taxonomy.tabular",header=TRUE,row.names=1,sep="\t")

Pollock456.seqtab<-read.table("Pollock.Set456.Galaxy3317.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Pollock456.taxonomy<-read.table("Pollock.Set456.Galaxy3319.taxonomy.tabular",header=TRUE,row.names=1,sep="\t")

Pollock789.seqtab<-read.table("Pollock.Set789.Galaxy5254.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
Pollock789.taxonomy<-read.table("Pollock.Set789.Galaxy5255.taxonomy.tabular",header=TRUE,row.names=1,sep="\t")


Pollock.metadata<-read.csv("Pollock.2018.Metadata.SraRunTable.txt",header=TRUE,row.names=1)

#are taxa in rows?
head(Pollock123.seqtab) #yes!

#do we need to subset the metadata? -- yes, I tried without and it did not work
OTU <- otu_table(Pollock123.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Pollock.metadata)
TAX <- tax_table(as.matrix(Pollock123.taxonomy))
Pollock123.ps<-phyloseq(OTU, SAMPLEDATA, TAX) #this fails because some samples not included

Pollock123.ps<-phyloseq(OTU, TAX)

OTU <- otu_table(Pollock456.seqtab,taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(Pollock456.taxonomy))
Pollock456.ps<-phyloseq(OTU, TAX)
OTU <- otu_table(Pollock789.seqtab,taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(Pollock789.taxonomy))
Pollock789.ps<-phyloseq(OTU, TAX)

Pollock2018.ps<-merge_phyloseq(Pollock123.ps,Pollock456.ps,Pollock789.ps)

sample_names(SAMPLEDATA)
sample_names(Pollock2018.ps)

head(sample_names(Pollock2018.ps))

names2edit<-sample_names(Pollock2018.ps)
View(Pollock.metadata) #the column alias is almost what we want
names2edit<-sample_names(Pollock2018.ps)

names2edit<-sample_names(Pollock2018.ps)
head(names2edit)
newNames<-gsub("_\\w+_\\w+_\\w+.fastq","",names2edit)
head(newNames)

sample_names(Pollock2018.ps)<-newNames

rownames(Pollock.metadata)<-Pollock.metadata$alias

SAMPLEDATA <- sample_data(Pollock.metadata)

Pollock2018.ps<-merge_phyloseq(Pollock2018.ps, SAMPLEDATA)
length(sample_names(Pollock2018.ps))
length(sample_names(SAMPLEDATA))
subset.Pollock.metadata<-Pollock.metadata[which(row.names(Pollock.metadata) %in% newNames),]

#restart
Pollock2018.ps<-merge_phyloseq(Pollock123.ps,Pollock456.ps,Pollock789.ps)
Pollock2018.ps

sample_names(SAMPLEDATA)
sample_names(Pollock2018.ps)

colnames(Pollock.metadata)
View(Pollock.metadata)

names2edit<-sample_names(Pollock2018.ps)
head(names2edit)
newNames<-gsub("_\\w+_\\w+_\\w+.fastq","",names2edit)
head(newNames)

sample_names(Pollock2018.ps)<-newNames

head(rownames(Pollock.metadata))
rownames(Pollock.metadata)<-Pollock.metadata$alias
SAMPLEDATA <- sample_data(Pollock.metadata)

Pollock2018b.ps<-merge_phyloseq(Pollock2018.ps, SAMPLEDATA)
Pollock2018b.ps

Pollock.metadata2<-read.table("gcmp16S_map_r25.txt",header=TRUE,row.names=1,sep="\t")
colnames(Pollock.metadata2)
head(rownames(Pollock.metadata2))
head(colnames(Pollock.metadata2))
subset.Pollock.metadata2<-Pollock.metadata2[which(row.names(Pollock.metadata2) %in% newNames),]
View(subset.Pollock.metadata2)
Pollock2018.ps<-merge_phyloseq(Pollock123.ps,Pollock456.ps,Pollock789.ps)

sample_names(Pollock2018.ps)<-newNames

sample_data(Pollock2018.ps)<-subset.Pollock.metadata2
Pollock2018combined2.ps<-Pollock2018.ps
saveRDS(Pollock2018combined2.ps,"Pollock2018combined2.ps")
Pollock2018combined2.ps

# how to remove Chloroplasts
Pollock.noChloroplast.ps<-subset_taxa(Pollock2018.ps,order != "Chloroplast")
Pollock.noChloroplast.ps

head(Pollock.noChloroplast.ps@tax_table)
myname<-grep("Mitochondria",Pollock.noChloroplast.ps@tax_table,value=TRUE)
myname

Pollock.noChloroMit.ps<-subset_taxa(Pollock.noChloroplast.ps,family != "Mitochondria")
Pollock.noChloroMit.ps

Pollock.filter.ps<-filter_taxa(Pollock.noChloroMit.ps, function(x){sum(x > 0) > 5}, prune = TRUE) 
Pollock.filter.ps
Pollock.filter.ps<-prune_samples(sample_sums(Pollock.filter.ps) > 999, Pollock.filter.ps) # prune with low sample number
Pollock.filter.ps
Pollock.filter.ps<-prune_taxa(taxa_sums(Pollock.filter.ps) > 49, Pollock.filter.ps)

sample_sums(Pollock.noChloroMit.ps)
Pollock.filter.ps<-prune_samples(sample_sums(Pollock.noChloroMit.ps) > 999, Pollock2020.noChloroMit.ps) # prune with low sample number
Pollock.filter.ps
Pollock.filter.ps<-prune_taxa(taxa_sums(Pollock.filter.ps) > 49, Pollock.filter.ps)
Pollock.filter.ps
min(sample_sums(Pollock.filter.ps)) #2724

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
dna <- Biostrings::DNAStringSet(taxa_names(Pollock.filter.ps))
names(dna) <- taxa_names(Pollock.filter.ps)
Pollock2020.filter.ps <- merge_phyloseq(Pollock.filter.ps, dna)
taxa_names(Pollock2020.filter.ps) <- paste0("Pollock2020.ASV", seq(ntaxa(Pollock2020.filter.ps)))
Pollock2020.filter.ps
head(taxa_names(Pollock2020.filter.ps))

#transform to relative abundance
Pollock.rel.ps = transform_sample_counts(Pollock2020.filter.ps, function(x) x/sum(x)) #changed to Pollock2020.ps without filter 
#extract vegan format otu table and metadata
Pollock.otu<-psotu2veg(Pollock.rel.ps)
Pollock.metadata<-pssd2veg(Pollock.rel.ps)




#with the mitochondria and chloroplasts
#8/13/21
hist(taxa_sums(Pollock2018combined2.ps),breaks=10000,xlim=c(20,200))
# filter taxa out that do not occur in at least 5 samples:
Pollock2018.filter.ps<-filter_taxa(Pollock2018combined2.ps, function(x){sum(x > 0) > 5}, prune = TRUE)
Pollock2018.filter.ps
Pollock2018.filter.ps<-prune_samples(sample_sums(Pollock2018.filter.ps) > 999, Pollock2018.filter.ps)
Pollock2018.filter.ps
Pollock2018.filter.ps<-prune_taxa(taxa_sums(Pollock2018.filter.ps) > 49, Pollock2018.filter.ps)

plot_bar(Pollock2018.filter.ps, "host_genus_id", "Abundance", "phylum")
plot_bar(Pollock.2018.filter.ps, "phylum", "Abundance", "host_genus_id")


library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(Pollock2018.filter.ps))
names(dna) <- taxa_names(Pollock2018.filter.ps)
Pollock2020.filter.ps <- merge_phyloseq(Pollock2018.filter.ps, dna)
taxa_names(Pollock2018.filter.ps) <- paste0("Pollock2018.ASV", seq(ntaxa(Pollock2018.filter.ps)))
Pollock2018.filter.ps
head(taxa_names(Pollock2018.filter.ps))

#transform to relative abundance
Pollock2018.rel.ps = transform_sample_counts(Pollock2018.filter.ps, function(x) x/sum(x))
#extract vegan format otu table and metadata
Pollock.otu<-psotu2veg(Pollock2018.rel.ps)
Pollock.metadata<-pssd2veg(Pollock2018.rel.ps)

sum(taxa_sums(Pollock2018combined2.ps)>200) #3013
sum(taxa_sums(Pollock2018combined2.ps)>9) #40495
sum(taxa_sums(Pollock2018combined2.ps)>49) #32839

sum(sample_sums(Pollock2018combined2.ps)>999) #421

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
#transform filtered data to relative abundance
sample_sums(Pollock2018.filter.ps)
Pollock2018.rel.ps = transform_sample_counts(Pollock2018.filter.ps, function(x) x/sum(x))
sample_sums(Pollock2018.rel.ps)

#save as RDS
saveRDS(Pollock2018.rel.ps,"Pollock2018.rel.ps.RDS")

Pollock.otu<-psotu2veg(Pollock2018.rel.ps)
Pollock.metadata<-pssd2veg(Pollock2018.rel.ps)

#Bray-Curtis Dissimilarity
library(vegan)
comm.bc.dist <- vegdist(Pollock.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis dissimilarity")

colnames(Pollock.metadata)
# big effects for both host_species and temperature (10% of variation at least)


ord<-metaMDS(Pollock.otu)
ordiplot(ord, display = "sites", type = "points")

vectorData<-Pollock.metadata[,c("temperature","salinity","depth")]

ef<- envfit(ord, vectorData, permu = 99)
ef
ef2<-vectorfit(ord, vectorData, permu=99) #prefered,but need to filter out the rows of unknown or NA

plot(ord, display = "sites")
plot(ef, p.max = 0.1)

ef3<-vectorfit(ord, vectorData$temperature, permu=99, na.rm=TRUE )
View(vectorData)

#removing no data
x <- Pollock.metadata[complete.cases(Pollock.metadata), ]
str(x)
x <- subset.Pollock.metadata2[complete.cases(subset.Pollock.metadata2), ]
str(x)
View (subset.Pollock.metadata2)
subset.Pollock.metadata2[!complete.cases(subset.Pollock.metadata2),]
subset.Pollock.metadata3 <- na.omit(subset.Pollock.metadata2)
View(subset.Pollock.metadata3)

#8/13/21
simperOutput<-simper(Pollock.otu,Pollock.metadata$temperature,perm=99)
adonis(comm.bc.dist ~ host_species + temperature, data = Pollock.metadata, permutations = 999)
adonis (comm.bc.dist~temperature, data= Pollock.metadata,perumations =999 )
#9/13
library(viridis)
phyloseq::plot_bar(Pollock2018.rel.ps, fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ temperature, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H

Pollock.alpha.ps= subset_taxa(Pollock2018.rel.ps, class=="Alphaproteobacteria") 
Pollock.alpha.ps #1779 taxa
title = "plot_bar; Pollock Alphaproteobacteria Abundance"
plot_bar(Pollock.alpha.ps, "temperature", "Abundance", title=title) #highest abundance in Dijbouti (3.6 abundance)

Pollock.gamma.ps= subset_taxa(Pollock2018.rel.ps, class=="Gammaproteobacteria") 
Pollock.gamma.ps #4212 taxa
title = "plot_bar; Pollock Gammaproteobacteria Abundance"
plot_bar(Pollock.gamma.ps, "temperature", "Abundance", title=title) #highest abundance in French Polynesia (37 abundance)

Pollock.endoz.ps= subset_taxa(Pollock2018.rel.ps, family=="Endozoicomonadaceae") 
Pollock.endoz.ps #235 taxa
title = "plot_bar; Pollock Endozoicomonadaceae Abundance"
plot_bar(Pollock.endoz.ps, "temperature", "Abundance", title=title) #highest abundance in French Polynesia (34 abundance)

Pollock.rhizo.ps= subset_taxa(Pollock2018.rel.ps, order=="Rhizobiales") 
Pollock.rhizo.ps #203 taxa
title = "plot_bar; Pollock Rhizobiales Abundance"
plot_bar(Pollock.rhizo.ps, "temperature", "Abundance", title=title) #highest abundance in New Caledonia (0.34 abundance)

Pollock.vibrio.ps= subset_taxa(Pollock2018.rel.ps, family=="Vibrionaceae") 
Pollock.vibrio.ps #134 taxa
title = "plot_bar; Pollock Vibrionaceae Abundance"
plot_bar(Pollock.vibrio.ps, "temperature", "Abundance", title=title) #highest abundance in New Caledonia (0.34 abundance)

adiv <- data.frame(
  "Observed" = specnumber(Pollock.otu),
  "Chao1" = phyloseq::estimate_richness(Pollock.metadata, measures = "Chao1"),
  "InvSimpson" = diversity(Pollock.otu, index="invsimpson"),
  "Treatment" = phyloseq::sample_data(Pollock.metadata)$temperature)
head(adiv)
library(dplyr)
library(tidyr)
adiv %>%
  gather(key = metric, value = value, c("Observed", "InvSimpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "InvSimpson"))) %>%
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1))+ 
  theme(axis.text.y= element_text(size= 12))

#diversity
Pollock.comm<-Pollock.otu
Pollock.metadata<-Pollock.metadata
colnames(Pollock.metadata)
# compare species Simpson diversity between host species
Pollock.invsimpson.comm <- diversity(Pollock.comm, index = "invsimpson")
Pollock.shannon <- diversity(Pollock.comm, index = "shannon")
boxplot(Pollock.invsimpson.comm ~ Pollock.metadata$temperature, xlab = "Temperature", ylab = "InvSimpson Index")
boxplot(Pollock.shannon ~ Pollock.metadata$temperature, xlab = "Temperature", ylab = "Shannon Index")

#making heatmap
Pollock.comm<-Pollock.otu
Pollock.relAbundComm <- decostand(Pollock.comm, method = "total") #transforms to relative abundance

# check relative abundance in each sample (each should now be 1)
print(apply(Pollock.relAbundComm, 1, sum)) #yay
Pollock.relAbundComm[,1:5]
head(colSums(Pollock.relAbundComm)) #note, need to resort

Pollock.sortedRelAbundComm<-Pollock.relAbundComm[,order(-colSums(relAbundComm))]
Pollock.top50comm<-Pollock.sortedRelAbundComm[,1:50]
dim(Pollock.sortedRelAbundComm)
dim(Pollock.top50comm)
Pollock.top50comm[1:5,1:5]

head(colSums(Pollock.top50comm))

class(Pollock.top50comm)
Pollock.top50comm<-as.data.frame(Pollock.top50comm)
mydata<-Pollock.top50comm
mydata$names<-rownames(Pollock.top50comm)
head(mydata)

comm.melted<-melt(mydata)
class(comm.melted)
head(comm.melted)

# rough view #use reshape2 package
ggplot(comm.melted, aes(x=names, y=variable, fill=value),low="#000033", high="#CCFF66") + geom_tile()

tax_table(Pollock.rel.ps)
Pollock.myx <- subset_taxa(Pollock.rel.ps, family== "Myxococcaceae")
Pollock.myx
Pollock.myx@tax_table
Pollock.myx@tax_table[4:6]
colnames(Pollock.myx@sam_data)
title="Myxococcota by Temperature and Genus"
plot_bar(Pollock.myx, "temperature","Abundance", title=title)

tax_table(Pollock.rel.ps)
Pollock.vibrio <- subset_taxa(Pollock.rel.ps, family== "Vibrionaceae")
Pollock.vibrio
Pollock.vibrio@tax_table
Pollock.vibrio@tax_table[4:6]
colnames(Pollock.vibrio@sam_data)
title="Vibrio by Temperature and Genus"
plot_bar(Pollock.vibrio, "temperature","Abundance","genus", title=title)

tax_table(Pollock.rel.ps)
Pollock.endoz <- subset_taxa(Pollock.rel.ps, family== "Endozoicomonadaceae")
Pollock.endoz
Polloc.endoz@tax_table
Pollock.endoz@tax_table[4:6]
colnames(Pollock.endoz@sam_data)
title="Endozoicomonas by Temperature and Genus"
plot_bar(Pollock.endoz, "temperature","Abundance", "genus", title=title)

Pollock.comm<-Pollock.otu
Pollock.relAbundComm <- decostand(Pollock.comm, method = "total") #transforms to relative abundance

# check relative abundance in each sample (each should now be 1)
print(apply(Pollock.relAbundComm, 1, sum)) #yay
Pollock.relAbundComm[,1:5]
head(colSums(Pollock.relAbundComm)) #note, need to resort

Pollock.sortedRelAbundComm<-Pollock.relAbundComm[,order(-colSums(Pollock.relAbundComm))]
Pollock.top50comm<-Pollock.sortedRelAbundComm[,1:50]
dim(Pollock.sortedRelAbundComm)
dim(Pollock.top50comm)
Pollock.top50comm[1:5,1:5]

head(colSums(Pollock.top50comm))
head(colSums(Pollock.top50comm))

class(Pollock.top50comm)
Pollock.top50comm<-as.data.frame(Pollock.top50comm)
mydata<-Pollock.top50comm
mydata$names<-rownames(Pollock.top50comm)
head(mydata)

library(reshape2)
comm.melted<-melt(mydata)
class(comm.melted)
head(comm.melted)

# rough view #use reshape2 package
ggplot(comm.melted, aes(x=names, y=variable, fill=value)) + geom_tile()

Pollock.rel.ps@tax_table[3:13]