#ChungZiegler
#8/5/21
getwd()
setwd("/Users/jchung800/Desktop/R/Ziegler") 

# inspected Galaxy output with BBedit, it is tab-delimited
Ziegler.seqtab<-read.table("Ziegler.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
head(Ziegler.seqtab)

Ziegler.taxonomy<-read.table("Ziegler.taxonomy.tabular",header=TRUE,row.names=1,sep="\t")

colnames(Ziegler.seqtab)
# we need the metadata file to relate the colnames of the seqtab to the samples in the paper
# SRR5047388 is which sample?

# here is the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309854/#S1
# here is the bioproject: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=319637
# on the bioproject page: send to: file: format: runinfo
# downloaded run info as .csv, but need to annotate from paper
# inspected Galaxy output with BBedit, it is tab-delimited

# here is the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309854/#S1
# here is the bioproject: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=319637
# on the bioproject page: send to: file: format: runinfo
# downloaded run info as .csv, but need to annotate from paper

Ziegler.metadata<-read.csv("Ziegler.SraRunInfo.Edited.csv",header=TRUE,row.names=1)
row.names(Ziegler.metadata)

library(phyloseq)
OTU <- otu_table(Ziegler.seqtab,taxa_are_rows=TRUE)
SAMPLEDATA <- sample_data(Ziegler.metadata)
TAX <- tax_table(as.matrix(Ziegler.taxonomy))
Ziegler.ps<-phyloseq(OTU, SAMPLEDATA, TAX)

Ziegler.ps

# Remove taxa not seen more than 3 times in at least 20% of the samples. 
# According to Phyloseq: This protects against an OTU with small mean & trivially large C.V.
# However,  lot of times, might have something in less than 20% of samples
# Often better to filter in other more conservative ways

Ziegler.filter.ps = filter_taxa(Ziegler.ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
Ziegler.filter.ps
Ziegler.filter.ps = filter_taxa(Ziegler.ps, function(x) sum(x > 9) > (0.02*length(x)), TRUE)
Ziegler.filter.ps
saveRDS(Ziegler.filter.ps, "Ziegler.filter.ps.RDS")

colnames(Ziegler.filter.ps@tax_table) #to see column names of the taxa, note: family not capitalized
head(Ziegler.filter.ps@tax_table)
Endoz.ps = subset_taxa(Ziegler.filter.ps, family=="Endozoicomonadaceae") #subset a single family of Bacteria
Endoz.ps
title = "plot_bar; Endozoicomonadaceae-only"
colnames(Ziegler.filter.ps@sam_data) #to see the column names
plot_bar(Endoz.ps, "timeTreat", "Abundance", title=title)

# convert the sequences to ASV numbers: note, this step is just for easy visualization
library(Biostrings)
Endoz.ps
dna <- Biostrings::DNAStringSet(taxa_names(Endoz.ps))
names(dna) <- taxa_names(Endoz.ps)
Endoz.ps <- merge_phyloseq(Endoz.ps, dna)
taxa_names(Endoz.ps) <- paste0("Ziegler.endoz.ASV", seq(ntaxa(Endoz.ps)))
Endoz.ps


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

Endoz.metadata<-pssd2veg(Endoz.ps)
Endoz.otu<-psotu2veg(Endoz.ps)

head(Endoz.otu)

#8/6/21
# install.packages('vegan')
library(vegan)
comm<-Endoz.otu
relAbundComm <- decostand(comm, method = "total") #transforms to relative abundance

# check relative abundance in each sample (each should now be 1)
print(apply(relAbundComm, 1, sum)) #yay
relAbundComm[,1:5]
head(colSums(relAbundComm)) #note, need to resort

sortedRelAbundComm<-relAbundComm[,order(-colSums(relAbundComm))]
top100comm<-sortedRelAbundComm[,1:100]
dim(sortedRelAbundComm)
dim(top100comm)
top100comm[1:5,1:5]

head(colSums(top100comm))

#making a plot
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("dplyr")
library(dplyr)

class(top100comm)
top100comm<-as.data.frame(top100comm)
mydata<-top100comm
mydata$names<-rownames(top100comm)
head(mydata)

comm.melted<-melt(mydata)
class(comm.melted)
head(comm.melted)

# rough view #use reshape2 package
ggplot(comm.melted, aes(x=names, y=variable, fill=value)) + geom_tile()

ggplot(comm.melted, aes(x=variable, y=names, fill=sqrt(value))) + 
  geom_tile(color="white") +
  scale_fill_gradient(low="white",high="darkblue") +
  labs(fill = "Square\nRoot\nRelative\nAbundance") +
  labs(x = "OTU") +
  labs(y = "Sample") +
  theme(axis.text.x = element_text(size=5,angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(size=7,color='black')) +
  theme(axis.text.y=element_text(size=12,color='black')) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))

dim(comm.melted)
ggplot(comm.melted[1:920,], aes(x=variable, y=names, fill=sqrt(value))) + 
  geom_tile(color="white") +
  scale_fill_gradient(low="white",high="darkblue") +
  labs(fill = "Square\nRoot\nRelative\nAbundance") +
  labs(x = "OTU") +
  labs(y = "Sample") +
  theme(axis.text.x = element_text(size=5,angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(size=7,color='black')) +
  theme(axis.text.y=element_text(size=12,color='black')) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))

# how to aggregate by treatments, etc.?
melted<-psmelt(Endoz.ps)
head(melted)
orig.dest.time.trt.mean<-aggregate(melted$Abundance,by=list(melted$OTU,melted$orig_dest_time_treatment),FUN=mean)
head(orig.dest.time.trt.mean)
# a major problem with the psmelt approach is that it does not have our standarized relative abundance
# let's look at it anyway, and notice that psmelt did not sort the OTUs the way we want
ggplot(orig.dest.time.trt.mean, aes(x=Group.1, y=Group.2, fill=sqrt(x))) + 
  geom_tile(color="white") +
  scale_fill_gradient(low="white",high="darkblue") +
  labs(fill = "Square\nRoot\nRelative\nAbundance") +
  labs(x = "OTU") +
  labs(y = "Sample") +
  theme(axis.text.x = element_text(size=5,angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(size=7,color='black')) +
  theme(axis.text.y=element_text(size=12,color='black')) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))

head(melted$orig_dest_time_treatment)
row.names(melted)
row.names(Endoz.metadata)
row.names(Endoz.otu)
colnames(Endoz.otu)
colnames(Endoz.metadata)
class(Endoz.metadata)
Endoz.temp<-Endoz.otu

Endoz.temp$orig_dest_time_treatment<-Endoz.metadata$orig_dest_time_treatment

Endoz.temp$orig_dest_time_treatment<-""
head(Endoz.temp)
Endoz.agg<-aggregate(Endoz.otu,by=list(Endoz.metadata$orig_dest_time_treatment),FUN='mean')
head(Endoz.agg)

### investigate Pollock set123
setwd("/Users/thacker35294/Downloads")

raw.seqtab<-read.table("rawPollock123.dada2_sequencetable",header=TRUE,row.names=1,sep="\t")
dim(raw.seqtab) #32548 x 127  #note cut from 150 in this set to 127 survive
sum(rowSums(raw.seqtab)<20) #9827, note, we already cut <10 #a bunch of these are likely bimeras

hist(rowSums(Ziegler.seqtab),breaks=100000, xlim=c(0,5000))
# how many are more than 200?
sum(rowSums(raw.seqtab)>200) #2205
hist(rowSums(raw.seqtab),breaks=10000, xlim=c(0,100))
hist(rowSums(raw.seqtab),breaks=1000, xlim=c(200,1200), ylim=c(0,200))
max(rowSums(raw.seqtab)) #4030
hist(rowSums(raw.seqtab),breaks=3800, xlim=c(230,4030), ylim=c(0,200))
sum(rowSums(Ziegler.seqtab)<50) #22216
# being drastic might greatly reduce Jonathan's dataset size and make interpretation easier
sum(rowSums(raw.seqtab)>49) #10332
sum(rowSums(raw.seqtab)>99) #4839
sum(rowSums(raw.seqtab)>499) #825

library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(Ziegler.ps))
names(dna) <- taxa_names(Ziegler.ps)
Ziegler2017.ps <- merge_phyloseq(Ziegler.ps, dna)
taxa_names(Ziegler2017.ps) <- paste0("Ziegler2017.ASV", seq(ntaxa(Ziegler2017.ps)))
Ziegler2017.ps

hist(taxa_sums(Ziegler2017.ps),breaks=100000)

temp<-taxa_sums(Ziegler2017.ps)
class(temp)
sum(temp)
sum(temp>9)

sum(taxa_sums(Ziegler2017.ps>9))

hist(taxa_sums(Ziegler2017.ps),breaks=10000,xlim=c(0,200))
Ziegler2017.filter.ps =  filter_taxa(Ziegler2017.ps, function(x) sum(x > 149), TRUE)
Ziegler2017.filter.ps = filter_taxa(Ziegler2017.ps, function(x) sum(x > 149) > (0.001*length(x)), TRUE)
#transform to relative abundance
Ziegler2017.rel.ps = transform_sample_counts(Ziegler.filter.ps, function(x) x/sum(x))
#extract vegan format otu table and metadata
Ziegler.otu<-psotu2veg(Ziegler2017.rel.ps)
Ziegler.metadata<-pssd2veg(Ziegler2017.rel.ps)
 
#Bray-Curtis Dissimilarity
library(vegan)
Ziegler.comm.bc.dist <- vegdist(Ziegler.otu, method = "bray")
# cluster communities using average-linkage algorithm
Ziegler.comm.bc.clust <- hclust(Ziegler.comm.bc.dist, method = "average")
plot(Ziegler.comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis dissimilarity")

colnames(Ziegler.metadata)
adonis(comm.bc.dist ~ orig_dest_time_treatment, data = Ziegler.metadata, permutations = 99)


adonis(comm.bc.dist ~ originPool * destinationPool * treatment * time, data = Ziegler.metadata, permutations = 99)

adonis(comm.bc.dist ~ originPool + destinationPool + treatment + time, data = Ziegler.metadata, permutations = 999)
simperOutput<-simper(Ziegler.otu,Ziegler.metadata$treatment,perm=9)

require("ggplot2")
p2 = plot_ordination(Ziegler2017.filter.ps, ordinate(Ziegler2017.filter.ps, "CCA"), type = "samples", color = "treatment")
p2 + geom_point(size = 5) + geom_polygon(aes(fill = treatment))


#8/23/21
colnames(Ziegler2017.rel.ps@sam_data) #to see the column names
#subsetting the Endozoicomonadaceae
Ziegler.endoz.ps=  subset_taxa(Ziegler2017.rel.ps, family=="Endozoicomonadaceae") 
Ziegler.endoz.ps #108 taxa
title = "plot_bar; Ziegler Endozoicomonadaceae-only"
plot_bar(Ziegler.endoz.ps, "orig_dest_time_treatment", "Abundance", title=title)
#subset the alphas 
Ziegler.alpha.ps=  subset_taxa(Ziegler2017.rel.ps, class=="Alphaproteobacteria") 
Ziegler.alpha.ps #1922 taxa
title = "plot_bar; Ziegler Alphaproteobacteria-only"
plot_bar(Ziegler.alpha.ps, "orig_dest_time_treatment", "Abundance", title=title)

#subset the gammas 
Ziegler.gamma.ps=  subset_taxa(Ziegler2017.rel.ps, class=="Gammaproteobacteria") 
Ziegler.gamma.ps #108 taxa
title = "plot_bar; Ziegler Gammaproteobacteria-only"
plot_bar(Ziegler.gamma.ps, "orig_dest_time_treatment", "Abundance", title=title)

library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(Ziegler.filter.ps))
names(dna) <- taxa_names(Ziegler.filter.ps)
Ziegler2019.filter.ps <- merge_phyloseq(Ziegler.filter.ps, dna)
taxa_names(Ziegler2019.filter.ps) <- paste0("Ziegler2019.ASV", seq(ntaxa(Ziegler2019.filter.ps)))
Ziegler2019.filter.ps
head(taxa_names(Ziegler2019.filter.ps))
#transform to relative abundance
Ziegler.rel.ps = transform_sample_counts(Ziegler.filter.ps, function(x) x/sum(x)) 
#extract vegan format otu table and metadata
Ziegler.otu<-psotu2veg(Ziegler.rel.ps)
Ziegler.metadata<-pssd2veg(Ziegler.rel.ps)

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
library(vegan)
comm.bc.dist <- vegdist(Ziegler.otu, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
plot(comm.bc.clust, cex = 0.7, ylab = "Bray-Curtis Dissimilarity")

phyloseq::plot_bar(Ziegler2017.rel.ps, fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H


phyloseq::plot_bar(Ziegler.rel.ps, fill = "family") + 
  geom_bar(aes(color = family, fill = family), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ orig_dest_time_treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_viridis(option = "H",discrete=TRUE) + #options A to H #viridis has different color sets and this function allows you to customize
  scale_fill_viridis(option = "H",discrete=TRUE) #options A to H



#Ordination
library(vegan) 
library(MASS)
Ziegler.dis <- Ziegler.comm.bc.dist
Ziegler.mds0 <- isoMDS(Ziegler.dis)
stressplot(Ziegler.mds0, Ziegler.dis)
ordiplot(Ziegler.mds0, type='t')
Ziegler.mds <- metaMDS(Ziegler.otu)
ef <- envfit(Ziegler.mds, smallerZiegler, permu=99)
dim(Ziegler.otu)
head(colSums(Ziegler.otu))
smallerZiegler= Ziegler.otu[,1:1000]
plot(ef)

plot(Ziegler.mds, display="sites")
plot(ef, p.max=0.1)

library(phyloseq)
Ziegler.pcoa <- pcoa(Ziegler.dis, correction= "none", rn=NULL)
plot(Ziegler.pcoa)
biplot(Ziegler.pcoa, Y=NULL, plot.aces= c(1,2), dir.axis1=1, dir.axis2=1)

#diveristy
Ziegler.comm<-Ziegler.otu
Ziegler.metadata<-Ziegler.metadata
colnames(Ziegler.metadata)
# compare species Simpson diversity between host species
Ziegler.invsimpson.comm <- diversity(Ziegler.comm, index = "invsimpson")
Ziegler.shannon <- diversity(Ziegler.comm, index = "shannon")
boxplot(Ziegler.invsimpson.comm ~ Ziegler.metadata$treatment, xlab = "Treatment", ylab = "InvSimpson Index")
boxplot(Ziegler.shannon ~ Ziegler.metadata$treatment, xlab = "Treatment", ylab = "Shannon Index")


comm<-Ziegler.otu
metadata<-Ziegler.metadata
colnames(metadata)
# compare species Simpson diversity between host species
invsimpson.comm <- diversity(comm, index = "invsimpson", MARGIN=1)

boxplot(invsimpson.comm ~ metadata$treatment, xlab = "treatment", ylab = "InvSimpson Index")
comm<-Yu2020.otu
metadata<-Yu2020.metadata
colnames(metadata)
# compare species Simpson diversity between host species
invsimpson.comm <- diversity(comm, index = "invsimpson", MARGIN= 1)
boxplot(invsimpson.comm ~ metadata$Treatment, xlab = "Treatment", ylab = "InvSimpson Index")