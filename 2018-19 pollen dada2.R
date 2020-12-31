#=============
#Getting ready
#=============
library(dada2); packageVersion("dada2")
path <- "/home/heuklang/R projects/2018_2019_pollen/raw"
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="1_[0-9]{4,4}", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="2_[0-9]{4,4}", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- gsub("/home/heuklang/R projects/2018_2019_pollen/raw/","",fnFs)
sample.names <- gsub(".fastq.gz","",sample.names)
sample.names <- gsub("_[0-9]_","",sample.names)
sample.names <- gsub("-","",sample.names)
sample.names <- gsub("2018","-2018",sample.names)
sample.names <- gsub("2019","-2019",sample.names)

#==============================
# Inspect read quality profiles
#==============================

plotQualityProfile(fnFs[1:24]) #280bp
plotQualityProfile(fnRs[1:24]) #200bp

#================
# Filter and trim
#================
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,200),
                     trimLeft=c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

?filterAndTrim
out

plotQualityProfile(filtFs[1:24])
plotQualityProfile(filtRs[1:24])

#======================
# Learn the Error Rates
#======================
errF <- learnErrors(filtFs, multithread=TRUE, 
                    nbases= 1e9)
errR <- learnErrors(filtRs, multithread=TRUE,
                    nbases= 1e9)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#=====================
# Sample infence
#===================
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#==================
#Merge paired reads
#==================
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#================
# Remove chimeras
#================
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


#==================================
# Track reads through the pipeline
#==================================
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#=================
# Assign taxonomy
#=================
tax <- assignTaxonomy(seqtab.nochim, "/home/heuklang/R projects/training/silva_nr_v132_train_set.fa.gz", multithread=T)
tax <- addSpecies(tax, "/home/heuklang/R projects/training/silva_species_assignment_v132.fa.gz")  
write.csv(tax, "/home/heuklang/R projects/2018_2019_pollen/no filtered table.csv")

results <- cbind(tax, t(seqtab.nochim))
results[is.na(results)] <- "Not assigned"
results <- as.data.frame(results)
filteredRs <- results[results$Kingdom == "Bacteria" | results$Kingdom == "Archaea",]
filteredRs <- filteredRs[filteredRs$Order != "Chloroplast",]
filteredRs <- filteredRs[filteredRs$Family != "Mitochondria",] 
write.csv(filteredRs, "/home/heuklang/R projects/2018_2019_pollen/Only bacteria.csv")
filteredRs <- read.csv("/home/heuklang/R projects/2018_2019_pollen/Only bacteria.csv")
filteredRs
#==================
# for  PiCrust2
#seqeunce fasta file output
seq <- as.vector(filteredRs$X)

#install. packages("foreach")
library(foreach)
#install.packages("doParallel")
library(doParallel)

registerDoParallel(cores=16) #최대 자신의 컴퓨터 코어 수에 맞게..., R 외에 여러가지 일을 하고 있으면 최대수-1을 추천
OTUsID <- foreach(i = 1:length(seq), .combine="rbind")%dopar%{
  OTUsN <- paste("OTU", i, sep="_")
  print(OTUsN)
}
OTUsID <- as.vector(OTUsID)

fasta <- foreach(i= 1:length(seq), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], seq[i])
  print(match)
}
fasta <- as.vector(fasta)
fasta <- as.data.frame(gsub("^OTU_", ">OTU_",fasta))
colnames(fasta) <- NA
fasta
write.table(fasta, file = "/home/heuklang/R projects/2018_2019_pollen/OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

#OTU counts tsv file output 
OTUtableRaw <- read.csv("/home/heuklang/R projects/2018_2019_pollen/Only bacteria.csv", header=T)
OTUtable <- OTUtableRaw[,-(1:8)]
head(OTUtable)
OTU_counts <- cbind(OTUsID,OTUtable)
write.table(OTU_counts, file = "/home/heuklang/R projects/2018_2019_pollen/OTUs_counts.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#==================
# rarefaction curve
#==================
library(ggplot2)
library(iNEXT)
install.packages("gridExtra")
library(gridExtra)

filtered_reads <- read.csv("/home/heeildo/R works/2018_2019_pollen/Only bacteria.csv")
filtered_reads[9:32]
iOF <- filtered_reads[9:32]

rareiOF_18apple <- iNEXT(iOF[4:6], q=0, datatype = "abundance")
rareiOF_18kiwi <- iNEXT(iOF[10:12], q=0, datatype = "abundance")
rareiOF_18peach <- iNEXT(iOF[16:18], q=0, datatype = "abundance")
rareiOF_18pear <- iNEXT(iOF[22:24], q=0, datatype = "abundance")
rareiOF_19apple <- iNEXT(iOF[1:3], q=0, datatype = "abundance")
rareiOF_19kiwi <- iNEXT(iOF[7:9], q=0, datatype = "abundance")
rareiOF_19peach <- iNEXT(iOF[13:15], q=0, datatype = "abundance")
rareiOF_19pear <- iNEXT(iOF[19:21], q=0, datatype = "abundance")

Replication <- as.character(c(1:3,2))
Method <- c("interpolated","interpolated","extrapolated","extrapolated")
df_legend<- data.frame(Replication,Method)
df_legend$Method <- factor(Method, labels = c("interpolated","extrapolated"))
leg<-ggplot(df_legend,aes(Replication, Method, shape=Replication, linetype=Method, group=Method))+
  geom_point(size=4.5)+
  geom_line(size=1.2)+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.key.width = unit(1,"cm"))

legend <- cowplot::get_legend(leg)

rarefaction_18apple <- ggiNEXT(rareiOF_18apple, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#F7756C","#F7756C","#F7756C"))+
  scale_fill_manual(values=c("#F7756C","#F7756C","#F7756C"))+
  ggtitle("2018 apple pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_18kiwi <- ggiNEXT(rareiOF_18kiwi, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#7BAD00","#7BAD00","#7BAD00"))+
  scale_fill_manual(values=c("#7BAD00","#7BAD00","#7BAD00"))+
  ggtitle("2018 kiwi pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_18peach <- ggiNEXT(rareiOF_18peach, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#00BFC4","#00BFC4","#00BFC4"))+
  scale_fill_manual(values=c("#00BFC4","#00BFC4","#00BFC4"))+
  ggtitle("2018 peach pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_18pear <- ggiNEXT(rareiOF_18pear, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#C67BFF","#C67BFF","#C67BFF"))+
  scale_fill_manual(values=c("#C67BFF","#C67BFF","#C67BFF"))+
  ggtitle("2018 pear pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_19apple <- ggiNEXT(rareiOF_19apple, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#F7756C","#F7756C","#F7756C"))+
  scale_fill_manual(values=c("#F7756C","#F7756C","#F7756C"))+
  ggtitle("2019 apple pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_19kiwi <- ggiNEXT(rareiOF_19kiwi, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#7BAD00","#7BAD00","#7BAD00"))+
  scale_fill_manual(values=c("#7BAD00","#7BAD00","#7BAD00"))+
  ggtitle("2019 kiwi pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_19peach <- ggiNEXT(rareiOF_19peach, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#00BFC4","#00BFC4","#00BFC4"))+
  scale_fill_manual(values=c("#00BFC4","#00BFC4","#00BFC4"))+
  ggtitle("2019 peach pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))

rarefaction_19pear <- ggiNEXT(rareiOF_19pear, type=1)+
  xlab("Reads per samples")+
  ylab("OTUs number")+
  scale_color_manual(values=c("#C67BFF","#C67BFF","#C67BFF"))+
  scale_fill_manual(values=c("#C67BFF","#C67BFF","#C67BFF"))+
  ggtitle("2019 pear pollen")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits=c(0,150))+
  scale_x_continuous(limits=c(0,130000),breaks=c(0:4*30000))    

rarefaction <- grid.arrange(rarefaction_18apple, rarefaction_18kiwi, rarefaction_18peach, rarefaction_18pear, 
                            rarefaction_19apple, rarefaction_19kiwi, rarefaction_19peach, rarefaction_19pear, 
                            nrow=2)
rarefaction_and_legend  <- grid.arrange(rarefaction, legend, nrow=2,
                                        heights= c(15,1))

#==============
# Diversity
#==============

install.packages('devtools')
library(devtools)  # Load the devtools package
install_github("umerijaz/microbiomeSeq")  # Install the package

library(microbiomeSeq)
install.packages('microbiomeSeq')

library(microbiomeSeq)  #load the package
library(phyloseq)
#====================
# Bargraph, before remove plant
#====================
before_remove_plant <- read.csv("/home/heeildo/R works/2018_2019_pollen/Only bacteria.csv")
head(before_remove_plant)
OTU <- otu_table(as.matrix(t(before_remove_plant[9:32])), taxa_are_rows = FALSE) #OTU table
TAX <- tax_table(as.matrix(before_remove_plant[2:8]))

sampleID <- colnames(before_remove_plant[9:32])
years <- gsub("^X","", sampleID)
years <- gsub("[0-9]$","", years)
years <- gsub("[A-Za-z]","", years)
years <- gsub("^","20",years)

crops <- gsub("X[0-9][0-9]", "", sampleID)
crops <- gsub("[0-9]", "", crops)

meta_table <- data.frame(sampleID,years,crops)
row.names(meta_table) <- meta_table$sampleID

SAM <- sample_data(meta_table)
phyloseq <- phyloseq(OTU, TAX, SAM)

physeq_non_remove <- phyloseq
physeq_non_remove <- taxa_level(physeq_non_remove, "Kingdom")
physeq_non_remove <- normalise_data(physeq_non_remove, norm.method = "relative")
plot_taxa(physeq_non_remove, grouping_column = "crops", method = "hellinger", number.taxa = 4, filename = NULL)+
  scale_y_continuous(labels = scales::percent)+
  xlab("")+
  ylab("Percantage of population")

#====================
# Bargraph, after remove plant
#====================
after_remove_plant <- read.csv("/home/heuklang/R projects/pollen/fungi/fungi.nochim.taxa.csv")
head(before_remove_plant)
OTU <- otu_table(as.matrix(t(after_remove_plant[9:32])), taxa_are_rows = FALSE) #OTU table
TAX <- tax_table(as.matrix(after_remove_plant[2:8]))

sampleID <- colnames(after_remove_plant[9:32])
years <- gsub("^X","", sampleID)
years <- gsub("[0-9]$","", years)
years <- gsub("[A-Za-z]","", years)
years <- gsub("^","20",years)

crops <- gsub("X[0-9][0-9]", "", sampleID)
crops <- gsub("[0-9]", "", crops)

meta_table <- data.frame(sampleID,years,crops)
row.names(meta_table) <- meta_table$sampleID

SAM <- sample_data(meta_table)
phyloseq <- phyloseq(OTU, TAX, SAM)

physeq <- phyloseq
physeq <- taxa_level(physeq, "Genus")
physeq <- normalise_data(physeq, norm.method = "relative")
?plot_taxa
plot_taxa(physeq, grouping_column = "crops", method = "chord", number.taxa = 10, filename = NULL)+
  scale_y_continuous(labels = scales::percent)+
  xlab("")+
  ylab("Percantage of population")


#=====
# NMDS
#=====
library(vegan)
library(grid)
library(gridExtra)
physeq <- phyloseq
ord.res <- ordination(physeq, which_distance = "bray", method = "NMDS", grouping_column = "crops", 
                      pValue.cutoff = 0.05)

plot.ordination(ord.res, method = "NMDS", extra_marginspace = 0.2, show.pvalues = T)

ord.res <- ordination(physeq, which_distance = "bray", method = "PCoA", grouping_column = "crops", 
                      pValue.cutoff = 0.05)
ord.res

plot.ordination(ord.res, method = "PCoA", extra_marginspace = 0.2, show.pvalues = T)

plot_cca
plot_cca2(physeq = physeq, grouping_column = "crops", pvalueCutoff = 0.01, 
          env.variables = NULL, num.env.variables = NULL, exclude.variables = "crops", 
          draw_species = T)

#======================
# plot_anovadiversity
#======================
physeq <- phyloseq
plot_anova_diversity(physeq, method = c("richness", "simpson", "shannon"), 
                     grouping_column = "crops", pValueCutoff = 0.05)


