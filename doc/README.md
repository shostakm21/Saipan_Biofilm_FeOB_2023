
# Title: "Saipan Dada2 Pipeline"
## Author: "Maggie Shostak"

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegan")
install.packages("devtools")
```

```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
library(rmarkdown)
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan) ; packageVersion("vegan") # 2.5.4
library(dbplyr)
library(microbiome)
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
```

# Format P-values
```{r}
formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}
```

```{r}
# Path needs to have fastq files unzipped!
path <- "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/fastq"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
```

```{r}
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
```

```{r}
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
```

## Merge Paired Reads
We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region.

```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
head(mergers[[1]])
```

**Most reads should pass the merging step! If that isn't the case, are you sure your truncated reads still overlap sufficiently?**

## Construct Sequence Table (ASV Table)
The sequence table is a `matrix` with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

  #451 
#171813 
```

The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. The lengths of the merged sequences all fall in the expected range for this amplicon.

## Remove chimeras
 
The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences.
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#[1]   141 73262
#[1] 0.9262261

# 97 67937
# 0.9241137
```

The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on factors including experimental procedures and sample complexity. Here chimeras make up about X% of the merged sequence variants, but when we account for the abundances of those variants we see they account for only about X% of the merged sequence reads.
**In some cases, most sequences will be chimric. But most reads should not be.**

## Track reads through the pipeline: Look at the number of reads that made it through each step in the pipeline
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#write.table(track, "/Users/maggieshostak/Desktop/Saipan_R_Studio/track_saipan_sequences_original.csv", sep=",", quote=F, col.names=NA)
write.table(track, "/Users/maggieshostak/Desktop/Saipan_R_Studio/track_saipan_sequences_rarefied.csv", sep=",", quote=F, col.names=NA)
```

## Assign Taxonomy
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/fastq/silva_nr99_v138.1_train_set.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to assign taxonomy to the sequence variants. The DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose. The assignTaxonomy function takes as input a set of sequences to be classified, and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.

We maintain formatted training fastas for the RDP training set, GreenGenes clustered at 97% identity, and the Silva reference database, and additional trainings fastas suitable for protists and certain specific environments have been contributed. For fungal taxonomy, the General Fasta release files from the UNITE ITS database can be used as is. To follow along, download the silva_nr_v128_train_set.fa.gz file, and species assignment file, and place it in the directory with the fastq files. https://benjjneb.github.io/dada2/training.html

The dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identitiy) is the only appropriate way to assign species to 16S gene fragments. Currently, species-assignment training fastas are available for the Silva and RDP 16S databases. To follow the optional species addition step, download the silva_species_assignment_v128.fa.gz file, and place it in the directory with the fastq files.

```{r}
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

```{r}
#write(asv_fasta, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_fasta_saipan.fa")
#write.table(asv_otu, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(asv_tax, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(otu_tax_table, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_tax_table_saipan.csv", sep=",", quote=F, col.names=NA)
```

# Rarefaction: Check for low sequence counts before processing through various analyses
```{r}
library(phyloseq) ; packageVersion("phyloseq") # 1.22.3
library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
library(viridis) ; packageVersion("viridis") # 0.5.1
library(ggpubr)
```

```{r}
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv")
otu_tab
```

```{r}
otu_counts <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts
```

```{r}
taxonomy <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv")
taxonomy
```

```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/data/metadata_saipan.csv")
metadata
```

# Rarefaction: ALREADY COMPLETED
```{r}
#deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab, colData = metadata, design = ~sample_id) 
#deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
#deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
#vst_trans_count_tab <- assay(deseq_counts_vst)
#euc_dist <- dist(t(vst_trans_count_tab))

#euc_clust <- hclust(euc_dist, method="ward.D2")
#plot(euc_clust) 
#euc_dend <- as.dendrogram(euc_clust, hang=0.1)
#dend_cols <- as.character(metadata$color[order.dendrogram(euc_dend)])
#labels_colors(euc_dend) <- dend_cols

#plot(euc_dend, ylab="VST Euc. dist.", width = 35, height = 35 )

#rarecurve(t(count_tab), step=500, col=metadata$color, lwd=2, ylab="ASVs", label=F, cex=0.6)
  #abline(v=(min(rowSums(t(count_tab)))), col = "red")
```

## OTU Relative Abundance Generation
```{r}
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_rel_abund_saipan.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_rel_abund_saipan.csv")
otu_rel_abund
```

# Stacked Barcharts
## All Samples
```{r}
## Phylum
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/phylum_stacked_barchart_saipan_all.tiff", width=20, height=10)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/class_stacked_barchart_saipan_all.tiff", width=25, height=10)

## Order
otu_rel_abund %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/order_stacked_barchart_saipan_all.tiff", width=55, height=10, limitsize = FALSE)

## Family
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/family_stacked_barchart_saipan_all.tiff", width=50, height=10, limitsize = FALSE)
```

# Non-multidimensional Scaling (NMDS Plots): All Samples
```{r}
df_meta <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/data/metadata_saipan.csv")
df_meta

df1_otu <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_1.csv")
df1_otu

df2_otu <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_2.csv")
df2_otu

df3_otu <- read.csv("//Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_3.csv")
df3_otu

df_otu_all <- list(df1_otu, df2_otu, df3_otu)
df_otu_all

write.table(df_otu_all,"/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_saipan.csv", sep=",", col.names=NA)

df_otu_all <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_saipan.csv")
df_otu_all

nmds_asv_otu_all<- inner_join(df_meta, df_otu_all, by="sample_id")
nmds_asv_otu_all

write.table(nmds_asv_otu_all, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/nmds_asv_otu_all.csv", sep=",", quote=F, col.names=NA)
```

```{r}
pc <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/nmds_asv_otu_all.csv")
pc
```

```{r}
#make community matrix: extract columns with ASV information
com <- pc[,6:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)
```

```{r}
#Run NMDS using Bray-Curtis distance
set.seed(123)
nmds <- metaMDS(m_com, distance="bray") #stress = 0.1244532 
nmds
plot(nmds)

#access the specific points data of the NMDS plot & scores
str(nmds)
nmds$points
scores(nmds)

#extract NMDS scores
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame
data.scores$sample_id = pc$sample_id
data.scores$location = pc$location
data.scores$depth = pc$depth
data.scores$sample_type = pc$sample_type

head(data.scores)
```

## All Samples: By Location
```{r}
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/NMDS_all_samples_location.tiff", width = 10, height = 10)
```

## All Samples: Bioiflm vs Sed vs Water
```{r}
xx1 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = sample_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Sample Type", y = "NMDS2")
xx1
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/NMDS_all_samples_sample_type.tiff", width = 10, height = 10)
```

# Statistical Analyses
## Anosim
```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/nmds_asv_otu_all.csv")
pc_ano
```

```{r}
# All Samples: Biof vs Sed vs Water by Location
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

#ANOSIM statistic R:   0.2 
      #Significance: 1e-04 
```

```{r}
# All Samples: Biof vs Sed vs Water by Sample Type
com_ano1 = pc_ano[,6:ncol(pc_ano)]
m_com_ano1 = as.matrix(com_ano1)
ano_all1 = anosim(m_com_ano1, pc_ano$sample_type, distance = "bray", permutations = 9999)
ano_all1

#ANOSIM statistic R: 0.9018 
      #Significance: 1e-04 
```

```{r}
# All Samples: Biof vs Sed vs Water by Depth
com_ano2 = pc_ano[,6:ncol(pc_ano)]
m_com_ano2 = as.matrix(com_ano2)
ano_all2 = anosim(m_com_ano2, pc_ano$depth, distance = "bray", permutations = 9999)
ano_all2

#ANOSIM statistic R: 0.1617 
      #Significance: 1e-04 
```

# Diversity Index Value Generation
```{r}
## Biofilm & Sediment & Water Samples
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv", header=T, row.names=1, check.names=FALSE)

## Transpose the data to have sample names on rows
otu.table.diver <- t(otu_tab)
otu.table.diver <- as.data.frame(otu.table.diver)
head(otu.table.diver)
```

```{r}
otu_counts <- Reduce(function(x, y) merge(x, y, all=TRUE), df_otu_all) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")

write.table(otu_counts, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_saipan_all.csv", sep=",", quote=F, col.names=NA)

otu_count_all %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

## Shannon H Diveristy, Species Richness & Pielou Evenness
```{r}
data(otu.table.diver)
H<-diversity(otu.table.diver)
H

richness <- specnumber(otu.table.diver)
richness

evenness <- H/log(richness)
evenness

alpha <- cbind(shannon = H, richness = richness, pielou = eveness, metadata)
write.csv(alpha, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/diversity_indices_saipan_all.csv")
head(alpha)
```

### Boxplots
```{r}
plot.shan <- ggplot(alpha, aes(x = location, y = shannon, colour = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Location_bio_sed_water.tiff")

plot.rich <-ggplot(alpha, aes(x = location, y = richness, colour = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Richness_Location_bio_sed_water.tiff")

plot.even <- ggplot(alpha, aes(x = location, y = pielou, colour = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Pielou's_Evenness_Location_bio_water.tiff")

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_bio_sed_water.tiff")
```
