
# Title: "Saipan Dada2 Pipeline"
## Author: "Maggie Shostak"

# Load all necessary packages: 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegan")
install.packages("devtools")
BiocManager::install("microbiome")
```

```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
library(BiocManager)
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
library(conflicted)
library(dplyr)
library(cowplot)
```

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

# DADA2 Pipeline
```{r}
# Path needs to have fastq files unzipped!
path <- "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/fastq"
list.files(path)
```

```{r}
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
```

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

```{r}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
```

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])
```

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Identified 8786 bimeras out of 43196 input sequences.
#[1]   134 34410
#[1] 0.9768214
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/track_sequences.csv", sep=",", quote=F, col.names=NA)
```

## Taxonomic Classification & ASV File Formatting
```{r}
taxa <- assignTaxonomy(seqtab.nochim,"/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/fastq/silva_nr99_v138.1_train_set.fa.gz")
```

```{r}
taxa.print <- taxa
rownames(taxa.print) <- NULL
#head(taxa.print)
```

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
write.table(asv_fasta, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_fasta_saipan.fa")
write.table(asv_otu, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_tax_saipan.csv", sep=",", quote=F, col.names=NA)
```

# Formating Files for Future Analysis
```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/data/metadata_saipan.csv")
metadata

otu_counts <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv")
taxonomy
```

# Phyloseq
```{r}
library(phyloseq)
library(Biostrings)
```

```{r}
# Read data into R
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv")
otu_tab

tax_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv")
tax_tab

samples_df <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/data/metadata_saipan.csv")
samples_df
```

```{r}
# Phyloseq objects need to have row.names
otu_tab <- otu_tab %>%
  tibble::column_to_rownames("ASV")

tax_tab <- tax_tab %>%
  tibble::column_to_rownames("ASV")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample_id")

# Transform OTU & Tax table into matrices
otu_tab <- as.matrix(otu_tab)
tax_tab <- as.matrix(tax_tab)

#Transform into Phyloseq Objects
ASV = otu_table(otu_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
sample_id = sample_data(samples_df)
  
ps <- phyloseq(ASV, TAX, sample_id)
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 33035 taxa and 97 samples ]
#sample_data() Sample Data:       [ 97 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 33035 taxa by 6 taxonomic ranks ]

# Visualize Data
sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)

# Bar graphs based on division
plot_bar(ps, fill = "Phylum")

plot_bar(ps, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Heatmaps
plot_heatmap(ps, method = "NMDS", distance = "bray")

ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.20) > 0, TRUE)
ps_abund

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6 taxa and 134 samples ]
#sample_data() Sample Data:       [ 134 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 6 taxa by 6 taxonomic ranks ]

otu_table(ps)[1:8, 1:5]

plot_heatmap(ps_abund, method = "NMDS", distance = "bray")
```

```{r}
# Alpha Diversity: Plot Chao1 richness estimator and Shannon diversity estimator
plot_richness(ps, measures=c("Chao1", "Shannon"))

# Regroup together from the same location
plot_richness(ps, measures=c("Chao1", "Shannon"), x="location", color="sample_type") +
  geom_boxplot()
ggsave("/Users/maggieshostak/Desktop/RStudio_Saipain_Data/results/Chao1_Shannon_Saipan_All_Location.tiff", width=15, height=10)

## Just Shannon
plot_richness(ps, measures=c("Shannon"), x="location", color="sample_type") +
  geom_boxplot()
ggsave("/Users/maggieshostak/Desktop/RStudio_Saipain_Data/results/Chao1_Shannon_Saipan_All_Location.tiff", width=15, height=10)
```

```{r}
# Ordination
ps.ord <- ordinate(ps, "NMDS", "bray")
plot_ordination(ps, ps.ord, type="taxa", color="Phylum", shape= "Location", 
                  title="Saipan ASVs")

plot_ordination(ps, ps.ord, type="taxa", color="Phylum", 
                  title="Sapain ASVs", label="Phylum") + 
                  facet_wrap(~Phylum, 3)
#ggsave("Ordination_Phylum_Saipan_All.tiff", width=25, height=10)

plot_ordination(ps, ps.ord, type="samples", color="Location", 
                  shape="depth", title="Samples") + geom_point(size=3)

plot_ordination(ps, ps.ord, type="split", color="Phylum", 
                  shape="Location", title="Biplot", label = "station") +  
  geom_point(size=3)
```

# Rarefaction
```{r}
#library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
#library(vegan) ; packageVersion("vegan") # 2.5.4
#library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
#library(viridis) ; packageVersion("viridis") # 0.5.1
#library(ggplot2)
#library(devtools)
#library(ggpubr)
```

```{r}
#count_tab <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan_original/data/asv_otu/asv_otu_saipan.csv", header=T, row.names=1, check.names=F, sep=",")
#count_tab

#tax_tab <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan_original/data/asv_tax/asv_tax_saipan.csv", header=T, row.names=1, check.names=F, sep=",")
#tax_tab

#metadata <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan_original/data/metadata/metadata_saipan.csv")
#metadata
```

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
```

```{r}
#rarecurve(t(count_tab), step=500, lwd=2, ylab="ASVs", label=T, cex=0.6)
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

# All Samples: Classified as Ship (Biofilm), Sediment, Water Only [SSW]
```{r}
metadata_SSW <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/metadata_SSW.csv")
metadata_SSW

otu_counts_SSW <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_SSW

taxonomy_SSW <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_tax_saipan.csv")
taxonomy_SSW
```

```{r}
otu_rel_abund_SSW <- inner_join(metadata_SSW, otu_counts_SSW, by="sample_id") %>%
  inner_join(., taxonomy_SSW, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_SSW

write.table(otu_rel_abund_SSW,"/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_rel_abund_SSW.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund_SSW <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_rel_abund_SSW.csv")
```

## Barchart
```{r}
## Phylum
otu_rel_abund_SSW %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, sample_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/phylum_stacked_barchart_SSW.tiff", width=20, height=15)

## Class
otu_rel_abund_SSW %>%
  filter(level=="Class") %>%
  group_by(sample_id, sample_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/class_stacked_barchart_SSW.tiff", width=30, height=20)

## Order
otu_rel_abund_SSW %>%
  filter(level=="Order") %>%
  group_by(sample_id, sample_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/order_stacked_barchart_SSW.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_SSW %>%
  filter(level=="Family") %>%
  group_by(sample_id, sample_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/family_stacked_barchart_SSW.tiff", width=55, height=20, limitsize = FALSE)
```

# Non-multidimensional Scaling (NMDS Plots): All Samples
```{r}
df_meta <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/metadata_saipan.csv")
df_meta

df_otu <-read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu.csv", row.names=1)
df_otu

nmds_asv_otu_all <- inner_join(df_meta, df_otu, by="sample_id")
nmds_asv_otu_all

write.table(nmds_asv_otu_all, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/nmds_asv_otu_all.csv", sep=",", quote=F, col.names=NA)
```

```{r}
pc <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/nmds_asv_otu_all.csv")
pc

#make community matrix: extract columns with ASV information
com <- pc[,7:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)

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

# ANOSIM statistic R: 0.2
      # Significance: 1e-04
```

```{r}
# All Samples: Biof vs Sed vs Water by Sample Type
com_ano1 = pc_ano[,6:ncol(pc_ano)]
m_com_ano1 = as.matrix(com_ano1)
ano_all1 = anosim(m_com_ano1, pc_ano$sample_type, distance = "bray", permutations = 9999)
ano_all1

# ANOSIM statistic R: 0.9018 
      # Significance: 1e-04 
```

```{r}
# All Samples: Biof vs Sed vs Water by Depth
com_ano2 = pc_ano[,6:ncol(pc_ano)]
m_com_ano2 = as.matrix(com_ano2)
ano_all2 = anosim(m_com_ano2, pc_ano$depth, distance = "bray", permutations = 9999)
ano_all2

# ANOSIM statistic R: 0.1617 
      # Significance: 1e-04 
```

# Ecological Distance Matrices
```{r}
library(tidyverse)

df1 <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_1.csv")
df1

df2 <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_2.csv")
df2

df3 <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_3.csv")
df3

df_otu <- list(df1, df2, df3)
df_otu

write.table(df_otu,"/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_saipan.csv", sep=",", col.names=NA)
```

```{r}
df_otu <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_saipan.csv")
df_otu <- subset(df_otu, select = -c(X))
df_otu

otu_count_all <- df_otu %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_all

write.table(otu_count_all, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_all_saipan.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_count_all <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_all_saipan.csv")

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

## Writing functions for Shannon, Richness & Simpson
```{r}
richness <- function(x){
  sum(x > 0)}

shannon <- function(x){
  relabund <- x[x>0]/sum(x)
  -sum(relabund * log(relabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) /(n*n-1))
}
```

## Diversity Metrics
```{r}
otu_count_all <- otu_count_all %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count_all, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_all_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

```{r}
diversity_metrics <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_all_diversity_metrics.csv")
diversity_metrics

write.table(diversity_metrics, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/saipan_all_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

## Plot Alpha Diversity Functions
```{r}
diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, color= sample_type)) +
  geom_boxplot() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/alpha_diversity_metrics_all_sample_type.tiff", width = 10, height = 15)

#Each point represents a sample, (n) Sum of Count, (X) Total number of sequences for each sample & (Y) Value of diversity metric
```

# Beta Diversity Metrics: Boxplots of Shannon, Evenness, and Pielou
```{r}
# Shannon
data(otu.table.diver.all)
H <- diversity(otu.table.diver.all)
H

richness <- specnumber(otu.table.diver.all)
richness

evenness <- H/log(richness)
evenness

metadata_all <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/metadata_saipan.csv")
metadata_all

alpha <- cbind(shannon = H, richness = richness, pielou = evenness, metadata_all)
write.csv(alpha, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/diversity_indices_bio_sed_water.csv")
head(alpha)
```

## Plot Diversity Boxplots: Sample Location
```{r}
plot.shan <- ggplot(alpha, aes(x = location, y = shannon, colour = location)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Location_bio_sed_water.tiff")

plot.rich <-ggplot(alpha, aes(x = location, y = richness, colour = location)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Richness_Location_bio_sed_water.tiff")

plot.even <- ggplot(alpha, aes(x = location, y = pielou, colour = location)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Pielou's_Evenness_Location_bio_water.tiff")

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_all.tiff")
```

## Plot Diversity Boxplots: Sample Type
```{r}
plot.shan <- ggplot(alpha, aes(x = sample_type, y = shannon, colour = sample_type)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Location_bio_sed_water.tiff")

plot.rich <-ggplot(alpha, aes(x = sample_type, y = richness, colour = sample_type)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Richness_Location_bio_sed_water.tiff")

plot.even <- ggplot(alpha, aes(x = sample_type, y = pielou, colour = sample_type)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Pielou's_Evenness_Location_bio_water.tiff")

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_sample_type_all.tiff")
```

## Plot Diversity Boxplots: Sample Depth
```{r}
plot.shan <- ggplot(alpha, aes(x = depth, y = shannon, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Location_bio_sed_water.tiff")

plot.rich <-ggplot(alpha, aes(x = depth, y = richness, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Richness_Location_bio_sed_water.tiff")

plot.even <- ggplot(alpha, aes(x = depth, y = pielou, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Pielou's_Evenness_Location_bio_water.tiff")

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_depth_all.tiff")
```

## Plot Diversity Boxplots: Metal Type
```{r}
plot.shan <- ggplot(alpha, aes(x = metal_type, y = shannon, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
#ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Location_bio_sed_water.tiff")

plot.rich <-ggplot(alpha, aes(x = metal_type, y = richness, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
#ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Richness_Location_bio_sed_water.tiff")

plot.even <- ggplot(alpha, aes(x = metal_type, y = pielou, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
#ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Pielou's_Evenness_Location_bio_water.tiff")

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_metal_type_all.tiff")
```
-----------------------------------------------------------------------
NEED TO MAKE EXCEL SHEETS!!
CODE NOT RUN YET 10/21/2024

## Diversity Metrics
```{r}
# Making OTU Table
df1_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_1_biof.csv")
df1_biof

df2_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_2_biof.csv")
df2_biof

df3_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_3_biof.csv")
df3_biof

df_otu_biof <- list(df1_biof, df2_biof, df3_biof)
df_otu_biof

write.table(df_otu_biof,"/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_biof_saipan.csv", sep=",", col.names=NA)
```

CODE NOT RUN YET 10/21/2024
```{r}
df_otu_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/df_otu_biof_saipan.csv")
df_otu_biof <- subset(df_otu, select = -c(X))
df_otu_biof

otu_count_biof <- df_otu_biof %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_biof

write.table(otu_count_biof, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_biof_saipan.csv", sep=",", quote=F, col.names=NA)
```

CODE NOT RUN YET 10/21/2024
```{r}
otu_count_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_biof_saipan.csv")

otu_count_biof %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

CODE NOT RUN YET 10/21/2024
```{r}
otu_count_biof <- otu_count_biof %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count_biof, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/otu_count_biof_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

CODE NOT RUN YET 10/21/2024
```{r}
diversity_metrics_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/diversity_metrics_biof.csv")
diversity_metrics_biof

write.table(diversity_metrics_biof, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/saipan_biof_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

CODE NOT RUN YET 10/21/2024
```{r}
diversity_metrics_biof %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, color= location)) +
  geom_boxplot() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/RStudio_Saipain_Data/results/alpha_diversity_metrics_biof.tiff", width = 10, height = 10)

#Each point represents a sample, (n) Sum of Count, (X) Total number of sequences for each sample & (Y) Value of diversity metric
```

CODE NOT RUN YET 10/21/2024
```{r}
otu_count_biof %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count),
            evenness = shannon/log(richness),
            n=sum(count)) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, color= sample_id,fill=sample_id)) +
  geom_point() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/alpha_diversity_metrics_biof.tiff", width = 30, height = 30)

#Each point represents a sample, (Y) Value of metric & (X) Total number of sequences for each sample
```

CODE NOT RUN YET 10/21/2024
```{r}
# Shannon, Richness and Evenness
data(otu.table.diver.biof)
H <- diversity(otu.table.diver.biof)
H

richness <- specnumber(otu.table.diver.biof)
richness

evenness <- H/log(richness)
evenness

metadata_biof <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/metadata_biof.csv")
metadata_biof

alpha_biof <- cbind(shannon = H, richness = richness, pielou = evenness, metadata_all)
write.csv(alpha, "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/diversity_indices_biof.csv")
head(alpha_biof)
```

CODE NOT RUN YET 10/21/2024
## Plot Diversity Boxplots: Biofilm Only - Sample Location
```{r}
plot.shan <- ggplot(alpha_biof, aes(x = location, y = shannon, colour = location)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan

plot.rich <-ggplot(alpha_biof, aes(x = location, y = richness, colour = location)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich

plot.even <- ggplot(alpha_biof, aes(x = location, y = pielou, colour = location)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_biof.tiff")
```

CODE NOT RUN YET 10/21/2024
## Plot Diversity Boxplots: Biofilm Only - Metal Type
```{r}
plot.shan <- ggplot(alpha_biof, aes(x = metal_type, y = shannon, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan

plot.rich <-ggplot(alpha_biof, aes(x = metal_type, y = richness, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich

plot.even <- ggplot(alpha_biof, aes(x = metal_type, y = pielou, colour = metal_type)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_biof_metaL_type.tiff")
```

CODE NOT RUN YET 10/21/2024
## Plot Diversity Boxplots: Biofilm Only - Depth
```{r}
plot.shan <- ggplot(alpha_biof, aes(x = depth, y = shannon, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan

plot.rich <-ggplot(alpha_biof, aes(x = depth, y = richness, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich

plot.even <- ggplot(alpha_biof, aes(x = depth, y = pielou, colour = depth)) +
geom_boxplot(size = 1) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/Shannon_Richness_Eveness_biof_depth.tiff")
```

CODE NOT RUN YET 10/21/2024
#SIMPER Analysis: 
NEED TO RUN THIS WHOLE CHUNK & GET WORKING!!!
```{r}
# Biofilm, Sediment and Water Samples
otu_table_all <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/asv_otu_saipan.csv", header=T, row.names=1, check.names=FALSE)

metadata <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/metadata_saipan.csv")

## Transpose the data to have sample names on rows
otu.table.diver.all <- t(otu_table_all)
otu.table.diver.all <- as.data.frame(otu.table.diver.all)
head(otu.table.diver.all)

otu.table.diver.mdf.all <- as.matrix.data.frame(otu.table.diver.all)
rownames(otu.table.diver.mdf.all) <- metadata$location

otu.table.diver.bray.all <- vegdist(otu.table.diver.mdf.all, method="bray")
otu.table.diver.bray.all

# Simper
simper_all <- simper(otu.table.diver.all, metadata$location, permutations=999)
options(max.print=10)

#summary(simper_all)
dput(simper_all, file = "simp_all.txt")
sim_all <- dget("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/simp_all.txt")
summary(sim_all)

simper <- simper(m_com, group=pc$location, permutations = 999)
simper
dput(simper, file = "/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/simp.txt")
```

ASK APRIL ABOUT ABUNDNACE CURVES!!
