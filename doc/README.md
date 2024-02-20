---
title: "Saipan_dada2_Pipeline"
author: "Maggie Shostak"
date: "`r Sys.Date()`"
output: html_document
---
# Load all necessary packages: 
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("microbiome")
#install.packages("vegan")
#library("vegan")

#install.packages("devtools")
#library("devtools")
```

```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
#library(decontam)
library(rmarkdown)
library(BiocGenerics)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
library(DECIPHER)
library(phangorn)
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
path <- "/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/fastq"
#list.files(path)
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
#FWD <- "GTGYCAGCMGCCGCGGTAA"
#REV <- "CCGYCAATTYMTTTRAGTTT"
#trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
```

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
#head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

```{r}
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
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "track_sequences.csv", sep=",", quote=F, col.names=NA)
```

## Taxonomic Classification & ASV File Formatting
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/fastq/silva_nr99_v138.1_train_set.fa.gz")
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
#write.table(asv_fasta, "asv_fasta_saipan.fa")
#write.table(asv_otu, "asv_otu_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(asv_tax, "asv_tax_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(otu_tax_table, "otu_tax_table_saipan.csv", sep=",", quote=F, col.names=NA)
```

# Sequence Alignment
```{r}
#library(DECIPHER)
#sequences <- getSequences(seqtab.nochim)
#names(sequences) <- sequences
#alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
```

```{r}
#library(phangorn)
#phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
#dist_matrix <- dist.ml(phang.align)
#treeNJ <- pml(treeNJ, data= dist_matrix)
#fit = pml(treeNJ, data = phang.align)
#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement= "stochastic", control=pml.control(trace=0))
```

```{r}
#map <- import_qiime_sample_data
```

# Formating Files for Future Analysis
```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata.csv")
metadata

otu_counts <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_tax_saipan.csv")
taxonomy

asv_tax <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_tax_saipan.csv")
asv_tax
```

# Phyloseq
```{r}
library(phyloseq)
library(Biostrings)
```

```{r}
# Read data into R
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan.csv")
otu_tab
tax_tab <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_tax_saipan.csv")
tax_tab
samples_df <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata.csv")
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
```

```{r}
# Transform OTU & Tax table into matrices
otu_tab <- as.matrix(otu_tab)
tax_tab <- as.matrix(tax_tab)
```

```{r}
#Transform into Phyloseq Objects
ASV = otu_table(otu_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
sample_id = sample_data(samples_df)
  
ps <- phyloseq(ASV, TAX, sample_id)
ps


#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 34762 taxa and 141 samples ]
#sample_data() Sample Data:       [ 141 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 34762 taxa by 6 taxonomic ranks ]
```

```{r}
# Visualize Data
sample_names(ps)
rank_names(ps)
sample_variables(ps)
```

```{r}
# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)
```

```{r}
# Bar graphs based on division
plot_bar(ps, fill = "Phylum")
```

```{r}
plot_bar(ps, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
```

```{r}
# Heatmaps
plot_heatmap(ps, method = "NMDS", distance = "bray")
```

```{r}
ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.20) > 0, TRUE)
ps_abund


#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6 taxa and 141 samples ]
#sample_data() Sample Data:       [ 141 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 6 taxa by 6 taxonomic ranks ]
```

```{r}
otu_table(ps)[1:8, 1:5]
```

```{r}
plot_heatmap(ps_abund, method = "NMDS", distance = "bray")
```

```{r}
# Alpha Diversity: Plot Chao1 richness estimator and Shannon diversity estimator
plot_richness(ps, measures=c("Chao1", "Shannon"))
```

```{r}
# Regroup together from the same location
plot_richness(ps, measures=c("Chao1", "Shannon"), x="location", color="location") +
  geom_boxplot()
#ggsave("Chao1_Shannon_Saipan_All_Location.tiff", width=15, height=10)
```

```{r}
## Just Shannon
plot_richness(ps, measures=c("Shannon"), x="location", color="location") +
  geom_boxplot()
ggsave("Chao1_Shannon_Saipan_All_Location.tiff", width=15, height=10)
```

```{r}
# Ordination
ps.ord <- ordinate(ps, "NMDS", "bray")
plot_ordination(ps, ps.ord, type="taxa", color="Phylum", shape= "Location", 
                  title="Saipan ASVs")
```

```{r}
plot_ordination(ps, ps.ord, type="taxa", color="Phylum", 
                  title="Sapain ASVs", label="Phylum") + 
                  facet_wrap(~Phylum, 3)
ggsave("Ordination_Phylum_Saipan_All.tiff", width=25, height=10)
```

```{r}
plot_ordination(ps, ps.ord, type="samples", color="Location", 
                  shape="depth", title="Samples") + geom_point(size=3)
```

```{r}
plot_ordination(ps, ps.ord, type="split", color="Phylum", 
                  shape="Location", title="Biplot", label = "station") +  
  geom_point(size=3)
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

#write.table(otu_rel_abund, "otu_rel_abund_saipan.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/otu_rel_abund_saipan.csv")
#otu_rel_abund
```

## Other Metadata Tables
```{r}
metadata_biof_sed_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_sed_water.csv")
#metadata_biof_sed_water

metadata_depth_biofilm <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_biofilm.csv")
metadata_depth_biofilm

metadata_depth_sediment <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_sediment.csv")
#metadata_depth_sediment

metadata_depth_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_water.csv")
#metadata_depth_water

metadata_location_biof_sed_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biof_sed_water.csv")
#metadata_location_biof_sed_water

metadata_location_biofilm <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biofilm.csv")
#metadata_location_biofilm

metadata_location_sediment <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_sediment.csv")
#metadata_location_sediment

#metadata_location_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_water.csv")
#metadata_location_water

metadata_metal_type <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_metal_type.csv")
#metadata_metal_type
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
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("phylum_stacked_barchart_saipan_all.tiff", width=20, height=10)

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
#ggsave("class_stacked_barchart_saipan_all.tiff", width=25, height=10)

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
ggsave("order_stacked_barchart_saipan_all.tiff", width=55, height=10, limitsize = FALSE)

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
ggsave("family_stacked_barchart_saipan_all.tiff", width=50, height=10, limitsize = FALSE)
```

### Biofilm Samples Only
```{r}
metadata_location_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biofilm.csv")
#metadata_location_biof_simper

otu_counts_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan_biof_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_biof_simper

taxonomy_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_saipan_biof_simper.csv")
#taxonomy_biof_simper
```

```{r}
otu_rel_abund_biof_simper <- inner_join(metadata_location_biof_simper, otu_counts_biof_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_biof_simper
#write.table(otu_rel_abund_biof_simper, "otu_rel_abund_saipan_biof_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## Phylum
otu_rel_abund_biof_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_biof_simper.tiff", width=25, height=15)

## Class
otu_rel_abund_biof_simper %>%
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
#ggsave("class_stacked_barchart_biof_simper.tiff", width=25, height=15)

## Order
otu_rel_abund_biof_simper %>%
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
ggsave("order_stacked_barchart_biof_simper.tiff", width=55, height=10, limitsize = FALSE)

## Family
otu_rel_abund_biof_simper %>%
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
ggsave("Family_stacked_barchart_biof_simper.tiff", width=55, height=10, limitsize = FALSE)
```

### Coronado Biofilm vs Deep Coronado Biofilm
```{r}
# Data Manipulation
metadata_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_coronados_comparisons.csv")
#metadata_biof_coronados_simper

otu_counts_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_coronados_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_biof_coronados_simper

taxonomy_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_coronados_simper.csv")
#taxonomy_biof_coronados_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_biof_coronados_simper <- inner_join(metadata_biof_coronados_simper, otu_counts_biof_coronados_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_coronados_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_biof_coronados_simper
#write.table(otu_rel_abund_biof_coronados_simper, "otu_rel_abund_biof_coronados_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_biof_coronados_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_biof_coronados_simper.tiff", width=20, height=15)

### Class
otu_rel_abund_biof_coronados_simper %>%
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
#ggsave("class_stacked_barchart_biof_coronados_simper.tiff", width=17, height=15)

# Order
otu_rel_abund_biof_coronados_simper %>%
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
ggsave("order_stacked_barchart_biof_coronados_simper.tiff", width=55, height=10, limitsize = FALSE)

#Family
otu_rel_abund_biof_coronados_simper %>%
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
ggsave("family_stacked_barchart_biof_coronados_simper.tiff", width=55, height=10, limitsize = FALSE)
```

### Tank 1 vs Tank 3
```{r}
# Data Manipulation
metadata_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_T1T3.csv")
#metadata_biof_T1T3_simper

otu_counts_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_T1T3_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_biof_T1T3_simper

taxonomy_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_T1T3_simper.csv")
#taxonomy_biof_T1T3_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_biof_T1T3_simper <- inner_join(metadata_biof_T1T3_simper, otu_counts_biof_T1T3_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_T1T3_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_biof_T1T3_simper
#write.table(otu_rel_abund_biof_T1T3_simper, "otu_rel_abund_biof_T1T3_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_biof_T1T3_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_biof_T1T3_simper.tiff", width=13, height=15)

### Class
otu_rel_abund_biof_T1T3_simper %>%
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
#ggsave("class_stacked_barchart_biof_T1T3_simper.tiff", width=13, height=15)

### Order
otu_rel_abund_biof_T1T3_simper %>%
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
ggsave("order_stacked_barchart_biof_T1T3_simper.tiff", width=55, height=10, limitsize = FALSE)
### Family
otu_rel_abund_biof_T1T3_simper %>%
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
ggsave("family_stacked_barchart_biof_T1T3_simper.tiff", width=55, height=10, limitsize = FALSE)
```

### Emily vs Coronado
```{r}
metadata_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_EmCoro.csv")
#metadata_biof_EmCoro_simper

otu_counts_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_EmCoro_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_biof_EmCoro_simper

taxonomy_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_EmCoro_simper.csv")
#taxonomy_biof_EmCoro_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_biof_EmCoro_simper <- inner_join(metadata_biof_EmCoro_simper, otu_counts_biof_EmCoro_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_EmCoro_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_biof_EmCoro_simper
#write.table(otu_rel_abund_biof_EmCoro_simper, "otu_rel_abund_biof_EmCoro_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_biof_EmCoro_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_biof_EmCoro_simper.tiff", width=13, height=15)

### Class
otu_rel_abund_biof_EmCoro_simper %>%
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
#ggsave("class_stacked_barchart_biof_EmCoro_simper.tiff", width=13, height=15)

### Order
otu_rel_abund_biof_EmCoro_simper %>%
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
ggsave("order_stacked_barchart_biof_EmCoro_simper.tiff", width=55, height=10, limitsize = FALSE)


### Family
otu_rel_abund_biof_EmCoro_simper %>%
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
ggsave("family_stacked_barchart_biof_EmCoro_simper.tiff", width=55, height=10, limitsize = FALSE)
```

### Diahatsu 1 vs Shoan Maru
```{r}
metadata_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_D1SM_simper.csv")
#metadata_biof_D1SM_simper

otu_counts_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_D1SM_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_biof_D1SM_simper

taxonomy_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_D1SM_simper.csv")
#taxonomy_biof_D1SM_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_biof_D1SM_simper <- inner_join(metadata_biof_D1SM_simper, otu_counts_biof_D1SM_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_T1T3_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_biof_D1SM_simper
#write.table(otu_rel_abund_biof_D1SM_simper, "otu_rel_abund_biof_D1SM_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_biof_D1SM_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_biof_D1SM_simper.tiff", width=13, height=15)

### Class
otu_rel_abund_biof_D1SM_simper %>%
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
#ggsave("class_stacked_barchart_biof_D1SM_simper.tiff", width=13, height=15)

### Order
otu_rel_abund_biof_D1SM_simper %>%
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
ggsave("order_stacked_barchart_biof_D1SM_simper.tiff", width=55, height=30, limitsize = FALSE)

## Family
otu_rel_abund_biof_D1SM_simper %>%
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
ggsave("family_stacked_barchart_biof_D1SM_simper.tiff", width=55, height=20, limitsize = FALSE)
```

## Aluminun vs Steel vs Iron
```{r}
metadata_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_metal_type.csv")
#metadata_metal_type_simper

otu_counts_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_metal_type_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_metal_type_simper

taxonomy_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_metal_type_simper.csv")
#taxonomy_metal_type_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_metal_type_simper <- inner_join(metadata_metal_type_simper, otu_counts_metal_type_simper, by="sample_id") %>%
  inner_join(., taxonomy_biof_T1T3_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_metal_type_simper
write.table(otu_rel_abund_metal_type_simper, "otu_rel_abund_metal_type_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_metal_type_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("phylum_stacked_barchart_biof_metal_type_simper.tiff", width=13, height=15)

### Class
otu_rel_abund_metal_type_simper %>%
  filter(level=="Class") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("class_stacked_barchart_biof_metal_type_simper.tiff", width=13, height=15)

## Order
otu_rel_abund_metal_type_simper %>%
  filter(level=="Order") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("order_stacked_barchart_biof_metal_type_simper.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_metal_type_simper %>%
  filter(level=="Family") %>%
  group_by(sample_id, metal_type, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(metal_type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=metal_type, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=metal_type, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("family_stacked_barchart_biof_metal_type_simper.tiff", width=55, height=20, limitsize = FALSE)
```

## Sediment Samples Only
```{r}
metadata_location_sed <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_sed.csv")
#metadata_location_sed

otu_counts_sed_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_location_sed_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_sed_simper

taxonomy_sed_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_location_sed_simper.csv")
#taxonomy_sed_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_sed_simper <- inner_join(metadata_location_sed, otu_counts_sed_simper, by="sample_id") %>%
  inner_join(., taxonomy_sed_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_sed_simper
#write.table(otu_rel_abund_sed_simper, "otu_rel_abund_sed_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_sed_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_sed_simper.tiff", width=25, height=15)

### Class
otu_rel_abund_sed_simper %>%
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
#ggsave("class_stacked_barchart_sed_simper.tiff", width=25, height=15)

## Order
otu_rel_abund_sed_simper %>%
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
ggsave("order_stacked_barchart_sed_simper.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_sed_simper %>%
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
ggsave("family_stacked_barchart_sed_simper.tiff", width=55, height=20, limitsize = FALSE)
```

### Coronado Sediment vs Deep Coronado Sediment
```{r}
metadata_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_sed_coronados_simper.csv")
#metadata_sed_coronados_simper

otu_counts_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_sed_coronados_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_sed_coronados_simper

taxonomy_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_sed_coronados_simper.csv")
#taxonomy_sed_coronados_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_sed_coronados_simper <- inner_join(metadata_sed_coronados_simper, otu_counts_sed_coronados_simper, by="sample_id") %>%
  inner_join(., taxonomy_sed_coronados_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_sed_coronados_simper
#write.table(otu_rel_abund_sed_coronados_simper, "otu_rel_abund_sed_coronados_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_sed_coronados_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_sed_coronados_simper.tiff", width=13, height=15)

### Class
otu_rel_abund_sed_coronados_simper %>%
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
#ggsave("class_stacked_barchart_sed_coronados_simper.tiff", width=13, height=15)

## Order
otu_rel_abund_sed_coronados_simper %>%
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
ggsave("order_stacked_barchart_sed_coronados_simper.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_sed_coronados_simper %>%
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
ggsave("family_stacked_barchart_sed_coronados_simper.tiff", width=55, height=20, limitsize = FALSE)
```

## Water Samples Only
```{r}
metadata_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_water.csv")
#metadata_water

otu_counts_water_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_location_water_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_water_simper

taxonomy_water_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_location_water_simper.csv")
#taxonomy_water_simper
```

```{r}
## Relative Abundance Generation
otu_rel_abund_water_simper <- inner_join(metadata_water, otu_counts_water_simper, by="sample_id") %>%
  inner_join(., taxonomy_water_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
#otu_rel_abund_water_simper
#write.table(otu_rel_abund_water_simper, "otu_rel_abund_water_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## GG Plots
### Phylum
otu_rel_abund_water_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("phylum_stacked_barchart_water_simper.tiff", width=20, height=15)

### Class
otu_rel_abund_water_simper %>%
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
#ggsave("class_stacked_barchart_water_simper.tiff", width=20, height=15)

## Order
otu_rel_abund_water_simper %>%
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
ggsave("order_stacked_barchart_water_simper.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_water_simper %>%
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
ggsave("family_stacked_barchart_water_simper.tiff", width=55, height=20, limitsize = FALSE)
```

## All Samples Depth Comparison
```{r}
metadata_depth <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_biofilm.csv")
#metadata_depth

otu_counts_depth_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_depth_biof.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
#otu_counts_depth_simper

taxonomy_depth_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_tax_depth_biof.csv")
#taxonomy_depth_simper
```

```{r}
otu_rel_abund_depth_simper <- inner_join(metadata_depth, otu_counts_depth_simper, by="sample_id") %>%
  inner_join(., taxonomy_depth_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_depth_simper
write.table(otu_rel_abund_depth_simper, "otu_rel_abund_saipan_depth_simper.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## Phylum
otu_rel_abund_depth_simper %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("phylum_stacked_barchart_depth_simper.tiff", width=27, height=15)

## Class
otu_rel_abund_depth_simper %>%
  filter(level=="Class") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
#ggsave("class_stacked_barchart_depth_simper.tiff", width=13, height=15)

## Order
otu_rel_abund_depth_simper %>%
  filter(level=="Class") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("order_stacked_barchart_depth_simper.tiff", width=55, height=20, limitsize = FALSE)

## Family
otu_rel_abund_depth_simper %>%
  filter(level=="Family") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("family_stacked_barchart_depth_simper.tiff", width=55, height=20, limitsize = FALSE)
```

# Non-multidimensional Scaling (NMDS Plots)
```{r}
df_meta <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_otu_categories.csv")
df_meta

nmds_asv_otu_all <- inner_join(df_meta, df_otu, by="sample_id")
nmds_asv_otu_all

write.table(nmds_asv_otu_all, "nmds_asv_otu_all.csv", sep=",", quote=F, col.names=NA)
```

```{r}
pc <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_asv_otu_all.csv")
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
#ggsave("NMDS_all_samples_location.tiff", width = 10, height = 10)
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
#ggsave("NMDS_all_samples_sample_type.tiff", width = 10, height = 10)
```

## All Samples: By Depth
```{r}
xx2 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = depth))+
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
 labs(x = "NMDS1", colour = "Depth", y = "NMDS2")
xx2
#ggsave("NMDS_all_samples_depth.tiff", width = 10, height = 10)
```

```{r}
df_meta_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_biof_only_categories.csv")
#df_meta_biof

df_otu_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_biof_")

nmds_asv_otu_biof <- inner_join(df_meta_biof, df_otu_biof, by="sample_id")
nmds_asv_otu_biof

#write.table(nmds_asv_otu_biof, "nmds_asv_otu_biof_only.csv", sep=",", quote=F, col.names=NA)
```

```{r}
pc1 <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_asv_otu_biof_only.csv")
pc1

#make community matrix: extract columns with ASV information
com1 <- pc1[,7:ncol(pc1)]
com1

#turn ASV information into a matrix
m_com1 <- as.matrix(com1)
```

```{r}
#Run NMDS using Bray-Curtis distance
set.seed(123)
nmds1 <- metaMDS(m_com1, distance="bray") #stress = 0.1910382 
nmds1
plot(nmds1)

#access the specific points data of the NMDS plot & scores
str(nmds1)
nmds1$points
scores(nmds1)

#extract NMDS scores
data.scores1 = as.data.frame(scores(nmds1)$sites)

#add columns to data frame
data.scores1$sample_id = pc1$sample_id
data.scores1$location = pc1$location
data.scores1$depth = pc1$depth
data.scores1$metal_type = pc1$metal_type

head(data.scores1)
```

## Biofilm Only: By Location
```{r}
xx3 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2)) +
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
 labs(x = "NMDS1", colour = "Location", y = "NMDS2")
xx3
#ggsave("NMDS_biof_only_location.tiff", width = 10, height = 10)
```

## Biofilm Only: By Metal Type
```{r}
xx4 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = metal_type, color=location))+
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
 labs(x = "NMDS1", colour = "Metal Type", y = "NMDS2")
xx4
ggsave("NMDS_biof_only_metal_type.tiff", width = 10, height = 10)
```

## Biofilm Only: By Depth
```{r}
xx5 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = depth))+
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
 labs(x = "NMDS1", colour = "Depth", y = "NMDS2")
xx5
#ggsave("NMDS_biof_only_depth.tiff", width = 10, height = 10)
```

## Biofilm Only: Coronado vs Deep Coronado
```{r}

```

## Biofilm Only: Emily vs Coronado vs Deep Coronado
```{r}

```

## Biofilm Only: Tank 1 vs Tank 3
```{r}

```

## Biofilm Only: Diahatsu 1 vs Shoan Maru
```{r}

```

## Deep Coronado Biof/Sed/Water vs Coronado Biof/Sed/Water
```{r}

```

## All Sediment: By Location
```{r}

```

## All Water: By Location
```{r}

```

# Statistical Analyses
## Anosim
```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/nmds_asv_otu_all.csv")
pc_ano
```

```{r}
# All Samples: Biof vs Sed vs Water by Location
com_ano = pc_ano[,6:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

# ANOSIM statistic R: 0.1327 
      # Significance: 1e-04
```

```{r}
# All Samples: Biof vs Sed vs Water by Sample Type
com_ano1 = pc_ano[,6:ncol(pc_ano)]
m_com_ano1 = as.matrix(com_ano1)
ano_all1 = anosim(m_com_ano1, pc_ano$sample_type, distance = "bray", permutations = 9999)
ano_all1

# ANOSIM statistic R: 0.8735 
      # Significance: 1e-04 
```

```{r}
# All Samples: Biof vs Sed vs Water by Depth
com_ano2 = pc_ano[,6:ncol(pc_ano)]
m_com_ano2 = as.matrix(com_ano2)
ano_all2 = anosim(m_com_ano2, pc_ano$depth, distance = "bray", permutations = 9999)
ano_all2

# ANOSIM statistic R: 0.221 
      # Significance: 1e-04 
```

```{r}
# Biofilm Samples Only by Location

```

# Ecological Distance Matrices
```{r}
#Diversity index value generating
otu_table<-read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan.csv",header=T,row.names=1,check.names=FALSE)
otu_table <- as.data.frame.matrix(otu_table)
otu_table

#Transpose the data to have sample names on rows
otu_table<-t(otu_table)
data(otu_table)
H<-diversity(otu_table)
simp<-diversity(otu_table, "simpson")
invsimp<-diversity(otu_table, "inv")

## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(otu_table, 2) - 1

## Fisher alpha
alpha <- fisher.alpha(otu_table)
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(otu_table)
J <- H/log(S)

## Plot all
pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="blue")
alpha
write.table(J, "/Users/maggieshostak/Desktop/pielou_evenness.txt", sep="\t")
```

```{r}
library(tidyverse)

df1 <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df1.csv")
#df1
df2 <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df2.csv")
#df2
df3 <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df3.csv")
#df3
df_otu <- list(df1, df2, df3)
df_otu

#write.table(df_otu,"df_otu.csv", sep=",", col.names=NA)
```

```{r}
df_otu <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df_otu.csv")
df_otu <- subset(df_otu, select = -c(X))
df_otu

otu_count_all <- df_otu %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_all

#write.table(otu_count, "otu_count.csv", sep=",", quote=F, col.names=NA)

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
otu_count_all %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            simpson = simpson(count),
            invsimp = 1/simpson,
            evenness = shannon/log(richness),
            n=sum(count))
```

## Plot Alpha Diversity Functions
```{r}
otu_count_all %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            simpson = simpson(count),
            invsimp = 1/simpson,
            evenness = shannon/log(richness),
            n=sum(count)) %>%
  pivot_longer(cols=c(richness, shannon, invsimp, simpson, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("alpha_diversity_metrics_all.tiff", width = 10, height = 10)

#Each point represents a sample, (Y) Value of metric & (X) Total number of sequences for each sample
```

# Ecological Distance Matrices: Biofilm Only
```{r}
## Biofilm Only
otu_table_bio <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan_biof_simper.csv", header=T, row.names=1, check.names=FALSE)
otu_table_bio

## Transpose the data to have sample names on rows
otu.table.diver.bio <- t(otu_table_bio)
otu.table.diver.bio <- as.data.frame(otu.table.diver.bio)
head(otu.table.diver.bio)
```

```{r}
df1_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df1_biof_only.csv")
df1_biof

df2_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df2_biof_only.csv")
df2_biof

df3_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df3_biof_only.csv")
df3_biof

df_otu_biof <- list(df1_biof, df2_biof, df3_biof)
df_otu_biof

#write.table(df_otu_biof,"df_otu_biof.csv", sep=",", col.names=NA)
```

```{r}
df_otu_biof <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/df_otu_biof.csv")
df_otu_biof <- subset(df_otu_biof, select = -c(X))
df_otu_biof

otu_count_biof <- df_otu_biof %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_biof

#write.table(otu_count_biof, "otu_count_biof.csv", sep=",", quote=F, col.names=NA)
```

```{r}
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

```{r}
otu_count_biof %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            simpson = simpson(count),
            invsimp = 1/simpson,
            evenness = shannon/log(richness),
            n=sum(count))
```

## Plot Alpha Diversity Functions
```{r}
otu_count_biof %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            simpson = simpson(count),
            invsimp = 1/simpson,
            evenness = shannon/log(richness),
            n=sum(count)) %>%
  pivot_longer(cols=c(richness, shannon, invsimp, simpson, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("alpha_diversity_metrics_biof.tiff", width = 10, height = 10)

#Each point represents a sample, (Y) Value of metric & (X) Total number of sequences for each sample
```

# SIMPER Analysis
```{r}
# Biofilm, Sediment and Water Samples
otu_table_all <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan.csv", header=T, row.names=1, check.names=FALSE)

## Transpose the data to have sample names on rows
otu.table.diver.all <- t(otu_table_all)
otu.table.diver.all <- as.data.frame(otu.table.diver.all)
head(otu.table.diver.all)

otu.table.diver.mdf.all <- as.matrix.data.frame(otu.table.diver.all)
rownames(otu.table.diver.mdf.all) <- metadata$location

otu.table.diver.bray.all <- vegdist(otu.table.diver.mdf.all, method="bray")
otu.table.diver.bray.all

# Simper
simper_all <- simper_all(otu.table.diver.all, metadata$location, permutations=999)
options(max.print=999999)
#summary(simper_all)
dput(simper_all, file = "simp_all.txt")
sim_all <- dget("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/simp_all.txt")
#summary(sim_all)

simper <- simper(m_com, group=pc$location, permutations = 999)
simper
dput(simper, file = "simp.txt")
```
