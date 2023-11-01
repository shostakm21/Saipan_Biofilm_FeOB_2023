Saipan Data Analysis Pipeline

# Load all necessary packages: 
```{r}
library(dada2)
library(phyloseq)
library(ggplot2)
library(ape)
library(microbiome)
library(tidyverse)
library(kableExtra)
library(phangorn)
library(DECIPHER)
library(reshape2)
library(treeio)
library(ggtree)
library(ggstance)
library(scales)
library(dplyr)
library(ggpattern)
library(vegan)
library(MASS)
library(ecodist)
library(scatterplot3d)
library(fso)
library(dbplyr)
library(svglite)
library(ggstatsplot)
library(grid)
library(labdsv)
library(indicspecies)
library(viridis)
library(ggplot2)
library('cowplot')
library(permute)
library(lattice)
library(breakaway)
library(dtplyr)
library(Rcpp)
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

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegan")
library("vegan")

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 
#change the ref argument to get other versions
```

# DADA2 Pipeline
```{r}
path <- "/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/fastq"
#list.files(path)
```

```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#list(sample.names)
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
```

```{r}
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
#head(track)
```

## Taxonomic Classification & ASV File Formatting
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/fastq/silva_nr99_v138.1_train_set.fa.gz")
```

```{r}
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
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
#write(asv_fasta, "asv_fasta_saipan.fa")
#write.table(asv_otu, "asv_otu_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(asv_tax, "asv_tax_saipan.csv", sep=",", quote=F, col.names=NA)
#write.table(otu_tax_table, "otu_tax_table_saipan.csv", sep=",", quote=F, col.names=NA)
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
otu_rel_abund
```

## Other Metadata Tables
```{r}
metadata_biof_sed_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata.csv")
metadata_biof_sed_water
```

```{r}
metadata_depth_biofilm <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_biofilm.csv")
metadata_depth_biofilm
```

```{r}
metadata_depth_sediment <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_sediment.csv")
metadata_depth_sediment
```

```{r}
metadata_depth_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth_water.csv")
metadata_depth_water
```

```{r}
metadata_location_biof_sed_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biof_sed_water.csv")
metadata_location_biof_sed_water
```

```{r}
metadata_location_biofilm <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biofilm.csv")
metadata_location_biofilm
```

```{r}
metadata_location_sediment <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_sediment.csv")
metadata_location_sediment
```

```{r}
metadata_location_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_water.csv")
metadata_location_water
```

```{r}
metadata_metal_type <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_metal_type.csv")
metadata_metal_type
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
ggsave("phylum_stacked_barchart_saipan_all.tiff", width=20, height=7)

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
ggsave("class_stacked_barchart_saipan_all.tiff", width=25, height=10)
```

### Biofilm Samples Only
```{r}
metadata_location_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_biofilm.csv")
metadata_location_biof_simper
```

```{r}
otu_counts_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan_biof_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_biof_simper
```

```{r}
taxonomy_biof_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_saipan_biof_simper.csv")
taxonomy_biof_simper
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
otu_rel_abund_biof_simper
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
#ggsave("phylum_stacked_barchart_biof_simper.tiff", width=13, height=15)

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
#ggsave("class_stacked_barchart_biof_simper.tiff", width=13, height=15)
```

### Coronado Biofilm vs Deep Coronado Biofilm
```{r}
# Data Manipulation
metadata_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_coronados_comparisons.csv")
metadata_biof_coronados_simper
```

```{r}
otu_counts_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_coronados_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_biof_coronados_simper
```

```{r}
taxonomy_biof_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_coronados_simper.csv")
taxonomy_biof_coronados_simper
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
otu_rel_abund_biof_coronados_simper
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
#ggsave("phylum_stacked_barchart_biof_coronados_simper.tiff", width=13, height=15)

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
#ggsave("class_stacked_barchart_biof_coronados_simper.tiff", width=13, height=15)
```

### Tank 1 vs Tank 3
```{r}
# Data Manipulation
metadata_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_T1T3.csv")
metadata_biof_T1T3_simper
```

```{r}
otu_counts_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_T1T3_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_biof_T1T3_simper
```

```{r}
taxonomy_biof_T1T3_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_T1T3_simper.csv")
taxonomy_biof_T1T3_simper
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
otu_rel_abund_biof_T1T3_simper
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
```

### Emily vs Coronado
```{r}
metadata_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_EmCoro.csv")
metadata_biof_EmCoro_simper
```

```{r}
otu_counts_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_EmCoro_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_biof_EmCoro_simper
```

```{r}
taxonomy_biof_EmCoro_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_EmCoro_simper.csv")
taxonomy_biof_EmCoro_simper
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
otu_rel_abund_biof_EmCoro_simper
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
```

### Diahatsu 1 vs Shoan Maru
```{r}
metadata_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_biof_D1SM_simper.csv")
metadata_biof_D1SM_simper
```

```{r}
otu_counts_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_biof_D1SM_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_biof_D1SM_simper
```

```{r}
taxonomy_biof_D1SM_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_biof_D1SM_simper.csv")
taxonomy_biof_D1SM_simper
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
otu_rel_abund_biof_D1SM_simper
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
```

## Aluminun vs Steel vs Iron
```{r}
metadata_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_metal_type.csv")
metadata_metal_type_simper
```

```{r}
otu_counts_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_metal_type_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_metal_type_simper
```

```{r}
taxonomy_metal_type_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_metal_type_simper.csv")
taxonomy_metal_type_simper
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
```

## Sediment Samples Only
```{r}
metadata_location_sed <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_location_sed.csv")
metadata_location_sed
```

```{r}
otu_counts_sed_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_location_sed_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_sed_simper
```

```{r}
taxonomy_sed_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_location_sed_simper.csv")
taxonomy_sed_simper
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
otu_rel_abund_sed_simper
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
ggsave("phylum_stacked_barchart_sed_simper.tiff", width=13, height=15)

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
ggsave("class_stacked_barchart_sed_simper.tiff", width=13, height=15)
```

### Coronado Sediment vs Deep Coronado Sediment
```{r}
metadata_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_sed_coronados.csv")
metadata_sed_coronados_simper
#NOT RUN
```

```{r}
otu_counts_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_sed_coronados_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_sed_coronados_simper
#NOT RUN
```

```{r}
taxonomy_sed_coronados_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_sed_coronados_simper.csv")
taxonomy_sed_coronados_simper
#NOT RUN
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
otu_rel_abund_sed_coronados_simper
write.table(otu_rel_abund_sed_coronados_simper, "otu_rel_abund_sed_coronados_simper.csv", sep=",", quote=F, col.names=NA)
#NOT RUN
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
ggsave("phylum_stacked_barchart_sed_coronados_simper.tiff", width=13, height=15)
#NOT RUN

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
ggsave("class_stacked_barchart_sed_coronados_simper.tiff", width=13, height=15)
#NOT RUN
```

## Water Samples Only
```{r}
metadata_water <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_water.csv")
metadata_water
#NOT RUN
```

```{r}
otu_counts_water_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_water_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_water_simper
#NOT RUN
```

```{r}
taxonomy_water_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_water_simper.csv")
taxonomy_water_simper
#NOT RUN
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
otu_rel_abund_water_simper
write.table(otu_rel_abund_water_simper, "otu_rel_abund_water_simper.csv", sep=",", quote=F, col.names=NA)
#NOT RUN
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
ggsave("phylum_stacked_barchart_water_simper.tiff", width=13, height=15)
#NOT RUN

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
ggsave("class_stacked_barchart_water_simper.tiff", width=13, height=15)
#NOT RUN
```

## All Samples Depth Comparison
```{r}
metadata_depth <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/metadata_depth.csv")
metadata_depth
#NOT RUN
```

```{r}
otu_counts_depth_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_saipan_depth_simper.csv") %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count")
otu_counts_depth_simper
#NOT RUN
```

```{r}
taxonomy_depth_simper <- read.csv("/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/asv_otu_tax_saipan_depth_simper.csv")
taxonomy_depth_simper
#NOT RUN
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
#NOT RUN
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
ggsave("phylum_stacked_barchart_depth_simper.tiff", width=13, height=15)
#NOT RUN

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
ggsave("class_stacked_barchart_depth_simper.tiff", width=13, height=15)
#NOT RUN
```

# NMDS Plots
## All Samples: Bioiflm Grouped vs Sed vs Water
```{r}

```

## All Samples: Biofilms Separate vs Sed vs Water
```{r}

```

## Biofilm Only: By Location
```{r}

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

## All Samples: By Depth
```{r}

```

# Statistical Analyses
## Anosim
```{r}

```




## Diversity Index Value Generation
```{r}

```


## Ecological Distances Matrices & Rarefaction
```{r}

```

### Shannons H Diveristy
```{r}

```


### Spp Richness
```{r}

```


### Pielou Evenness
```{r}

```

# SIMPER Analysis
```{r}

```

