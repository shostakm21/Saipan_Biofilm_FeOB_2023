Saipan Data Analysis Pipeline

# Load all Necessary Packages:
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
list.files(path)
```

```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
```

```{r}
plotQualityProfile(fnFs[1:2])
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

```
                                  reads.in reads.out
Shos1_S31_L001_R1_001.fastq.gz       43232     33472
Shos10_S44_L001_R1_001.fastq.gz      28882     23937
Shos100_S167_L001_R1_001.fastq.gz    54542     45358
Shos101_S177_L001_R1_001.fastq.gz    29016     23487
Shos102_S187_L001_R1_001.fastq.gz    41947     33891
Shos103_S196_L001_R1_001.fastq.gz    56944     47914
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

```
XXX
```

```{r}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
```

```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
head(mergers[[1]])
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
```

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Dissertation/RStudio_Saipan/Saipan/data/fastq/silva_nr99_v138.1_train_set.fa")
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
write(asv_fasta, "asv_fasta_saipan.fa")
write.table(asv_otu, "asv_otu_saipan.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_saipan.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table_saipan.csv", sep=",", quote=F, col.names=NA)
```

# Formating Files for Future Analysis
```{r}
metadata_all <- read.csv("****")
metadata_all

otu_counts <- read.csv("????????asv_otu_saipan.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv(">>>>>>>asv_tax_saipan.csv")
taxonomy

asv_tax <- read.csv("??????asv_tax_saipan.csv")
asv_tax
```

## OTU Relative Abundance Generation
```{r}
otu_rel_abund <- inner_join(metadata_all, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "otu_rel_abund_saipan.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_rel_abund <- read.csv("???????otu_rel_abund_saipan.csv")
otu_rel_abund
```

## Metadata Excel Data Tables





