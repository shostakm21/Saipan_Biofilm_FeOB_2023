

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
