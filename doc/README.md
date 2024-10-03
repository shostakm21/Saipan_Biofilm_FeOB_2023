---
Title: "Saipan Dada2 Pipeline"
Author: "Maggie Shostak"
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
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
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


