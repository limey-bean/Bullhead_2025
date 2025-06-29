---
title: "Differential Abundance Bullhead tissue types"
output: html_notebook
---

## This notebook attempts to answer the following question:

```         
Does the diversity/composition of bacterial communities differ between tissue types.
```

## Step 1. Load R packages

The packages should be installed with phyloseq or tidyverse.

```{r}
library(phyloseq)
library(tidyverse)
library(vegan)
library(DESeq2)
library(dplyr)
```

If there is a problem with the vegan or DESeq2 installation, do the following:

Install Vegan

```{r, eval=FALSE}

install.packages('vegan',
    repos = c('https://vegandevs.r-universe.dev','https://cloud.r-project.org'))

```

Install DESeq2

```{r, eval=FALSE}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

## Step 2. Get data into your machine.

```{r, eval=FALSE}
# an example for how to read wf-16S output file
file <- "~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/VIrus_hunting/otu_table_G.csv"
raw <- read.delim(file, sep=",")

# an example for how toread the metadatafile
meta <- "~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/mito_N1_N2/meta_data.csv"

samp_data <- read.delim(meta, sep=",")
```

## Step 3. Assign output directory path.

We need to tell R where to store the output. R will make the directory if it does not exist.

```{r}
#### Add the full path to the output directory for storing the analysis.  R will make the directory if it does not exist.
#out_dir <- "/Users/eguswa/Documents/UVM_Projects/long_reads_module/VBRN_Material_For_BI426F24/week_12/lab/analysis/echino_notinf-gonad/"

#If you do not want to assign a full path, a directory will be made in which ever directory you stored this R notebook.
out_dir <- "~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/VIrus_hunting/analysis/"
```

### 

## Step 4. Process the Files and make a Phyloseq object.

An object in R is a way of storing data. A phyloseq object links information across dataframes in a way that is meaningful for analysis for phyloseq.

Make a physeq object

```{r}
# make directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# make otu table by removing taxonomy
otu <- raw %>% select(-otu) 

# make taxonomy table by removing sample counts
tax <- raw %>% select(otu) 

# separate taxonomy
tax_sep <- tax #%>% separate(tax, c("Domain","Kingdom", "Phylum","Class", "Order", "Family", "Genus"), ";")

# read metadata
row.names(samp_data) <- samp_data$sample

# turn otu and tax dataframes to matricies - phyloseq takes matricies 

otum <- as.matrix(otu)
taxm <- as.matrix(tax_sep)


# format the otu table for phyloseq object
OTU = otu_table(otum, taxa_are_rows = TRUE)
# format the tax table for phyloseq object
TAX = tax_table(taxm)
# format the metadata table for phyloseq object
sample_data<- sample_data(samp_data)
# make the phyloseq object and check it out
physeq = phyloseq(OTU, TAX,sample_data)
```

### Check out the new objects.

The OTU variable stores a phyloseq object that gives the abundance of each OTU for each sample.

The TAX variable stores a phyloseq object that gives the full taxonomic path for each OTU. The ranks for each step in the path have their own column.

The sample-data variable stores a phyloseq object that stores the metadata

The phylseq variable stores a phyloseq object that combines the OTU, TAX, sample_data objects.

```{r}
# check out the new ojects.
OTU

TAX

sample_data

physeq
```

## Step 5. Subset the phyloseq object

We will subset by samples by gonad samples for echinostome and not infected.

```{r}
# we will subset the samples to only include echino tissue that is gonad and notif tissue that is gonad.
sub <- 'Molecule'

subset <- subset_samples(physeq, Pair == "Y" )

# get info on new object
subset

#save copy unfiltered for later
subset_div <- subset

#get new counts of reads per sample
sample_sums(subset)

# get a big summary of each sample
summary(subset@otu_table@.Data)

```

## Step 6. Differential Abundance

Let's look at the differential abundance (a statistical methods used to identifies differences in the abundance of microbial taxa between groups of samples) of microorganisms in gonad tissue by infection status.

Before we perform the analysis, using DESeq2, we will filter out OTUs with low counts. We will remove organisms with less than 10 counts across the data set, and that are present in fewer than a quarter of the samples. This will remove false positives from the data set, but could also remove true positives.

Differential abundance here is measured in log fold change (logFC), which is a way to look at the ratio of two values on a logarithmic scale, typically base 2. Positive fold change means something is more abundant that what it is being compared with, and a negative fold change means something is less abundant that what it is compared with.

We consider an OTU as differentially abundant if it as an alpha less than or equal to 0.05. Alpha (aka the significance level) is the probability that we are willing to accept that their results are incorrect. We are willing to accept a 1 in 50 chance that our observations were the result of a chance even in the experiment.

Because we ran a lot of tests in this differential abundance experiment, the value of alpha we are considering is the adjusted P-value. An adjusted P value is a measure of significance that is used to correct for multiple comparison testing.

In this step we will produce and save a table that includes the OTUs with significant differential abundance, the log fold change, and the p value and adjusted P value.

```{r}


sample_data(subset)$Tumor <- as.factor(sample_data(subset)$Tumor)

#sample_data(subset)$sample_type <- as.factor(sample_data(subset)$sample_type)

ds = phyloseq_to_deseq2(subset, ~ Tumor )
nrow(ds)
keep <- rowSums(counts(ds) >= 50) >= as.numeric((ncol(subset@otu_table)-1)/4)

ds <- ds[keep,]

nrow(ds)

ds = DESeq(ds)

alpha = 0.05

resultsNames(ds)

res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig


res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(subset)[rownames(res_sig), ], "matrix"))

#write_csv(as.data.frame(res_sig),paste0(out_dir,"Tumor_mol_sig_differential_abundance_50.csv"))


ggplot(res_sig, aes(x=otu, y=log2FoldChange, color=otu)) +
    geom_jitter(size=3, width = 0.2) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#ggsave(paste0(out_dir,"Tumor_mol_Differentially_Abundant_otu_50.pdf"))


```

In the chunk below we will save a figure of the Family level rank of diferentially abundant OTUs in gonad tissue between echino and no infection individuals.

```{r}

dat <- subset %>% tax_glom(taxrank = "otu") %>% psmelt()

plot1 <- dat %>% group_by(sample, Tumor) %>% 
        reframe( Homo = Abundance[otu == "Homo"]
                  ) %>%
        pivot_longer(cols = -c(sample, Tumor), names_to = "otu", values_to = "Abundance") %>% 

ggplot(aes(Tumor, Abundance, fill = Tumor)) + geom_violin() + geom_jitter(height = 0, width = 0.1) + facet_wrap(~ otu, nrow = 1) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


plot1

ggsave(paste0(out_dir,"DNA_tumor_sig_differential_abundance_plot_Pseudomonas.pdf"),plot=plot1)

Cupriavidus = Abundance[otu == "Cupriavidus"]
Rhizobium = Abundance[otu == "Rhizobium"]
Serratia = Abundance[otu == "Serratia"]
Streptococcus = Abundance[otu == "Streptococcus"]
Streptomyces = Abundance[otu == "Streptomyces"]
Pseudoalteromonas = Abundance[otu == "Pseudoalteromonas"]
Klebsiella = Abundance[otu == "Klebsiella"]
Fischerella = Abundance[otu=="Fischerella"]
Escherichia = Abundance[otu=="Escherichia"]
Akkermansia = Abundance[otu == "Akkermansia"]
Acinetobacter= Abundance[otu=="Acinetobacter"]
Proteiniphilum = Abundance[otu=="Proteiniphilum"]
Aeromonas = Abundance[otu == "Aeromonas"]
Cutibacterium = Abundance[otu == "Cutibacterium"] 
Pseudomonas = Abundance[otu == "Pseudomonas"]
Salmonella = Abundance[otu == "Salmonella"]
Homo = Abundance[otu == "Homo"]
Stigmatella = Abundance[otu=="Stigmatella"]
Citromicrobium = Abundance[otu=="Citromicrobium"]
Blastococcus = Abundance[otu=="Blastococcus"]
Proteiniphilum= Abundance[otu=="Proteiniphilum"]

```
