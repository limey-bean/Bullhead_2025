---
title: "R Notebook"
output: html_notebook
---

## Load R packages
```{r}
library(stringr)
library(dplyr)
library(reshape)
library(wesanderson)
library(tidyverse)
```


## Load data and modify header depeding of how the data were generated

```{r}


sv_sites <- read.table("path_to_file/*vcf.AD.txt", sep = "\t", header = TRUE, comment.char = "") 

# Uncomment below if using Sniffles
#colnames(sv_sites) <- colnames(sv_sites) %>%
#  str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
#  str_remove("\\.DR.+DV$")                # remove trailing ".AD"

# Uncomment below if using Delly 
#colnames(sv_sites) <- colnames(sv_sites) %>%
#  str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
#  str_remove("\\.RR.+RV$")                # remove trailing ".AD"

# if using Haplotype caller
colnames(sv_sites) <- colnames(sv_sites) %>%
  str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
  str_remove("\\.AD$")                # remove trailing ".AD"

sv_sites <- sv_sites 

```

## Make two color palettes

```{r}
pal4 <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 2, type = "continuous")
```


## Split allele depth by Reference sv_sites3 and Variant alleles sv_sites2
```{r}

sv_sites2 <- sv_sites  %>%
  mutate(across(
    -c(CHROM, POS, REF, ALT),  # adjust to actual column names to exclude
    ~ .x %>%
      str_replace_all("^\\.$", "0,0") %>%                     # replace "." with "0,0"
      str_split(",", simplify = TRUE) %>%                     # split on comma
      magrittr::extract(, 2) %>%                              # get second column
      as.numeric()
  )) 

sv_sites3 <- sv_sites  %>%
  mutate(across(
    -c(CHROM, POS, REF, ALT),  # adjust to actual column names to exclude
    ~ .x %>%
      str_replace_all("^\\.$", "0,0") %>%                     # replace "." with "0,0"
      str_split(",", simplify = TRUE) %>%                     # split on comma
      magrittr::extract(, 1) %>%                              # get second column
      as.numeric()
  )) 

```

## Find all sites were there are Tumor samples's have a call but Ref is 0 or NA

```{r}
# filter for all Ns to be 0, all Ts to be nonzero
allC_noHt <- sv_sites2 %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  )) %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  ))

allC_noH <- sv_sites2 %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  )) %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
         ) %>% 
  select(CHROM,POS,REF,ALT, COUNT) %>%
  mutate(ID = "allC")

allC_noH_counts <- allC_noH %>% count(COUNT)

allC_noH_reft <- sv_sites3 %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  ))   %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
  )

allC_noH_ref <- sv_sites3 %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  ))   %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
  ) %>% 
  select(CHROM,POS,REF,ALT, COUNT) %>%
  mutate(ID = "allC_ref")

allC_noH_ref_counts <- allC_noH_ref %>% count(COUNT)

```


## Find all sites were there are Ref samples's have a call but Tumor is 0 or NA

```{r}


#### normal

allN_noH <- sv_sites2 %>%
  filter(if_all(
    .cols =!ends_with("N") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  )) %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
  ) %>% 
  select(CHROM,POS,REF,ALT, COUNT) %>%
  mutate(ID = "allN")

allN_noH_counts <- allN_noH %>% count(COUNT)

allN_noH_ref <- sv_sites3 %>%
  filter(if_all(
    .cols = !ends_with("N") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  ))   %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
  ) %>% 
  select(CHROM,POS,REF,ALT, COUNT) %>%
  mutate(ID = "allN_ref")

allN_noH_ref_counts <- allN_noH_ref %>% count(COUNT)

```


## Count site with normal only or tumor only variants

```{r}


normalsites <- full_join(allN_noH_counts, allN_noH_ref_counts, by="COUNT") %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>% mutate(Normal_Samples_only = n.y + n.x) %>% select(COUNT, Normal_Samples_only) 
  
tumorsites <- full_join( allC_noH_ref_counts,allC_noH_counts, by="COUNT") %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>% mutate(Tumor_Samples_only = n.y + n.x) %>% select(COUNT, Tumor_Samples_only)

#tumorsites <- allC_noH_counts %>% mutate(Tumor_Samples_only = n) %>% select(COUNT, Tumor_Samples_only)
  
onlysite <- full_join(normalsites,tumorsites, by="COUNT") %>% select(COUNT,Tumor_Samples_only, Normal_Samples_only) #%>% add_row(COUNT=3,Tumor_Samples_only = 0, Normal_Samples_only = 0) %>% add_row(COUNT=10,Tumor_Samples_only = 0, Normal_Samples_only = 0)
onlysite[is.na(onlysite)] <- 0

only <- onlysite[order(as.numeric(as.character(onlysite$COUNT))), ]

site <- only %>% select(-COUNT)

```

## Make figures, and check that sites have counts for all numbers of samples.  If there is a zero count you may have to add rows to only site (e.g. commented out text above)

```{r}
meltData <- melt(site) 


meltData$shared_by <- c("01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples","09_samples","10_samples","01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples","09_samples","10_samples")




ggplot(meltData, aes(fill=shared_by, y=value, x=variable)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 



ggplot(meltData, aes(fill=shared_by, y=value, x=variable)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal2)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal2)  + geom_bar(position="dodge", stat="identity") 



```

