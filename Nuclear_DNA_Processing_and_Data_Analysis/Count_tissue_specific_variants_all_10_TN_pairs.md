---
title: "R Notebook"
output: html_notebook
---

```{r}
library(stringr)
library(dplyr)
library(reshape)
library(wesanderson)
library(tidyverse)
```

```{r}


#mutectT <- read.table("~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/mito_N1_N2/vcf_parse_2N/mito_2n_TN.vcf.AD.txt", sep = "\t", header = TRUE, comment.char = "") 

#mutectT <- read.table("longSV/sniffles_al_INV_AD.txt", sep = "\t", header = TRUE, comment.char = "") 

#mutectT <- read.table("latest_delly/2025_latest/merged_130M.PASS.inv_.9_10x.TN.ad.txt", sep = "\t", header = TRUE, comment.char = "") 

#mutectT <- read.table("mutect2-nuc/all_hap_mut.ad.txt", sep = "\t", header = TRUE, comment.char = "") 

#mutectT <- read.table("mutect2-nuc/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.ad.txt", sep = "\t", header = TRUE, comment.char = "") 


colnames(mutectT) <- colnames(mutectT) %>%
  str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
  str_remove("\\.DR.+DV$")                # remove trailing ".AD"

#colnames(mutectT) <- colnames(mutectT) %>%
#  str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
#  str_remove("\\.RR.+RV$")                # remove trailing ".AD"

#colnames(mutectT) <- colnames(mutectT) %>%
 # str_remove("^X\\.\\.?\\d+\\.") %>%  # remove leading "X.2." or "X..1."
  #str_remove("\\.AD$")                # remove trailing ".AD"

mutectT <- mutectT #%>% select(-FB8N,-FB8T,-JR9N, -JR9T,-SB14T,-SB14N,-HC4N,-HC4T)

```

```{r}
pal4 <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 2, type = "continuous")
```

```{r}

mutect2T <- mutectT  %>%
  mutate(across(
    -c(CHROM, POS, REF, ALT),  # adjust to actual column names to exclude
    ~ .x %>%
      str_replace_all("^\\.$", "0,0") %>%                     # replace "." with "0,0"
      str_split(",", simplify = TRUE) %>%                     # split on comma
      magrittr::extract(, 2) %>%                              # get second column
      as.numeric()
  )) 

mutect3T <- mutectT  %>%
  mutate(across(
    -c(CHROM, POS, REF, ALT),  # adjust to actual column names to exclude
    ~ .x %>%
      str_replace_all("^\\.$", "0,0") %>%                     # replace "." with "0,0"
      str_split(",", simplify = TRUE) %>%                     # split on comma
      magrittr::extract(, 1) %>%                              # get second column
      as.numeric()
  )) 

#mito
#mutect2T <- mutect2T%>%
  #mutate(across(.cols= c(5:24),.fns = function(x) ifelse(x <= 6,0,x)))

#mutect3T <- mutect3T%>%
  #mutate(across(.cols= c(5:24),.fns = function(x) ifelse(x <= 6,0,x)))


```

```{r}
# filter for all Ns to be 0, all Ts to be nonzero
allC_noHt <- mutect2T %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  )) %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  ))

allC_noH <- mutect2T %>%
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

allC_noH_reft <- mutect3T %>%
  filter(if_all(
    .cols = !ends_with("T") & !c("CHROM", "POS", "REF", "ALT"),  # Exclude CHROM, POS, REF, ALT manually
    .fns = ~ .x == 0 | is.na(.x)                                 # Keep only if zero or NA
  ))   %>%
  mutate(COUNT=rowSums(.[5:24]!=0
  #mutate(COUNT=rowSums(.[5:16]!=0
  )) %>%
  filter(COUNT != 0
  )

allC_noH_ref <- mutect3T %>%
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

```{r}


#### normal

allN_noH <- mutect2T %>%
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

allN_noH_ref <- mutect3T %>%
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

```{r, eval=FALSE}
pos1 <- allC_noHt %>% filter(COUNT > 0) %>% select("POS","COUNT")
pos2 <- allC_noH_reft %>% filter(COUNT > 0) %>% select("POS","COUNT")

pos <- rbind(pos1,pos2) %>% arrange

mutectTtumorsitesRef <- left_join(pos,mutect2T)
mutectTtumorsitesALT <- left_join(pos,mutect3T) 

Tref <- mutectTtumorsitesRef %>% select(c(6:25))
Tref <- as.matrix(Tref)
Talt <- mutectTtumorsitesALT %>% select(c(6:25))
Talt <- as.matrix(Talt)



Ttot <- Talt + Tref
TAltfreq <- Talt / Ttot
TReffreq <- Tref / Ttot

TAltfreq <- as.data.frame(TAltfreq)
TReffreq <- as.data.frame(TReffreq)


TAltfreq <- cbind(pos,as.data.frame(TAltfreq))
TReffreq <- cbind(pos,as.data.frame(TReffreq))





ref_int <- TReffreq %>% filter(POS == 1718 | POS == 4626 | POS == 12895)

tumorsitesfreq <- rbind(TAltfreq, ref_int)

write.csv(tumorsitesfreq, file="~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/mito_N1_N2/vcf_parse_2N/may2025-92/tumorsites_freq.csv")

tumAD <- left_join(pos,mutectT)

write_delim(tumAD, file="~/Documents/UVM_Projects/Bullhead_transmissible_cancer_adventure/mito_N1_N2/vcf_parse_2N/may2025-92/tumorsites_AD.tsv", col_names =TRUE, delim="\t")

```

```{r, eval=FALSE}
all_pos <- mutect2T %>% select("CHROM","POS")

ref_all <- mutect2T %>% select(c(5:24)) 
alt_all <- mutect3T %>% select(c(5:24)) 
ref_all <- as.matrix(ref_all) 
alt_all <- as.matrix(alt_all)

tot <- alt_all + ref_all 
Altfreq <- alt_all / tot 
Reffreq <- ref_all / tot

Altfreq <- cbind(all_pos,as.data.frame(Altfreq)) 
Reffreq <- cbind(all_pos,as.data.frame(Reffreq))

plot(Altfreq$SC7N ~ Altfreq$SC7T)
```

```{r}


normalsites <- full_join(allN_noH_counts, allN_noH_ref_counts, by="COUNT") %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>% mutate(Normal_Samples_only = n.y + n.x) %>% select(COUNT, Normal_Samples_only) 
  
tumorsites <- full_join( allC_noH_ref_counts,allC_noH_counts, by="COUNT") %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>% mutate(Tumor_Samples_only = n.y + n.x) %>% select(COUNT, Tumor_Samples_only)

#tumorsites <- allC_noH_counts %>% mutate(Tumor_Samples_only = n) %>% select(COUNT, Tumor_Samples_only)
  
onlysite <- full_join(normalsites,tumorsites, by="COUNT") %>% select(COUNT,Tumor_Samples_only, Normal_Samples_only) #%>% add_row(COUNT=3,Tumor_Samples_only = 0, Normal_Samples_only = 0) %>% add_row(COUNT=10,Tumor_Samples_only = 0, Normal_Samples_only = 0)
onlysite[is.na(onlysite)] <- 0

only <- onlysite[order(as.numeric(as.character(onlysite$COUNT))), ]

site <- only %>% select(-COUNT)

```

```{r}
meltData <- melt(site) 


meltData$shared_by <- c("01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples","09_samples","10_samples","01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples","09_samples","10_samples")

#meltData$shared_by <- c("01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples","01_sample","02_samples","03_samples","04_samples","05_samples","06_samples","07_samples","08_samples")


ggplot(meltData, aes(fill=shared_by, y=value, x=variable)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 



ggplot(meltData, aes(fill=shared_by, y=value, x=variable)) + scale_fill_manual(values = pal4)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal2)  + geom_bar(position="stack", stat="identity") 

ggplot(meltData, aes(fill=variable, y=value, x=shared_by)) + scale_fill_manual(values = pal2)  + geom_bar(position="dodge", stat="identity") 




RVNnsnp <- allN_noH %>% select(CHROM,POS)
RRNnsnp <-  allN_noH_ref %>% select(CHROM,POS)
Nsites <- rbind(RVNnsnp,RRNnsnp)


RVTnsnp <- allC_noH %>% select(CHROM,POS)
RRTnsnp <- allC_noH_ref %>% select(CHROM,POS)
Tsites <- rbind(RVTnsnp,RRTnsnp)

filter_these_sites <- rbind(Nsites,Tsites) %>% group_by(CHROM) %>% distinct( POS, .keep_all = TRUE)


allsites_t <- full_join(allC_noH, allC_noH_ref, join_by(CHROM,POS,REF,ALT))

allsites_n <- full_join(allN_noH, allN_noH_ref, join_by(CHROM,POS,REF,ALT))

allsites_t_n <- full_join(allsites_t, allsites_n, join_by(CHROM,POS,REF,ALT))


#write.table(filter_these_sites, "mutect2-nuc/clean_som_PASS-6.txt", quote = FALSE, row.names = FALSE, sep ="\t")

```

```{r}
pos_all <- allsites_t_n %>% select("CHROM","POS")

rsites <- left_join(pos_all,mutect2T) 
asites <- left_join(pos_all,mutect3T)

ref_all <- rsites %>% select(c(5:24)) 
alt_all <- asites %>% select(c(5:24)) 
ref_all <- as.matrix(ref_all) 
alt_all <- as.matrix(alt_all)

tot <- alt_all + ref_all 
Altfreq <- alt_all / tot 
Reffreq <- ref_all / tot

Altfreq <- cbind(pos_all,as.data.frame(Altfreq)) 
Reffreq <- cbind(pos_all,as.data.frame(Reffreq))

plot(Altfreq$SC7N ~ Altfreq$SC7T)
plot(Reffreq$SC7N ~ Reffreq$SC7T)

```

