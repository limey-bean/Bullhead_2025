---
title: "Phase 2N Mito vcf by Allele Frequencies"
author: "Emily Curd"
output: pdf_document
---

```{r setup}

library(vcfR)
library(stringr)
library(dplyr)
library(tidyverse)
library(hoardeR)
library(bedr)
library("Biostrings")
library(vtable)

```

### Set paths / variables
```{r}
outpath <- "path to output dir"

dir.create(outpath, showWarnings = FALSE)

vcf_file <- "path_to_input/*_TN_mito_snps-and-indels_clean_0.9_10x.vcf.gz"

mito_fasta <- "path_to_input/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa"

```

### Parse vcf
```{r}
vcf <- read.vcfR(vcf_file, verbose = FALSE )

G <- vcfR2tidy(vcf, format_fields = c("GT", "AD", "DP"))
base <- G$fix %>% select(POS, REF)

m <- G$gt %>% mutate(across('gt_GT', str_replace, '\\|', '/')) %>% mutate(across('gt_GT_alleles', str_replace, '\\|', '/'))

freq_allele <- base %>% select(POS)

f <- G$gt$gt_AD

```

### Make a dataframe from the mito genome fasta file.    
```{r}

fastaFile <- readDNAStringSet(mito_fasta)
sequence = paste(fastaFile)
REFa <- strsplit(sequence, "")[[1]] 
# Create a data frame
REF_seq_df <- as.data.frame(REFa, stringsAsFactors = FALSE)
# View the result
REF_seq_df <- REF_seq_df %>% mutate(POS = 1:n()) %>% select(POS,REFa) %>% filter(!POS %in% c(15815,15674,11734,11733,11732,11731,11730,11729,11728,11727,11726,11725,11724,11723,11722,2767,2766,2765,2764,2763,945,944)) # filter ref specific indels

```

### Uncomment the individuals list you are processing
```{r}
## short read DNA samples
individuals <- c("FB1","FB15","FB5","HC1","JR13","JR14","LB1","MA15","MA18","SC11","SML11","SML15","SML3","BS1","BS4","FB11","FB12","FB13","FB16","FB18","FB2","FB21","FB3","FB8N","HC10","HC11","HC12","HC14","HC17","HC18","HC19","HC27","HC2N","HC3","HC4N","HC4T","HC7","HL3","HL4","JR1","JR10","JR11","JR12","JR18","JR19","JR2","JR3","JR4","JR9N","LB2","MA10","MA11","MA12","MA13","MA14","MA1N","MA5","MA7","MA8","MA9","SB1","SB10","SB12","SB13","SB14N","SB15","SB16","SB2","SB3","SB4","SB5","SB9","SC1","SC12","SC13","SC16","SC17","SC3","SC4","SC5","SC6N","SC7N","SML12","SML13","SML14","SML16","SML1N","SML2","SML4N","SML5","SML7","SP2","SP3")

## short read RNA samples
#individuals <- c("vt1_02","vt1_02","vt1_03","vt1_05","vt1_07","vt2_02","vt2_05","vt2_05","vt2_06","vt2_10","vt2_14","vt2_16","vt2_18","vt2_23","vt2_23","vt2_28","vt2_29R","vt2_30","vt2_30","vt2_32","vt2_37","vt2_40")

## Long read DNA samples
#individuals <- c("BS1")#,"BS4","FB8N","FB8T","HC14","HC2N","HC2T","HC4N","HC4T","HL3","HL4","JR11","JR19","JR9N","JR9T","LB1","LB2","MA18","MA1N","MA1T","SB14N","SB14T","SC1","SC6N","SC6T","SC7N","SC7T","SML13","SML1N","SML1T","SML4N","SML4T","SP2","SP3")
```

### Determine if the sample is Homozygous, Heterozygous and clean (e.g. Variant Allele Freqeuncy is <40 or >60), or Heterozygous and messy (e.g. Variant Allele Freqeuncy is 40 > or < 60).  If it is messy and a tumor sample you can use the normal to phase mitochondria.  If you do not have a paired sample (e.g. RNAseq) you have to make some judgement calls.  It it is homozygous or clean heterozygous, it will pop out a fasta or two fasta files for each sample.

```{r}
hom_list <- "is_hom"
het_list <- "is_het"
messy <- "is_messy"

for (ind in individuals){  
  print(ind)

  major_alleleN <- paste0(ind,"_Hap1")
  minor_alleleN <- paste0(ind,"_Hap2")
  major_allele_freq <- paste0(ind,"_Hap1_freq")
  minor_allele_freq <- paste0(ind,"_Hap2_freq")
  
  n <- m %>% filter(Indiv==ind)
  n2 <- n %>% separate_wider_delim(gt_GT_alleles, delim = "/" , names = c("a0", "a1")) 
  
  n3 <- n2 %>% separate_wider_delim(gt_AD, delim = "," , names = c("a0_AD", "a1_AD", "other_AD"),too_few = "align_start", too_many= "merge")
  
  n4 <- n3 %>% mutate(a0_AF = as.numeric(a0_AD) / gt_DP) %>% mutate(a1_AF = as.numeric(a1_AD) / gt_DP) %>% mutate(a0_AF = replace_na(a0_AF, 1))  %>% mutate(a1_AF = replace_na(a1_AF, 0))
  
  n5 <- n4 %>% mutate(!!major_alleleN := case_when(a0_AF == 1 ~ a0,
                                                a1_AF == 1 ~ a1,
                                                a0_AF < 1 & a0_AF >= .5  ~ a0,
                                                a1_AF < 1 & a1_AF >= .5 ~ a1,
                                                a0_AF <= 0.5 & a1_AF <= 0.5 & a0_AF < a1_AF ~ a1,
                                                a0_AF <= 0.5 & a1_AF <= 0.5 & a1_AF < a0_AF ~ a0)) %>% 
    mutate(!!minor_alleleN := case_when(a0_AF == 1 ~ a1,
                                     a1_AF == 1 ~ a0,
                                     a0_AF < 1 & a0_AF >= .5  ~ a1,
                                     a1_AF < 1 & a1_AF >= .5 ~ a0,
                                     a0_AF <= 0.5 & a1_AF <= 0.5 & a0_AF < a1_AF ~ a0, 
                                     a0_AF <= 0.5 & a1_AF <= 0.5 & a1_AF < a0_AF ~ a1)) %>%
    mutate(agree = case_when( !!as.name(major_alleleN) == !!as.name(minor_alleleN) ~ "Y",
                              !!as.name(major_alleleN) != !!as.name(minor_alleleN) & abs(a1_AF - a0_AF) > 0.10  ~ "N",
                              !!as.name(major_alleleN) != !!as.name(minor_alleleN) & abs(a1_AF - a0_AF) < 0.10  ~ "U")) 
  
  
   if ((sum(str_detect(n5$agree, 'N')) == 0 & (sum(str_detect(n5$agree, 'U')) == 0) == TRUE)) {
    print(paste0(ind," is homozygous"))
     hom_list <- c(hom_list,ind)
     n_allele <- n5 %>% select(POS, !!major_alleleN)
     hom_seq <- left_join(REF_seq_df,n_allele)
     hom_seq <- hom_seq %>% mutate(new_fasta = case_when(!is.na(!!as.name(major_alleleN)) 
                                                         ~ !!as.name(major_alleleN),TRUE ~ REFa))
     
     seq_name <- paste0(">L",ind, "_Hap1_2N_vcf")
     new_seq <- as.character(hom_seq$new_fasta)
     new_seq <- paste(new_seq,collapse="")
     fasta <- character(2)
     fasta[c(TRUE, FALSE)] <- seq_name
     fasta[c(FALSE, TRUE)] <- new_seq
  
     writeLines(fasta, file.path(paste0(outpath, ind, "_Hap1_2N_vcf.fasta")))
                
     allele_freq <- n5 %>% mutate( !!major_allele_freq := case_when(a0_AF > .93 ~ a0_AF + a1_AF,
                                                        a1_AF > .93 ~ a0_AF + a1_AF
                                                        )) %>% 
                           select(POS, !!as.name(major_allele_freq))  
     
     freq_allele <- left_join(freq_allele, allele_freq)
     
     
   } else if ((sum(str_detect(n5$agree, 'N')) > 0 & (sum(str_detect(n5$agree, 'U')) == 0) == TRUE)) {
     
     het_list <- c(het_list,ind)
     print(paste0(ind," has hets"))
     n_alleles <- n5 %>% select(POS, !!major_alleleN, !!minor_alleleN) 
     het_seq <- left_join(REF_seq_df,n_alleles)
     het_seq <- het_seq %>% mutate(hap1_fasta = case_when(!is.na(!!as.name(major_alleleN)) 
                                                         ~ !!as.name(major_alleleN),TRUE ~ REFa)) %>% mutate(hap2_fasta = case_when(!is.na(!!as.name(minor_alleleN)) 
                                                         ~ !!as.name(minor_alleleN),TRUE ~ REFa))
     
     seq_name1 <- paste0(">L",ind, "_Hap1_2N_vcf")
     new_seq1 <- as.character(het_seq$hap1_fasta)
     new_seq1 <- paste(new_seq1,collapse="")
     fasta1 <- character(2)
     fasta1[c(TRUE, FALSE)] <- seq_name1
     fasta1[c(FALSE, TRUE)] <- new_seq1
     writeLines(fasta1, file.path(paste0(outpath, ind, "_Hap1_2N_vcf.fasta")))
  
     
     seq_name2 <- paste0(">L",ind, "_Hap2_2N_vcf")
     new_seq2 <- as.character(het_seq$hap2_fasta)
     new_seq2 <- paste(new_seq2,collapse="")
     fasta2 <- character(2)
     fasta2[c(TRUE, FALSE)] <- seq_name2
     fasta2[c(FALSE, TRUE)] <- new_seq2
     writeLines(fasta2, file.path(paste0(outpath, ind, "_Hap2_2N_vcf.fasta")))
     
     
     allele_freq <- n5 %>%  mutate(!!major_allele_freq := case_when(
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) == a1 & a0_AF > a1_AF ~ a0_AF,
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) == a1 & a0_AF < a1_AF ~ a1_AF,
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) != a1 ~ a0_AF,
      !!as.name(major_alleleN) != a0 & !!as.name(major_alleleN) == a1 ~ a1_AF,
    )) %>%
    mutate(!!minor_allele_freq := case_when(
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) == a1 & a0_AF > a1_AF ~ a0_AF,
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) == a1 & a0_AF < a1_AF ~ a1_AF,
      !!as.name(major_alleleN) == a0 & !!as.name(major_alleleN) != a1 ~ a1_AF,
      !!as.name(major_alleleN) != a0 & !!as.name(major_alleleN) == a1 ~ a0_AF)) %>% 
                           select(POS, !!as.name(major_allele_freq), !!as.name(minor_allele_freq))
     
     freq_allele <- left_join(freq_allele, allele_freq)  
     
  
   } else {
     
     messy <- c(messy,ind)
    
     print(paste0("Run subtract-normal code on sample: ", ind ))
     
   }
  
}



```



### Subtract host reads from tumor tissue using normal brain tissue - use for messy paired tumors or any paired tumors. It will pop out a fasta or two fasta files for each sample. Note - none of the brain samples in the BBH study had SNVs and it might behave strangley if your host sample has an SNV.

```{r}
individualsT <- c("FB8","HC2","JR9","MA1","SB14","SC6","SC7","SML1","SML4","HC4")

for (ind in individualsT){  
  print(ind)

  indN <- paste0(ind,"N")
  indT <- paste0(ind,"T")
  
  
  freqN <- paste0(indT, "_normal_haplotype") 
  freqTN <- paste0(indT, "_TN_haplotype") 
  
  n <- m %>% filter(Indiv==indN)
  n2 <- n %>% separate_wider_delim(gt_GT_alleles, delim = "/" , names = c("a0", "a1")) 
  
  n3 <- n2 %>% separate_wider_delim(gt_AD, delim = "," , names = c("a0_AD", "a1_AD", "other_AD"),too_few = "align_start", too_many= "merge")
  
  n4 <- n3 %>% mutate(a0_AF = as.numeric(a0_AD) / gt_DP) %>% mutate(a1_AF = as.numeric(a1_AD) / gt_DP)
  
  n5 <- n4 %>% mutate(brain_major_allele = case_when(a0_AF == 1 ~ a0,
                                                a1_AF == 1 ~ a1,
                                                a0_AF < 1 & a0_AF > .5  ~ a0,
                                                a1_AF < 1 & a1_AF > .5 ~ a1)) %>% 
    mutate(brain_minor_allele = case_when(a0_AF == 1 ~ a1,
                                     a1_AF == 1 ~ a0,
                                     a0_AF < 1 & a0_AF > .5  ~ a1,
                                     a1_AF < 1 & a1_AF > .5 ~ a0)) %>%
    mutate(agree = case_when(brain_major_allele == brain_minor_allele ~ "Y",
                             brain_major_allele != brain_minor_allele ~ "N")) 
  
  
  n6 <- n5 %>% select(c(POS,brain_major_allele, brain_minor_allele))
  
  
  t <- m %>% filter(Indiv==indT)
  
  t2 <- t %>% separate_wider_delim(gt_GT_alleles, delim = "/" , names = c("a0", "a1","a_other"),too_few = "align_start", too_many= "merge")
  
  t3 <- t2 %>% separate_wider_delim(gt_AD, delim = "," , names = c("a0_AD", "a1_AD", "other_AD"),too_few = "align_start", too_many= "merge")
  
  
  t4 <- t3 %>% mutate(a0_AF = as.numeric(a0_AD) / gt_DP) %>% mutate(a1_AF = as.numeric(a1_AD) / gt_DP)
  
  t5 <- t4 %>% mutate(major_alleleT = case_when(a0_AF == 1 ~ a0,
                                                a1_AF == 1 ~ a1,
                                                a0_AF < 1 & a0_AF >= .5  ~ a0,
                                                a1_AF < 1 & a1_AF > .5 ~ a1)) %>% 
    mutate(minor_alleleT = case_when(a0_AF == 1 ~ a1,
                                     a1_AF == 1 ~ a0,
                                     a0_AF < 1 & a0_AF >= .5  ~ a1,
                                     a1_AF < 1 & a1_AF > .5 ~ a0)) %>%
    mutate(agree = case_when(major_alleleT == minor_alleleT ~ "Y",
                             major_alleleT != minor_alleleT ~ "N"))
  
  t6 <- t5 %>% select(c(POS, major_alleleT, minor_alleleT))
  
  
  tn_alleles <- full_join(t6,n6) %>%
  mutate(agreeN = case_when(brain_major_allele == brain_minor_allele ~ "Y",
                            brain_major_allele != brain_minor_allele ~ "N")) %>%
    mutate(agreeT = case_when(major_alleleT == minor_alleleT ~ "Y",
                              major_alleleT != minor_alleleT ~ "N")) %>%
    mutate(agreeNt = case_when(brain_major_allele == minor_alleleT ~ "Y",
                               brain_major_allele != minor_alleleT ~ "N")) %>% 
    mutate(agreeTn = case_when(major_alleleT == brain_minor_allele ~ "Y",
                               major_alleleT != brain_minor_allele ~ "N"))
  
  
  subtract_alleles <- full_join(t5,n6, by="POS")
  
  
  tn_bean_counting <- subtract_alleles %>%
    mutate(N = case_when(
      brain_major_allele == a0 & brain_major_allele == a1 ~ a0,
      brain_major_allele == a0 & brain_major_allele != a1 ~ a0,
      brain_major_allele != a0 & brain_major_allele == a1 ~ a1,
    )) %>%
    mutate(TN = case_when(
      brain_major_allele == a0 & brain_major_allele == a1 ~ a0,
      brain_major_allele == a0 & brain_major_allele != a1 ~ a1,
      brain_major_allele != a0 & brain_major_allele == a1 ~ a0,
    )) %>%
    mutate(!!as.name(freqN) := case_when(
      brain_major_allele == a0 & brain_major_allele == a1 & a0_AF > a1_AF ~ a0_AF,
      brain_major_allele == a0 & brain_major_allele == a1 & a0_AF < a1_AF ~ a1_AF,
      brain_major_allele == a0 & brain_major_allele != a1 ~ a0_AF,
      brain_major_allele != a0 & brain_major_allele == a1 ~ a1_AF,
    )) %>%
    mutate(!!as.name(freqTN) := case_when(
      brain_major_allele == a0 & brain_major_allele == a1 & a0_AF > a1_AF ~ a0_AF,
      brain_major_allele == a0 & brain_major_allele == a1 & a0_AF < a1_AF ~ a1_AF,
      brain_major_allele == a0 & brain_major_allele != a1 ~ a1_AF,
      brain_major_allele != a0 & brain_major_allele == a1 ~ a0_AF,
    )) 
  
  tn_m_alleles <- tn_bean_counting %>%  mutate(major_allele_is_N = case_when(major_alleleT == N ~ "Y",
                              major_alleleT != N ~ "N")) 
  tn_m_allele_freq <- tn_bean_counting %>% select(POS, !!as.name(freqN), !!as.name(freqTN))
  
  freq_allele <- left_join(freq_allele, tn_m_allele_freq)
  
  TN <- tn_m_alleles %>% select(POS,TN,N)
  
  seq_df_N <- left_join(REF_seq_df,n6)
  seq_df_N <- seq_df_N %>% mutate(new_fasta = case_when(!is.na(brain_major_allele) ~ brain_major_allele, TRUE ~ REFa))
  
  seq_df_N_TNN <-left_join(seq_df_N,TN)
  
  seq_df_N_TNN <- seq_df_N_TNN %>% mutate(TNfasta = case_when(!is.na(TN) ~ TN,
                                                        TRUE ~ REFa)) %>%
                          mutate(hostfasta = case_when(!is.na(N) ~ N,
                               TRUE ~ REFa))
       
       seq_name <- paste0(">",indT, "_normal_Hap1_2N_vcf")
       new_seq <- as.character(seq_df_N_TNN$hostfasta)
       new_seq <- paste(new_seq,collapse="")
       fasta1 <- character(2)
       fasta1[c(TRUE, FALSE)] <- seq_name
       fasta1[c(FALSE, TRUE)] <- new_seq
       writeLines(fasta1, file.path(paste0(outpath, indT, "_normal_Hap1_2N_vcf.fasta")))
  
  
       seq_name <- paste0(">",indT, "_TN_Hap2_2N_vcf")
       new_seq <- as.character(seq_df_N_TNN$TNfasta)
       new_seq <- paste(new_seq,collapse="")
       fasta2 <- character(2)
       fasta2[c(TRUE, FALSE)] <- seq_name
       fasta2[c(FALSE, TRUE)] <- new_seq
  
       writeLines(fasta2, file.path(paste0(outpath, indT, "_TN_Hap2_2N_vcf.fasta")))
  
}  

```
### save output

```{r}
write_csv(freq_allele, paste0("path to output directory/ALL_2N_vcf_ref_allele_frequencies.csv"))
```

### review the output and detemine if there are methodological artifacts. If there are not paired samples for messy heterozygous loci make decisions about calling SNVs.  Modify the fastas if there is a valid reason to do so.
