######################################
# Generate VCF phylogenetic trees in R 
######################################

rm(list = ls())
options(java.parameters="-Xmx50G")

### load libraries

library(fastreeR)
library(ape)
library(rJava)

### increase java resources
.jinit()
num_gigs_ram_available = .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory") / 1e9 
paste("You have ", round(num_gigs_ram_available, 2), "GB memory available.", sep = "")

### indicate vcf file
vcfFile <- "path_to_clean_vcf_dataset/*_0.9_10x.vcf.gz"

### make tree from VCF
vcftree<- vcf2tree(
    vcfFile,
    threads = 10,
    ignoreMissing = FALSE,
    onlyHets = FALSE,
    ignoreHets = FALSE 
  )
  
vcftree

### format output as tree file
tree <- ape::read.tree(text = vcftree)

### import meta data
tree_labels <- read.csv("path_to_sample/meta.csv",  header = TRUE)

## format of the metadatafile:
#              name sample          Location            Tumor Pair          site MultiMito tumor_type percent_mito  libtype TumorP oldcall   call 
# 1  BS1_Hap1_2N_vcf    BS1     Lake_Bomoseen Normal_Haplotype <NA> Lake_Bomoseen         N       <NA>          100       SL   <NA>     BS1 BS1-SL
# 2  BS4_Hap1_2N_vcf    BS4     Lake_Bomoseen Normal_Haplotype <NA> Lake_Bomoseen         N       <NA>          100       SL   <NA>     BS4 BS4_SL
# 3  FB1_Hap1_2N_vcf  FB1.2 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         Y       <NA>           75 ShortDNA   <NA>   FB1.2  FB1.2
# 4  FB1_Hap2_2N_vcf  FB1.1 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         Y       <NA>           25 ShortDNA   <NA>   FB1.1  FB1.1
# 5 FB11_Hap1_2N_vcf   FB11 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         N       <NA>          100 ShortDNA   <NA>    FB11   FB11
# 6 FB12_Hap1_2N_vcf   FB12 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         N       <NA>          100 ShortDNA   <NA>    FB12   FB12
# ....

## join the tree data with metadata
x <- left_join(as_tibble(tree), tree_labels, by = c("label" = "name"))  %>% mutate(name=label)


## make tree including clade lables.  Will have to change clade labels for data set.  Will need to find the nodes for each set of samples by uncommenting #geom_label(aes(label=node))


p9 <- ggtree(tree) %<+% x  + 
  geom_tippoint(aes(color = TumorP))  + coord_cartesian(clip = 'off') + 
  geom_tiplab(aes(label=sample, color=TumorP, fontface="bold"),size=3, align=TRUE, linesize=.1) +
  geom_nodelab( aes(label=label, subset=label > 79), hjust = 1.44, nudge_y = 0.65, size = 3) +
  geom_nodelab( aes(label=label, subset=label ==100 ), hjust = 1.44, nudge_y = 0.65, size = 3) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(10, 200, 10, 10))  + 
  scale_color_manual(values=c("#3A9AB2","#F11B00","grey40")) + theme(legend.position="none") +
  geom_cladelab(node=100, label="Lake Memphramegog HF Tumor Samples", align=TRUE, 
                textcolor='black', barsize=2,barcolor="#3A9AB2", fontsize =2.5,offset=.03) + 
  geom_cladelab(node=106, label="Reference Fish", align=TRUE, 
                textcolor='black', barcolor="#000000", barsize=2,fontsize =2.5,offset=.03) + 
  geom_cladelab(node=113, label="Lake Memphramegog Normal Samples", align=TRUE, barsize=2,
                textcolor='black', barcolor="#F11B00",  fontsize =2.5,offset=.03) +
  theme(axis.text.x = element_text(face="bold", size = 10)) 

p9
