######################################
# Generate phased mitochondrial phylogenetic trees
######################################

##########
# concatonate the fasta files or interest into one mutlifasta file
##########

cat path_to_fastas/*.fasta > mito_reads.fasta


#########
# Align phased mitochondrial fasta files using mafft
#########

mafft mito_reads.fasta > mito_reads.align.fasta

# view the alignment in an alignment viewer to find method specific SNPs or indels
# make HL4 the first samples in the fasta file.

#########
# Generate phylogenetic tree from an iqtree2 alignment
#########

iqtree -s mito_reads.align.fasta -alrt 1000 -bb 1000  -safe 

# generated the figures using the output mito_reads.align.fasta.contree


########
# get SNP sites from the mutlifasta file
########

snp-sites -o test.snpsites_aln.tsv test.aln.fasta 
snp-sites -v -o test.snpsites_vcf.tsv test.aln.fasta 

# Manually combind the outputs so that the sample names are column names and position in the mito genome are row names. Change "HL_Hap1_2N_vcf" 4 to "ref" so that you have a text file that looks like the following:

#	ref	LBS1_Hap1_2N_vcf	LBS4_Hap1_2N_vcf	LBS4_Hap2_2N_vcf	LFB8N_Hap1_2N_vcf	LFB8T_normal_Hap1_2N_vcf
#29	A	A	A	A	A	A
#364	G	G	G	G	G	G
#584	A	A	A	A	A	A
#627	G	G	G	G	G	G
#663	C	C	C	C	C	C
#848	G	G	A	A	G	G
#938	T	T	T	T	T	T
#1035	C	C	C	C	C	C

########
# Run the following in R
########

library(scico)
library(ggplot2)
library(ggmsa)
library(fastreeR)
library(tidyverse)
library(dplyr)
library(ggnewscale)
library(aplot)
library("wesanderson")
library(scales)

# load datafiles for the tree, meta data, and snps
tree <- read.tree(file="path_to/mito_reads.align.fasta.contree")
tree_labels <- read.csv("path_to/meta.csv",  header = TRUE)
snps <- read.delim("path_to/snps.txt", header = TRUE)

##format of the metadatafile:
#              name sample          Location            Tumor Pair          site MultiMito tumor_type percent_mito  libtype TumorP oldcall   call 
# 1  BS1_Hap1_2N_vcf    BS1     Lake_Bomoseen Normal_Haplotype <NA> Lake_Bomoseen         N       <NA>          100       SL   <NA>     BS1 BS1-SL
# 2  BS4_Hap1_2N_vcf    BS4     Lake_Bomoseen Normal_Haplotype <NA> Lake_Bomoseen         N       <NA>          100       SL   <NA>     BS4 BS4_SL
# 3  FB1_Hap1_2N_vcf  FB1.2 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         Y       <NA>           75 ShortDNA   <NA>   FB1.2  FB1.2
# 4  FB1_Hap2_2N_vcf  FB1.1 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         Y       <NA>           25 ShortDNA   <NA>   FB1.1  FB1.1
# 5 FB11_Hap1_2N_vcf   FB11 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         N       <NA>          100 ShortDNA   <NA>    FB11   FB11
# 6 FB12_Hap1_2N_vcf   FB12 Lake_Memphremagog Normal_Haplotype <NA>     Fitch_bay         N       <NA>          100 ShortDNA   <NA>    FB12   FB12
# ....

# join the tree data with metadata
x <- left_join(as_tibble(tree), tree_labels, by = c("label" = "name"))  %>% mutate(name=label)

# format SNP data
snps <- snps %>% remove_rownames() %>% column_to_rownames(var="X")
gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
  x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]
snp_data <- snp_data 

# Build tree
p <- ggtree(tree) %<+% x  + 
  geom_tippoint(aes(color = TumorP, shape=libtype),size=2)  + coord_cartesian(clip = 'off') + geom_tiplab(aes(label=call, color=TumorP, fontface="bold"),size=2.4, align=TRUE, linesize=.1) +
  geom_nodelab( aes(label=label, subset=label > 79), hjust = 1.44, nudge_y = 0.65, size = 3) +
  geom_nodelab( aes(label=label, subset=label ==100 ), hjust = 1.44, nudge_y = 0.65, size = 3) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(10, 200, 10, 10))  + 
  theme(axis.text.x = element_text(face="bold", size = 10)) +
  scale_color_manual(values=c("#3A9AB2","#3A9AB2","orange1","#F11B00","#3A9AB2","#3A9AB2","orange1")) + theme(legend.position="none") 
p

# add snp panel to tree
p1 <- p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = TumorP), shape = 16) + 
  theme_tree2(plot.margin=margin(5, 200, 5, 5),)   + 
  theme(axis.text.x = element_text(face="bold", size = 10))

p1


######################################
# Generate VCF mitochondrial phylogenetic trees in R - short reads only
######################################

rm(list = ls())
options(java.parameters="-Xmx50G")

library(fastreeR)

library(rJava)
.jinit()
num_gigs_ram_available = .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory") / 1e9 
paste("You have ", round(num_gigs_ram_available, 2), "GB memory available.", sep = "")

vcfFile <- "path_to/short_TN_mito_snps-and-indels_clean_0.9_10x.vcf.gz"
 
vcftree<- vcf2tree(
    vcfFile,
    threads = 10,
    ignoreMissing = FALSE,
    onlyHets = FALSE,
    ignoreHets = FALSE 
  )
  
vcftree
