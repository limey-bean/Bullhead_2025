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


tree <- read.tree(file="path_to/mito_reads.align.fasta.contree")



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
