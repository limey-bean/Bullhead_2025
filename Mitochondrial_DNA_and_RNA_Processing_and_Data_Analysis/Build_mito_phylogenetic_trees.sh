######################################
# Generate phased mitochondrial phylogenetic trees
######################################

##########
# first concatonate the fasta files or interest into one mutlifasta file
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
