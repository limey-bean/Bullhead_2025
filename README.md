# Bullhead 2025 Data Analysis
A repository for scripts used to generate the 2025 Brown Bullhead data analysis

**Note:** I ran these analyses on the High Performance Compute Cluster at the University of Vermont’s Vermont Advanced Computing Center (VACC; RRID:SCR_017762).  I do not include the slurm commands or the arrays used to process multiple samples simultaneously. Presented in this repository are the base commands used to run software or modify outputs from software.

## Raw Data Processing
Contains scripts to process sequences generated for
1. Short Read RNA and DNA
2. Long Read DNA

## Genome Build
Contains scripts to build the 
1. Mitochondrial Genome
2. Draft Nuclear Genome:
   * hybrid assembly with PEGASUS
   * annotation funnanotate 

## Mitochondrial DNA Processing and Data Analysis
### Mitochondrial SNV Calling
Phylogenetic Tree Based on Genotype
### Mitogenome resonstruction
Phased Mitogenome Phylogenetic tree

## Nuclear DNA Processing and Data Analysis
### SNV Calling
HaplotypeCaller
Mutect2
### SV Calling
Short Read DNA
Long Read DNA

## Human TCGA Scripts
