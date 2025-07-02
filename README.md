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
Contains scripts to
1. Call Mitochondrial SNV
   * Short read data
   * Long read data
2. Phase Mitochondria
3. Build phylogenetic trees using
  * SNVs
  * Mitochondrial reads

## Nuclear DNA Processing and Data Analysis
Contains scripts to
1. Call Nuclear SNVs
   * HaplotypeCaller
   * Mutect2
2. Call Nuclear SVs
   * Delly for Short Read DNA
   * Sniffles for Long Read DNA
3. Count tissue specific variants
4. Generate Phylogenetic trees from VCF

## Human TCGA Scripts
Contains
1. Scripts to parse Human melanoma SNVs.
2. Data files required to run the script.

## Viral and mitochondrial scripts
Contains scripts to run
1. nf-core mag
2. differential abundance analysis
