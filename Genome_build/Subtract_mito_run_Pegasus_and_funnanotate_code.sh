##########
# Remove Mitochondrial reads from the Long reads used to generate the Genome
##########

#####
# assign variables and make new directories
#####
sample=<sample name>

mito_ref=<path to the mitochondrial reference genome>
dir_path=<path to longreads minus mitochondria output directory>
sample_long_reads_bam=<path to the long reads bam file for the sample>

mkdir -p ${dir_path}/${sample}/cleaned_fastq

#####
# use dorado to trim reads, just in case...
#####

dorado trim ${sample_long_reads_bam} > ${dir_path}/${sample}/${sample}_trimmed.bam

#####
# use dorado to align treimmed reads to the mitochondrial reference genome
#####

dorado aligner ${mito_ref} ${dir_path}/${sample}/${sample}_trimmed.bam  > ${dir_path}/${sample}/${sample}_aligned.bam

#####
# remove mitochondrial samples from the Bam file
#####

samtools view -f 4  -b ${dir_path}/${sample}/${sample}_aligned.bam > ${dir_path}/${sample}/${sample}_no_mito.bam

#####
# Sort reads in mitochondria free Bam file
#####

samtools sort ${dir_path}/${sample}/${sample}_no_mito.bam -o ${dir_path}/${sample}/${sample}_no_mito_sorted.bam

#####
# Reaname sample feilds in bam file then index bam file
#####

samtools addreplacerg -r 'ID:UVM' -r 'PL:Illumina' -r 'SM:'${sample_name} -o ${dir_path}/${sample}/${sample}_mito_only.bam ${dir_path}/${sample}/${sample}_no_mito_sorted.bam

samtools index ${dir_path}/${sample}/${sample}_no_mito_only.bam

#####
# Remove superfluous files
#####

rm ${dir_path}/${sample}/${sample}_no_mito.bam
rm ${dir_path}/${sample}/${sample}_aligned.bam
rm ${dir_path}/${sample}/${sample}_no_mito_sorted.bam

#####
# Convert the bam file to a fastq file
#####

bedtools bamtofastq -i ${dir_path}/${sample}/${sample}_no_mito_only.bam -fq ${dir_path}/${sample}/cleaned_fastq/HL4_no_mito.fastq.gz



######################################
# Hybrid Genome Assembly Using PEGASUS
# Pegasus requires singularity or apptainer
# module load singularity
# command to run Pegasus. For more details see: https://github.com/jaxlub/PEGASUS

~/PEGASUS/pegasus.sh \
            -n ~/HL4/Long_reads/HL4_no_mito.fastq.gz \
            -s1 <path to concatonated short_reads1>  \
            -s2 <path to concatonated short_reads2> \
            -phv ~/PHVindexes\ # centrifuge database for removing contamination
            -b ~/actinopterygii_odb10\ # Busco database for fish
            -qg ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.gff \ # genome of closely related species - for additional scaffolding - we did not use results because A melas has a different number of chromosomes thatn A. nebulosa
            -qf ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.fna \
            -r ~/HL4RefGen.fasta


######################################
# Genome Annotation using Funnanotate

