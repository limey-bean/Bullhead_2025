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
            -n ${dir_path}/${sample}/Long_reads/HL4_no_mito.fastq.gz \
            -s1 <path to directory that contains the concatonated short_reads for the *_R1_001.fastq.gz and *_R2_001.fastq.gz>  \
            -phv ~/PHVindexes\ # centrifuge database for removing contamination
            -b ~/actinopterygii_odb10\ # Busco database for fish
            -qg ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.gff \ # genome of closely related species - for additional scaffolding - we did not use results because A melas has a different number of chromosomes thatn A. nebulosa
            -qf ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.fna \
            -r ${dir_path}/${sample}/HL4_genome_no_mitochondria.fasta


######################################
# Genome Annotation using Funnanotate
# 

#add genome to the following directory with the genome name specified below ${genome}.fasta

cd ${dir_path}/${sample}/

genomedir=${dir_path}/${sample}/
genome=HL4_genome_no_mitochondria
RNAseq_path1R1=<path to clean RNAseq data sample 1 R1>
RNAseq_path1R2=<path to clean RNAseq data sample 1 R2>
RNAseq_path2R1=<path to clean RNAseq data sample 2 R1>
RNAseq_path2R2=<path to clean RNAseq data sample 2 R2>
....

# I used RNAseq samples vt1_03, vt2_10, vt1_07, vt2_23, vt1_02, vt2_06, vt2_32, vt2_29R, vt2_05

#######################
# clean up names in genome file
#######################

# remove artifacts in contig file
sed 's/JAAGNN01000//g' ${genomedir}${genome}.fasta > ${genomedir}${genome}.1.fasta
sed 's/_polished//g' ${genomedir}${genome}.1.fasta > ${genomedir}${genome}.2.fasta

# use funannotate's cleanup function  remove all contigs < 1000 bp
funannotate \
   clean -i ${genomedir}${genome}.2.fasta --minlen 1000 -o ${genomedir}${genome}.cleaned.fa

# sort contigs by length
seqkit sort -l --reverse --two-pass ${genomedir}${genome}.cleaned.fa > ${genomedir}${genome}.cleaned.sorted.fa

# rename reads
seqkit replace -p .+ -r HL4_{nr} --nr-width 4 ${genomedir}${genome}.cleaned.sorted.fa -o ${genomedir}${genome}.cleaned.sorted.renamed.fa

# remove superfluous files
rm ${genomedir}${genome}.1.fasta
rm ${genomedir}${genome}.2.fasta
rm ${genomedir}${genome}.cleaned.fa
rm ${genomedir}${genome}.cleaned.sorted.fa


#######################
# mask assembly....
#######################

funannotate \
  mask -i ${genomedir}${genome}.cleaned.sorted.renamed.fa  --cpus 20 -o ${genomedir}${genome}_masked_assembly.fa

######################
# train gene model with rnaseq data
######################


 funannotate \
    train -i ${genomedir}${genome}_masked_assembly.fa -o ${genome}Fun_out \
    --left ${RNAseq_path1R1} ${RNAseq_path2R1} ... \
    --right ${RNAseq_path1R2} ${RNAseq_path2R2} ... \
    --stranded RF --jaccard_clip --cpus 20

######################
# predict genes
######################

 funannotate \
    predict -i ${genomedir}${genome}_masked_assembly.fa  -o ${genome}Fun_out -s "HL4_no_scaffold_masked_assembly" \
    --cpus 20 --busco_seed_species zebrafish  --busco_db ${path_toFunannotate_db}/actinopterygii --organism  other --optimize_augustus --database ${path_toFunannotate_db}

######################
# update model
######################

funannotate \
    update -i ${genome}Fun_out --cpus 20

######################
# annotate final 
######################

# I ran the ${genome}Fun_out/predict_results/HL4_no_scaffold_masked_assembly.proteins.fa through the eggnog_mapper annotation webserver (http://eggnog-mapper.embl.de) and added the output to:
mkdir ${genome}Fun_out/eggnog_mapper_out/

funannotate \
     annotate -i ${genome}Fun_out --cpus 20 --eggnog ${genome}Fun_out/eggnog_mapper_out/out.emapper.annotations \
     --isolate HL4 --busco_db ${path_toFunannotate_db}/actinopterygii -d ${path_toFunannotate_db} --species "Ameiurus nebulosus"


#######
# Index ref genome with samtools and gatk
######

samtools faidx <path to>/A_nebulosus_HL4.scaffolds.fa
samtools index <path to>/A_nebulosus_HL4.scaffolds.fa 
gatk CreateSequenceDictionary -R <path to>/A_nebulosus_HL4.scaffolds.fa
