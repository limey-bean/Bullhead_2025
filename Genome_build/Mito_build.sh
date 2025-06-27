#######
# set varaibles for path to short reads used to make the genome and the sample name
#######

path=<path to>/concatonated_short_reads/ 
sample_name=HL4


#######
# Assemble mitogenome with getOrganelle
#######

get_organelle_from_reads.py -1 ${path}/${sample_name}_all_R1_001.fastq.gz \
			    -2 ${path}/${sample_name}_all_R2_001.fastq.gz \
           -R 10 -k 21,45,65,85,105 -F animal_mt -o <path to>/getorganelle/${sample_name}_mt_out --target-genome-size 17000  

#######
# Standardize origin and annotate mtgrasp
#######

mtgrasp_standardize.py -i <path to>/getorganelle/HL4_mt_out/HL4_complete_mitogenome.fasta -o ${sample_name}_stand  -c 2 -p ${sample_name}_ann -a -mp <path to>/miniforge3/bin

# the output name will resemble the following: HL4_ann.final-mtgrasp_v1.1.8-assembly.fa


#######
# Index ref genome with samtools and gatk
######

samtools faidx <path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa 
samtools index <path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa 
gatk CreateSequenceDictionary -R <path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa 

     
