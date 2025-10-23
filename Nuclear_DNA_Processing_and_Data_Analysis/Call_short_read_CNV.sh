###########################
# Call short reads Nuclear CNVs
###########################

You need to make a mappability file before running this script.  

The following is Copied from https://github.com/dellytools/delly, you will need to modify it for the brown bullhead genome:

dicey chop sacCer3.fa
bwa index sacCer3.fa
bwa mem sacCer3.fa read1.fq.gz read2.fq.gz | samtools sort -@ 8 -o srt.bam -
samtools index srt.bam 
dicey mappability2 srt.bam 
gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 


########
# set paths and variables
########

ref=/gpfs1/home/e/g/eguswa/scratch/bullhead/Funannotate/HL4/no_mitochondria/HL4_no_mitochondriaFun_out/annotate_results/Ictalurus_punctatus_HL4.scaffolds.fa
map=/gpfs1/home/e/g/eguswa/scratch/bullhead/Funannotate/HL4/no_mitochondria/HL4_no_mitochondriaFun_out/annotate_results/map.fa.gz

# we downsampled to 130M reads
sample==<name of sample>
ref=/path_to/HL4.scaffolds.fa
map=/path_to_map_file/map.fa.gz
input=/path_to/downsampled_130M/${sample}_down_sorted.bam
output=/path_to/CNV_down130M/

mkdir -p ${output}
cd ${output}

########
# Run Delly2
########

delly cnv -u -g ${ref} --mappability ${map}  -c ${sample}.cov.gz -o ${sample}_cnv.bcf ${input}

