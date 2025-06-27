####################
# Call SNVs for short reads DNA and RNA
####################

#####
# assign variables and make directories
#####


sample=<sample name>
mito_ref=<path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa
dir_path=<path to long read mitochondrial mapping output>
raw_sample_path=${merged_bam}/bam/All_${sample}.bam
workdir=${dir_path}/Mito_Variant_Calling/

mkdir -p ${workdir}/comb_and_genotype_VFC_long_dna/



sample_path=/gpfs1/home/e/g/eguswa/scratch/bullhead/short_reads/fastp/
path=/gpfs1/home/e/g/eguswa/scratch/bullhead/is_it_transmissible_cancer/
workdir=${path}/mitochondria/short_mito_mapping/${sample_name}/

mkdir -p ${workdir}
cd ${workdir}

module load singularity
singularity_path=/gpfs1/home/e/g/eguswa/scratch/Containers/gatk_latest.sif
ref_genome=/gpfs1/home/e/g/eguswa/scratch/bullhead/getorganelle/HL4_mt_out/stand/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa

my_job_header

fin=${sample_path}/${sample_name}_clean_R1_001.fastq.gz
rin=${sample_path}/${sample_name}_clean_R2_001.fastq.gz
map_name=${workdir}/${sample_name}_mapped

#######
# Index ref genome
######

if [[ -f /gpfs1/home/e/g/eguswa/scratch/bullhead/getorganelle/HL4_mt_out/stand/HL4_ann.final-mtgrasp_v1.1.8-assembly.dict ]]; then
   echo "genome appears to be indexed"
else
   bwa index ${ref_genome}
   singularity exec  ${singularity_path} \
     samtools faidx ${ref_genome}
   singularity exec  ${singularity_path} \
     gatk CreateSequenceDictionary -R ${ref_genome}
fi

if [[ -f ${ref_genome}.fai ]]; then
   echo "genome appears to be indexed (fai)"
else
  singularity exec  ${singularity_path} \
   samtools index ${ref_genome}
fi


########
# downsample
########

#if [[ -f ${fin} ]]; then
#  echo "forward reads have been downsampled"
#else
#  seqkit head -n 30000000 ${sample_path}/${sample_name}_all_R1_001.fastq.gz  > ${fin}
#fi

#if [[ -f ${rin} ]]; then
#  echo "reverse reads have been downsampled"
#else
#  seqkit head -n 30000000 ${sample_path}/${sample_name}_all_R2_001.fastq.gz  > ${rin}
#fi




########
# map to HL4 mitochondria
########

if [[ -f ${map_name}.sam ]]; then
  echo "Reads were mapped"
else
  minimap2 -t 1 -ax sr ${ref_genome} ${fin} ${rin} > ${map_name}.sam
fi

if [[ -f ${map_name}.bam ]]; then
  echo "Reads were converted to bam"
else
  singularity exec  ${singularity_path} \
  samtools view -S -b ${map_name}.sam -o ${map_name}.bam
fi


if [[ -f ${map_name}.bam ]]; then
  rm ${map_name}.sam
else
  echo "Conversion from sam to bam failed."
  exit
fi


if [[ -f ${map_name}.sorted.bam ]]; then
  echo "Bam Reads were sorted"
else
  singularity exec  ${singularity_path} samtools sort ${map_name}.bam -o ${map_name}.sorted.bam
fi


if [[ -f ${map_name}.sorted.bam ]]; then
  rm ${map_name}.bam
else
  echo "Sorting bam failed."
  exit
fi

if [[ -f ${map_name}.sorted.un.bam ]]; then
  echo "Bam Reads were renamed"
else
  singularity exec  ${singularity_path} samtools addreplacerg -r 'ID:UVM' -r 'PL:Illumina' -r 'SM:'${sample_name} -o ${map_name}.sorted.un.bam ${map_name}.sorted.bam

fi


if [[ -f ${map_name}.sorted.un.bam ]]; then
  rm ${map_name}.sorted.bam
else
  echo "Renaming bam feilds failed."
  exit
fi

if [[ -f ${map_name}.mito_only.bam ]]; then
  echo "Only Mito alignments retaind in Bam file"
else
  singularity exec  ${singularity_path} samtools view -F 0x04  -b ${map_name}.sorted.un.bam > ${map_name}.mito_only.bam
  samtools index ${map_name}.mito_only.bam
fi


if [[ -f ${map_name}.mito_only.bam ]]; then
  rm ${map_name}.sorted.un.bam
else
  echo "Renaming bam feilds failed."
  exit
fi


if [[ -f ${sample_name}_gatk_N1.vcf ]]; then
  echo "VCF N=1 was generated"
else
  singularity exec  ${singularity_path} gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=16" HaplotypeCaller -R ${ref_genome} -I ${map_name}.mito_only.bam  -O ./${sample_name}_gatk_N1.vcf -ERC GVCF -ploidy 1
fi


if [[ -f ${sample_name}_gatk_N2.vcf ]]; then
  echo "VCF N=2 was generated"
else
  singularity exec  ${singularity_path} gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=16" HaplotypeCaller -R ${ref_genome} -I ${map_name}.mito_only.bam  -O ./${sample_name}_gatk_N2.vcf -ERC GVCF -ploidy 2
fi

