####################
# Map fastq reads to Mitochondria
####################

#####
# Long Reads
#####

sample=<sample name>
mito_ref=/gpfs1/home/e/g/eguswa/scratch/bullhead/getorganelle/HL4_mt_out/stand/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa

my_job_header

dir_path=/gpfs1/home/e/g/eguswa/scratch/bullhead/is_it_transmissible_cancer/mitochondria/long_mito_mapping/
mkdir -p ${dir_path}/${sample}

raw_sample_path=/netfiles02/mbsr/Projects/VFWBullhead/raw_concatonated_long_reads-bam_and_fastq/longreads_merge_call/${sample}/${sample}_dorodo/bam/All_${sample}.bam


if [[ -f ${dir_path}/${sample}/${sample}_trimmed.bam ]]; then
  echo "Bam file trimmed of adapters"
else
  /gpfs1/home/e/g/eguswa/scratch/bullhead/longreads_merge_call/dorado-0.7.2-linux-x64/bin/dorado trim ${raw_sample_path} > ${dir_path}/${sample}/${sample}_trimmed.bam
fi

if [[ -f ${dir_path}/${sample}/${sample}_aligned.bam ]]; then
  echo "Bam file alighed to mito ref"
else
  /gpfs1/home/e/g/eguswa/scratch/bullhead/longreads_merge_call/dorado-0.7.2-linux-x64/bin/dorado aligner ${mito_ref} ${dir_path}/${sample}/${sample}_trimmed.bam  > ${dir_path}/${sample}/${sample}_aligned.bam
fi

if [[ -f ${dir_path}/${sample}/${sample}_mito.bam ]]; then
  echo "Only Mito alignments retiand in Bam file"
else
  singularity exec ${singularity_path} samtools view -F 0x04  -b ${dir_path}/${sample}/${sample}_aligned.bam > ${dir_path}/${sample}/${sample}_mito.bam
fi

if [[ -f ${dir_path}/${sample}/${sample}_mito_sorted.bam ]]; then
  echo "Bam Reads were sorted"
else
  singularity exec ${singularity_path} samtools sort ${dir_path}/${sample}/${sample}_mito.bam -o ${dir_path}/${sample}/${sample}_mito_sorted.bam
fi


if [[ -f ${dir_path}/${sample}/${sample}_mito_sorted.bam ]]; then
  rm ${dir_path}/${sample}/${sample}_trimmed.bam
#  rm ${dir_path}/${sample}/${sample}_mito.bam
  rm ${dir_path}/${sample}/${sample}_aligned.bam
else
  echo "Sorting bam failed."
  exit
fi


if [[ -f ${dir_path}/${sample}/${sample}_mito_only.bam ]]; then
  echo "Bam Reads were renamed"
else
  singularity exec ${singularity_path} samtools addreplacerg -r 'ID:UVM' -r 'PL:Illumina' -r 'SM:'${sample} -o ${dir_path}/${sample}/${sample}_mito_only.bam ${dir_path}/${sample}/${sample}_mito_sorted.bam

  singularity exec ${singularity_path} samtools index ${dir_path}/${sample}/${sample}_mito_only.bam

fi

if [[ -f ${dir_path}/${sample}/${sample}_mito_only.bam.bai ]]; then
  rm ${dir_path}/${sample}/${sample}_mito_sorted.bam
else
  echo "Indexing bam failed."
  exit
fi

if [[ -f ${sample_name}_gatk_N2.vcf ]]; then
  echo "VCF N=2 was generated"
else
  singularity exec  ${singularity_path} gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=16" HaplotypeCaller -R ${mito_ref} -I ${dir_path}/${sample}/${sample}_mito_only.bam  -O ${dir_path}/${sample}/${sample}_mito_only.vcf -ERC GVCF -ploidy 2
fi

