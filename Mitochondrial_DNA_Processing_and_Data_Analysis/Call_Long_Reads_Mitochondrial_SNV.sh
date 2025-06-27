####################
# Call SNVs for long reads 
####################


# assign variables and make directories

sample=<sample name>
mito_ref=<path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa
dir_path=<path to long read mitochondrial mapping output>
raw_sample_path=${merged_bam}/bam/All_${sample}.bam

mkdir -p ${dir_path}/${sample}

#####
# Map fastq reads to Mitochondria using dorado
# Use dorado to trim addapter -just in case- and then map to the mito reference
#####

dorado trim ${raw_sample_path} > ${dir_path}/${sample}/${sample}_trimmed.bam
dorado aligner ${mito_ref} ${dir_path}/${sample}/${sample}_trimmed.bam  > ${dir_path}/${sample}/${sample}_aligned.bam

#######
# Use samtools to subset mapped reads, then sort, update sample identifiers, and index the bam file
#######

samtools view -F 0x04  -b ${dir_path}/${sample}/${sample}_aligned.bam > ${dir_path}/${sample}/${sample}_mito.bam
samtools sort ${dir_path}/${sample}/${sample}_mito.bam -o ${dir_path}/${sample}/${sample}_mito_sorted.bam
samtools addreplacerg -r 'ID:UVM' -r 'PL:ONT' -r 'SM:'${sample} -o ${dir_path}/${sample}/${sample}_mito_only.bam ${dir_path}/${sample}/${sample}_mito_sorted.bam
samtools index ${dir_path}/${sample}/${sample}_mito_only.bam


#######
# remove superfluous files
#######

rm ${dir_path}/${sample}/${sample}_trimmed.bam
rm ${dir_path}/${sample}/${sample}_mito.bam
rm ${dir_path}/${sample}/${sample}_aligned.bam
rm ${dir_path}/${sample}/${sample}_mito_sorted.bam

#######
# Call Variants with GATK's Haplotype Caller  
#######
gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=16" HaplotypeCaller -R ${mito_ref} -I ${dir_path}/${sample}/${sample}_mito_only.bam  -O ${dir_path}/${sample}/${sample}_mito_only.vcf -ERC GVCF -ploidy 2



