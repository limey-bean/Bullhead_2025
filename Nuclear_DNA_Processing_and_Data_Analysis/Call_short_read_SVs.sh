####################
# Call Structural Variants in short reads using Delly
# down sampled to 130M reads to include all paired samples and minimize effect of detection using variable number of reads.
####################


sample=<sample_name>
ref_genome=/path_to/HL4.scaffolds.fa

################### Run as a loop of array for all samples #############################

input=/path_to/${sample}/${sample}_mapped.genome_only.bam
output=/path_to/delly/downsampled_130M/

mkdir -p ${output}
cd ${output}

echo ${sample} 
COUNT=$(${samtools} view -@ 10 -c ${input})
echo ${COUNT}
FRACTION=$(awk -v COUNT=${COUNT} "BEGIN {print 130000000 / COUNT}") # use awk to calculate fraction
echo ${FRACTION}

samtools view -b -h -@ 5 -s ${FRACTION} ${input} -o ${output}/${sample}_down.bam

samtools sort -o ${output}/${sample}_down_sorted.bam -@ 5 ${output}/${sample}_down.bam

samtools index -b -@ 10 ${output}/${sample}_down_sorted.bam

samtools view -@ 5 -c ${output}/${sample}_down_sorted.bam >> ${output}/final.counts.txt
echo ${sample} >> ${output}/final.counts.txt
rm ${output}/${sample}_down.bam

delly call -o ${output}/${sample}.delly.bcf -q 30 -g ${ref} ${output}/${sample}_down_sorted.bam

bcftools view ${output}/${sample}.delly.bcf  > ${output}/${sample}.delly.vcf

################### Run as a loop of array for all samples #############################

bcftools merge ${output}/*.delly.vcf -Oz -o ${output}/merged_130M.PASS.vcf

vcftools --gzvcf  ${output}/merged_130M.PASS.vcf --max-missing 0.9 --minDP 10 --recode --out ${output}/merged_130M.PASS_.9_10x

bcftools view -s FB8N,HC2N,HC4N,JR9N,MA1N,SB14N,SC6N,SC7N,SML1N,SML4N,FB8T,JR9T,HC4T,HC2T,MA1T,SML1T,SML4T,SC6T,SC7T,SB14T ${path}/${output}merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.TN.vcf


bcftools view -i 'INFO/SVTYPE="INS"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.del.vcf





