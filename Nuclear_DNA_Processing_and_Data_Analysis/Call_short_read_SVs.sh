####################
# Call Structural Variants in short reads using Delly
# down sampled to 130M reads to include all paired samples and minimize effect of detection using variable number of reads.
####################


sample=<sample_name>
ref_genome=/path_to/HL4.scaffolds.fa


################### Run as a loop of array for all samples #############################

#######
# set variables make directories
#######

input=/path_to/${sample}/${sample}_mapped.genome_only.bam
output=/path_to/delly/downsampled_130M/

mkdir -p ${output}
cd ${output}

#####
# Get fraction of bam file that catches ~130 M reads per sample
#####

echo ${sample} 
COUNT=$(${samtools} view -@ 10 -c ${input})
echo ${COUNT}
FRACTION=$(awk -v COUNT=${COUNT} "BEGIN {print 130000000 / COUNT}") # use awk to calculate fraction
echo ${FRACTION}

## for loop to deal with samples with slightly less depth

if [[ ${FRACTION} > 1 ]]; then
  
	cp ${input}  ${output}/${sample}_down_sorted.bam
        cp ${input}.bai  ${output}/${sample}_down_sorted.bam.bai
       echo ${sample} >> ${output}/final.counts
  ${samtools} view -@ 5 -c ${output}/${sample}_down_sorted.bam >> ${output}/final.counts

else

  echo "my_variable is not greater than 1"

#####
# Using samtools downsample to 130M reads, sort, index and give final counts for bam files 
#####

  ${samtools} view -b -h -@ 5 -s ${FRACTION} ${input} -o ${output}/${sample}_down.bam
 
  ${samtools} sort -o ${output}/${sample}_down_sorted.bam -@ 5 ${output}/${sample}_down.bam

  ${samtools} index -b -@ 10 ${output}/${sample}_down_sorted.bam

  echo ${sample} >> ${output}/final.counts

  ${samtools} view -@ 5 -c ${output}/${sample}_down_sorted.bam >> ${output}/final.counts

  rm ${output}/${sample}_down.bam

fi


########
# Call SVs with Delly
########

delly call -o ${output}/${sample}.delly.bcf -q 30 -g ${ref} ${output}/${sample}_down_sorted.bam

########
# Merge SVs sites with Delly merge
###########

delly merge -o sites.bcf *.bcf

########
# Call SVs with Delly PASS sites
########

delly call -v ${output}sites.bcf -o ${output}/${sample}.delly.PASS.bcf -q 30 -g ${ref} ${output}/${sample}_down_sorted.bam


########
# Merge all sample BCFs with bcftools
########

bcftools merge ${output}/*.delly.PASS.bcf -m id -Oz -o ${output}/merged_130M_PASS.bcf

########
# Convert .bcf to vcf with bcftools
########

bcftools view ${output}/merged_130M.PASS.bcf  > ${output}/merged_130M.PASS.vcf

########
# Split VCF by INS, INV, DEL, DUP
########

bcftools view -i 'INFO/SVTYPE="INS"' ${output}/merged_130M.PASS.vcf > ${output}/merged_130M.PASS.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${output}/merged_130M.PASS.vcf > ${output}/merged_130M.PASS.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${output}/merged_130M.PASS.vcf > ${output}/merged_130M.PASS.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${output}/merged_130M.PASS.vcf > ${output}/merged_130M.PASS.del.vcf



########
# Remove sites from VCF with low quality calls and less than 10X or fewer than 90% of the samples will reads for that site
########

vcftools --vcf  ${output}/merged_130M.PASS.ins.vcf --remove-filtered-all --max-missing 0.9 --recode --out ${output}/merged_130M.PASS.ins_.9
vcftools --vcf  ${output}/merged_130M.PASS.inv.vcf --remove-filtered-all --max-missing 0.9 --recode --out ${output}/merged_130M.PASS.inv_.9
vcftools --vcf  ${output}/merged_130M.PASS.dup.vcf --remove-filtered-all --max-missing 0.9 --recode --out ${output}/merged_130M.PASS.dup_.9
vcftools --vcf  ${output}/merged_130M.PASS.del.vcf --remove-filtered-all --max-missing 0.9 --recode --out ${output}/merged_130M.PASS.del_.9



########
# if needed Subset VCF for only paired Tumor Normal (TN) samples
########

bcftools view -S <path_to_sample_list.txt>  ${path}/${output}merged_130M.PASS_.9.recode.vcf > ${output}/merged_130M.PASS_.9_TN.vcf


########
# Get Allele Depth for VCFs
########


bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%RR,%RV]\n' ${output}/merged_130M.PASS.ins_.9.TN.vcf > ${output}/merged_130M.PASS.ins_.9.TN.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%RR,%RV]\n' ${output}/merged_130M.PASS.inv_.9.TN.vcf > ${output}/merged_130M.PASS.inv_.9.TN.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%RR,%RV]\n' ${output}/merged_130M.PASS.dup_.9.TN.vcf > ${output}/merged_130M.PASS.dup_.9.TN.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%RR,%RV]\n' ${output}/merged_130M.PASS.del_.9.TN.vcf > ${output}/merged_130M.PASS.del_.9.TN.AD.txt


