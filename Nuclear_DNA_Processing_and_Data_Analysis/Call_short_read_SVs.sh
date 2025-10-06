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

#####
# Using samtools downsample to 130M reads, sort, index and give final counts for bam files 
#####

samtools view -b -h -@ 5 -s ${FRACTION} ${input} -o ${output}/${sample}_down.bam
samtools sort -o ${output}/${sample}_down_sorted.bam -@ 5 ${output}/${sample}_down.bam
samtools index -b -@ 10 ${output}/${sample}_down_sorted.bam

samtools view -@ 5 -c ${output}/${sample}_down_sorted.bam >> ${output}/final.counts.txt
echo ${sample} >> ${output}/final.counts.txt
rm ${output}/${sample}_down.bam

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
# Remove sites from VCF with low quality calls and less than 10X or fewer than 90% of the samples will reads for that site
########

vcftools --vcf  ${output}/merged_130M.vcf --remove-filtered-all --max-missing 0.9 --minDP 10 --recode --out ${output}/merged_130M.PASS_.9_10x


########
# Split VCF by INS, INV, DEL, DUP
########

bcftools view -i 'INFO/SVTYPE="INS"' ${output}/merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.all.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${output}/merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.all.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${output}/merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.all.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${output}/merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.all.del.vcf




########
# Subset VCF for only paired Tumor Normal (TN) samples
########

bcftools view -s FB8N,HC2N,HC4N,JR9N,MA1N,SB14N,SC6N,SC7N,SML1N,SML4N,FB8T,JR9T,HC4T,HC2T,MA1T,SML1T,SML4T,SC6T,SC7T,SB14T ${path}/${output}merged_130M.PASS_.9_10x.recode.vcf > ${output}/merged_130M.PASS_.9_10x.TN.vcf


########
# Split TN only VCF by INS, INV, DEL, DUP
########

bcftools view -i 'INFO/SVTYPE="INS"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${output}/merged_130M.PASS_.9_10x.TN.vcf > ${output}/merged_130M.PASS_.9_10x.TN.del.vcf


########
# Get Allele Depth for VCFs
########


bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.TN.ins.vcf > ${output}/merged_130M.PASS_.9_10x.TN.ins.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.TN.inv.vcf > ${output}/merged_130M.PASS_.9_10x.TN.inv.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.TN.dup.vcf > ${output}/merged_130M.PASS_.9_10x.TN.dup.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.TN.del.vcf > ${output}/merged_130M.PASS_.9_10x.TN.del.vcf.AD.txt

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.all.ins.vcf > ${output}/merged_130M.PASS_.9_10x.all.ins.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.all.inv.vcf > ${output}/merged_130M.PASS_.9_10x.all.inv.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.all.dup.vcf > ${output}/merged_130M.PASS_.9_10x.all.dup.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${output}/merged_130M.PASS_.9_10x.all.del.vcf > ${output}/merged_130M.PASS_.9_10x.all.del.vcf.AD.txt


