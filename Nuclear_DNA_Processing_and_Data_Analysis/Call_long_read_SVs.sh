###########################
# Call long reads Nuclear SVs using Sniffles
###########################

########
# set paths and variables
########

sample=<sample name>

out=/path_to_output_for/sniffles/

datadir=/path_to/map_long/
ref_genome=/path_to/HL4.scaffolds.fa

mkdir ${out}

#######
# Run Sniffles on each sample individually 
#######

sniffles --input ${datadir}/${sample}/${sample}_aligned_sorted_filtered.bam --snf ${out}/${sample}.snf


#######
# Merge Sniffles output from each individual  and make a vcf file
#######

sniffles --input ${out}/${sample1}.snf ${out}/${sample2}.snf ... ${out}/${sampleN}.snf --vcf sniffles_all.vcf


########
# Remove sites from VCF with less than 10X or fewer than 90% of the samples will reads for that site
########

vcftools --gzvcf  ${out}/sniffles_all.vcf --max-missing 0.9 --minDP 10 --recode --out ${out}/sniffles_all_.9_10x


########
# Split VCF by INS, INV, DEL, DUP
########

bcftools view -i 'INFO/SVTYPE="INS"' ${out}/sniffles_all_.9_10x.recode.vcf > ${out}/sniffles_all_.9_10x.all.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${out}/sniffles_all_.9_10x.recode.vcf > ${out}/sniffles_all_.9_10x.all.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${out}/sniffles_all_.9_10x.recode.vcf > ${out}/sniffles_all_.9_10x.all.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${out}/sniffles_all_.9_10x.recode.vcf > ${out}/sniffles_all_.9_10x.all.del.vcf




########
# Subset VCF for only paired Tumor Normal (TN) samples
########

bcftools view -s FB8N,HC2N,HC4N,JR9N,MA1N,SB14N,SC6N,SC7N,SML1N,SML4N,FB8T,JR9T,HC4T,HC2T,MA1T,SML1T,SML4T,SC6T,SC7T,SB14T ${path}/${out}sniffles_all_.9_10x.recode.vcf > ${out}/sniffles_all_.9_10x.TN.vcf


########
# Split TN only VCF by INS, INV, DEL, DUP
########

bcftools view -i 'INFO/SVTYPE="INS"' ${out}/sniffles_all_.9_10x.TN.vcf > ${out}/sniffles_all_.9_10x.TN.ins.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${out}/sniffles_all_.9_10x.TN.vcf > ${out}/sniffles_all_.9_10x.TN.inv.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${out}/sniffles_all_.9_10x.TN.vcf > ${out}/sniffles_all_.9_10x.TN.dup.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${out}/sniffles_all_.9_10x.TN.vcf > ${out}/sniffles_all_.9_10x.TN.del.vcf


########
# Get Allele Depth for VCFs
########


bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.TN.ins.vcf > ${out}/sniffles_all_.9_10x.TN.ins.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.TN.inv.vcf > ${out}/sniffles_all_.9_10x.TN.inv.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.TN.dup.vcf > ${out}/sniffles_all_.9_10x.TN.dup.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.TN.del.vcf > ${out}/sniffles_all_.9_10x.TN.del.vcf.AD.txt

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.all.ins.vcf > ${out}/sniffles_all_.9_10x.all.ins.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.all.inv.vcf > ${out}/sniffles_all_.9_10x.all.inv.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.all.dup.vcf > ${out}/sniffles_all_.9_10x.all.dup.vcf.AD.txt
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${out}/sniffles_all_.9_10x.all.del.vcf > ${out}/sniffles_all_.9_10x.all.del.vcf.AD.txt




