######################
# Use somatic SNV caller Mutect2 to identify tumor specific SNVs
# use to compare with GATK diploid haplotyping and genotyping workflow
######################

sample=<paired_sample_name> # no T or N

path=/path_to_main_output_dir/short_read_genome_mapping/mutect2_with_germline_N_som_T/
workdir=${path}/${sample}/

singularity_path='apptainer -q exec /gpfs1/home/e/g/eguswa/scratch/Containers/gatk_latest.sif'
ref_genome=/path_to/HL4.scaffolds.fa

mkdir -p ${workdir}

#########
# Use Mutect to call somatic variants
########

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=5" \
  Mutect2 \
    -R ${ref_genome} \
    -I /path_to/short_read_genome_mapping/${sample}T/${sample}T_mapped.genome_only.bam \
    -I /path_to/short_read_genome_mapping/${sample}N/${sample}N_mapped.genome_only.bam \
    -normal ${sample}N --genotype-germline-sites true \
    --callable-depth 5 \
    -O ${workdir}${sample}T_mutect2_tumor_only.vcf.gz

#########
# Use FilterMutectCalls to filter somatic variants
########

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=5" \
  FilterMutectCalls \
    -R ${ref_genome} \
    -V ${workdir}/${sample}T_mutect2_tumor_only.vcf.gz \
    -O ${workdir}/${sample}T_mutect2_tumor_Filtered.vcf.gz

#########
# Use bcftools to find bialleleic pass variants for each tumor sample
########

samples=(
FB8
HC2
HC4
JR9
SB14
SC7
SML1
MA1
SML4
SC6
)

for s in ${samples[@]};
do
  echo ${s}
  bcftools view -s ${s}T ${workdir}/${sample}/${sample}T_mutect2_tumor_Filtered.vcf.gz -o ${workdir}/${sample}/${sample}T_mutect2_Filtered.vcf.gz
  bcftools view --max-alleles 2 ${workdir}/${sample}/${sample}T_mutect2_Filtered.vcf.gz > ${workdir}/${sample}/${sample}T_mutect2_Filtered_no_poly.vcf
  bcftools view -f 'PASS' ${workdir}/${sample}/${sample}T_mutect2_Filtered_no_poly.vcf  > ${workdir}/${sample}/${sample}T_mutect2_PASS.vcf
  bgzip ${workdir}/${sample}/${sample}T_mutect2_PASS.vcf
  tabix ${workdir}/${sample}/${sample}T_mutect2_PASS.vcf.gz
done

#########
# Merge Tumor pass vcfs, extract list of pass variant sites from tumor samples
########

bcftools merge ${workdir}/*T_mutect2_tumor_Filtered.vcf.gz -Oz -o Tumor_only_PASS.vcf.gz

bcftools query -H -f '%CHROM %POS\n' Tumor_only_PASS.vcf.gz > mutect2_pass_somatic_sites.txt


#########
# Filter TN from the GATK haplotyped and genotyped vcf to get the variants for the brain samples as well as the Tumor samples.
########

bcftools view -s FB8N,HC2N,HC4N,JR9N,MA1N,SB14N,SC6N,SC7N,SML1N,SML4N,FB8T,HC2T,HC4T,JR9T,MA1T,SB14T,SC6T,SC7T,SML1T,SML4T /path_to/HL4_0001-0100_snps_indels_clean_0.9_10x.vcf.gz -o ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels.vcf.gz


vcftools --gzvcf ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels.vcf.gz --max-missing 0.9 --minDP 10 --positions mutect2_pass_somatic_sites.txt --recode --out ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS

mv ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.recode.vcf.gz ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.vcf.gz

bgzip ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.vcf.gz
tabix ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.vcf.gz

######
# Use BCFtools to extract allele frequency
#####

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.vcf.gz > ${path}/Haplotyped_TN_pairs_HL4_0001-0100_snps_indels_90per_10x_PASS.vcf.AD.txt

