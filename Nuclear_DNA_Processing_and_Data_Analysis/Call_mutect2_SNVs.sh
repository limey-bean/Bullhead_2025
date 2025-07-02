######################
# Use somatic SNV caller Mutect2 to identify tumor specific SNVs
# use to compare with GATK diploid haplotyping and genotyping workflow
######################

sample=<paired_sample_name> # no T or N

path=/path_to_main_output_dir/
workdir=${path}/short_read_genome_mapping/mutect2_with_germline_N_som_T/${sample}/

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
  bcftools view -s ${s}T ${path}${s}/${s}T_mutect2_tumor_Filtered.vcf.gz -o ${path}${s}/${s}T_mutect2_Filtered.vcf.gz
  bcftools view --max-alleles 2 ${path}${s}/${s}T_mutect2_Filtered.vcf.gz > ${path}${s}/${s}T_mutect2_Filtered_no_poly.vcf
  bcftools view -f 'PASS' ${path}${s}/${s}T_mutect2_Filtered_no_poly.vcf  > ${path}${s}/${s}T_mutect2_PASS.vcf
  bgzip ${path}/${s}/${s}T_mutect2_PASS.vcf
  tabix ${path}/${s}/${s}T_mutect2_PASS.vcf.gz
done

#########
# Extract list of pass variant sites to filter all paired brain and tumor samples
########

bcftools merge ${workdir}/${sample}T_mutect2_tumor_Filtered.vcf.gz -Oz -o ten_mutect_TN_merged_files.vcf.gz


reheader.txt 
FB8N	FB8Nm
HC2N	HC2Nm
HC4N	HC4Nm
JR9N	JR9Nm
SB14N	SB14Nm
SC6N	SC6Nm
SC7N	SC7Nm
SML1N	SML1Nm
SML4N	SML4Nm

bcftools reheader --samples reheader.txt -o ten_mutect_TN_merged_reheader.vcf.gz ten_mutect_TN_merged_files.vcf.gz
tabix ten_mutect_TN_merged_reheader.vcf.gz

bcftools view ten_mutect_TN_merged_reheader.vcf.gz --regions HL4_0001,HL4_0002,HL4_0003,HL4_0004,HL4_0005,HL4_0006,HL4_0007,HL4_0008,HL4_0009,HL4_0010,HL4_0011,HL4_0012,HL4_0013,HL4_0014,HL4_0015,HL4_0016,HL4_0017,HL4_0018,HL4_0019,HL4_0020,HL4_0021,HL4_0022,HL4_0023,HL4_0024,HL4_0025,HL4_0026,HL4_0027,HL4_0028,HL4_0029,HL4_0030,HL4_0031,HL4_0032,HL4_0033,HL4_0034,HL4_0035,HL4_0036,HL4_0037,HL4_0038,HL4_0039,HL4_0040,HL4_0041,HL4_0042,HL4_0043,HL4_0044,HL4_0045,HL4_0046,HL4_0047,HL4_0048,HL4_0049,HL4_0050,HL4_0051,HL4_0052,HL4_0053,HL4_0054,HL4_0055,HL4_0056,HL4_0057,HL4_0058,HL4_0059,HL4_0060,HL4_0061,HL4_0062,HL4_0063,HL4_0064,HL4_0065,HL4_0066,HL4_0067,HL4_0068,HL4_0069,HL4_0070,HL4_0071,HL4_0072,HL4_0073,HL4_0074,HL4_0075,HL4_0076,HL4_0077,HL4_0078,HL4_0079,HL4_0080,HL4_0081,HL4_0082,HL4_0083,HL4_0084,HL4_0085,HL4_0086,HL4_0087,HL4_0088,HL4_0089,HL4_0090,HL4_0091,HL4_0092,HL4_0093,HL4_0094,HL4_0095,HL4_0096,HL4_0097,HL4_0098,HL4_0099,HL4_0100 > ten_mutect_TN_merged_reheader1-100.vcf

bgzip ten_mutect_TN_merged_reheader1-100.vcf
tabix ten_mutect_TN_merged_reheader1-100.vcf.gz 

