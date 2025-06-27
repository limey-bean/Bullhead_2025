####################
# Call SNVs for long reads 
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


#######
# Make a database of haplotyped calls GATK's GenomicsDBImport
#######

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" GenomicsDBImport \
-R ${ref_genome} \
-V ${dir_path}/BS4/BS4_mito_only.vcf \
-V ${dir_path}/FB8N/FB8N_mito_only.vcf \
-V ${dir_path}/FB8T/FB8T_mito_only.vcf \
-V ${dir_path}/HC14/HC14_mito_only.vcf \
-V ${dir_path}/HC2N/HC2N_mito_only.vcf \
-V ${dir_path}/HC2T/HC2T_mito_only.vcf \
-V ${dir_path}/HC4N/HC4N_mito_only.vcf \
-V ${dir_path}/HC4T/HC4T_mito_only.vcf \
-V ${dir_path}/HL3/HL3_mito_only.vcf \
-V ${dir_path}/HL4/HL4_mito_only.vcf \
-V ${dir_path}/JR11/JR11_mito_only.vcf \
-V ${dir_path}/JR19/JR19_mito_only.vcf \
-V ${dir_path}/JR9N/JR9N_mito_only.vcf \
-V ${dir_path}/JR9T/JR9T_mito_only.vcf \
-V ${dir_path}/LB1/LB1_mito_only.vcf \
-V ${dir_path}/LB2/LB2_mito_only.vcf \
-V ${dir_path}/MA18/MA18_mito_only.vcf \
-V ${dir_path}/MA1N/MA1N_mito_only.vcf \
-V ${dir_path}/MA1T/MA1T_mito_only.vcf \
-V ${dir_path}/SB14N/SB14N_mito_only.vcf \
-V ${dir_path}/SB14T/SB14T_mito_only.vcf \
-V ${dir_path}/SC1/SC1_mito_only.vcf \
-V ${dir_path}/SC6N/SC6N_mito_only.vcf \
-V ${dir_path}/SC6T/SC6T_mito_only.vcf \
-V ${dir_path}/SC7N/SC7N_mito_only.vcf \
-V ${dir_path}/SC7T/SC7T_mito_only.vcf \
-V ${dir_path}/SML13/SML13_mito_only.vcf \
-V ${dir_path}/SML1N/SML1N_mito_only.vcf \
-V ${dir_path}/SML1T/SML1T_mito_only.vcf \
-V ${dir_path}/SML4N/SML4N_mito_only.vcf \
-V ${dir_path}/SML4T/SML4T_mito_only.vcf \
-V ${dir_path}/SP2/SP2_mito_only.vcf \
-V ${dir_path}/SP3/SP3_mito_only.vcf \
--genomicsdb-workspace-path ${workdir}/ \
-L HL4_annotated # name of the mitochondiral fasta header


#######
# Genotype samples with GATK's GenotypeGVCFs
#######

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" GenotypeGVCFs \
   -R ${ref_genome} \
   -V gendb:////${dir_path}/Mito_Variant_Calling_0525/ \
   -O ${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito.vcf.gz

#######
# Select SNPs from samples with GATK's SelectVariants
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
 SelectVariants \
    -V /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito.vcf.gz \
    -select-type SNP \
    -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps.vcf.gz

#######
# Select INDELS from samples with GATK's SelectVariants
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
    SelectVariants \
        -V /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito.vcf.gz \
        -select-type INDEL \
        -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_indel.vcf.gz


#######
# Hard Filter SNPs with GATK's VariantFiltration
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      VariantFiltration \
          -V /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2.5" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps_filtered.vcf.gz --verbosity ERROR

#######
# Hard Filter INDELS with GATK's VariantFiltration
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      VariantFiltration \
      -V /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_indel.vcf.gz \
      -filter "QD < 2.0" --filter-name "QD2" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "FS > 200.0" --filter-name "FS200" \
      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
      -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_indel_filtered.vcf.gz \
      --verbosity ERROR

#######
# Merge filtered SNP and INDEL vcfs back together
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      MergeVcfs -I /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps_filtered.vcf.gz \
      -I /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_indel_filtered.vcf.gz \
      -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps-and-indels_filtered.vcf.gz

#######
# Extract the PASS variants from the merged SNP+INDEL vcf
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      SelectVariants -R ${ref_genome} -V /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps-and-indels_filtered.vcf.gz \
      -O /${workdir}/comb_and_genotype_VFC_long_dna/long_TN_mito_snps-and-indels_clean.vcf.gz  --exclude-filtered

