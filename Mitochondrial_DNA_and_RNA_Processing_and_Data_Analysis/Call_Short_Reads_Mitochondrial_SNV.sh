####################
# Call SNVs for short reads DNA and RNA
####################

#####
# assign variables and make directories
#####

sample=<sample name DNA or RNA>
sample_path=<path to fastp output directory>
mito_ref=<path to>/HL4_ann.final-mtgrasp_v1.1.8-assembly.fa
dir_path=<path to short read mitochondrial mapping output>
raw_sample_path=${merged_bam}/bam/All_${sample}.bam
workdir=${dir_path}/Mito_Variant_Calling_short/
fin=${sample_path}/${sample_name}_R1_001.fastq.gz
rin=${sample_path}/${sample_name}_R2_001.fastq.gz
map_name=${workdir}/${sample_name}_mapped


mkdir -p ${workdir}/comb_and_genotype_VFC_short_rna_and_dna/
cd ${workdir}



########
# map to HL4 mitochondria using minimap2
########

minimap2 -t 1 -ax sr ${ref_genome} ${fin} ${rin} > ${map_name}.sam

########
# Use samtools to convert Reads to bam, then sort, update identifiers, subset by mapped reads, and index
########

samtools view -S -b ${map_name}.sam -o ${map_name}.bam
samtools sort ${map_name}.bam -o ${map_name}.sorted.bam
samtools addreplacerg -r 'ID:UVM' -r 'PL:Illumina' -r 'SM:'${sample_name} -o ${map_name}.sorted.un.bam ${map_name}.sorted.bam
samtools view -F 0x04  -b ${map_name}.sorted.un.bam > ${map_name}.mito_only.bam
samtools index ${map_name}.mito_only.bam


########
# remove superfluous  reads
########

rm ${map_name}.sorted.bam
rm ${map_name}.sorted.un.bam
rm ${map_name}.bam

#######
# Call Variants with GATK's Haplotype Caller  
#######

gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=16" HaplotypeCaller -R ${ref_genome} -I ${map_name}.mito_only.bam  -O ./${sample_name}_gatk_N2.vcf -ERC GVCF -ploidy 2

#######
# Make a database of haplotyped calls GATK's GenomicsDBImport
#######

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" GenomicsDBImport \
-R ${ref_genome} \
-V ${workdir}/BS1/BS1_gatk.vcf \
-V ${workdir}/BS4/BS4_gatk.vcf \
-V ${workdir}/FB11/FB11_gatk.vcf \
-V ${workdir}/FB12/FB12_gatk.vcf \
-V ${workdir}/FB13/FB13_gatk.vcf \
-V ${workdir}/FB15/FB15_gatk.vcf \
-V ${workdir}/FB16/FB16_gatk.vcf \
-V ${workdir}/FB18/FB18_gatk.vcf \
-V ${workdir}/FB1/FB1_gatk.vcf \
-V ${workdir}/FB21/FB21_gatk.vcf \
-V ${workdir}/FB2/FB2_gatk.vcf \
-V ${workdir}/FB3/FB3_gatk.vcf \
-V ${workdir}/FB5/FB5_gatk.vcf \
-V ${workdir}/FB8N/FB8N_gatk.vcf \
-V ${workdir}/FB8T/FB8T_gatk.vcf \
-V ${workdir}/HC10/HC10_gatk.vcf \
-V ${workdir}/HC11/HC11_gatk.vcf \
-V ${workdir}/HC12/HC12_gatk.vcf \
-V ${workdir}/HC14/HC14_gatk.vcf \
-V ${workdir}/HC17/HC17_gatk.vcf \
-V ${workdir}/HC18/HC18_gatk.vcf \
-V ${workdir}/HC19/HC19_gatk.vcf \
-V ${workdir}/HC1/HC1_gatk.vcf \
-V ${workdir}/HC27/HC27_gatk.vcf \
-V ${workdir}/HC2N/HC2N_gatk.vcf \
-V ${workdir}/HC2T/HC2T_gatk.vcf \
-V ${workdir}/HC3/HC3_gatk.vcf \
-V ${workdir}/HC4N/HC4N_gatk.vcf \
-V ${workdir}/HC4T/HC4T_gatk.vcf \
-V ${workdir}/HC7/HC7_gatk.vcf \
-V ${workdir}/HL3/HL3_gatk.vcf \
-V ${workdir}/HL4/HL4_gatk.vcf \
-V ${workdir}/JR10/JR10_gatk.vcf \
-V ${workdir}/JR11/JR11_gatk.vcf \
-V ${workdir}/JR12/JR12_gatk.vcf \
-V ${workdir}/JR13/JR13_gatk.vcf \
-V ${workdir}/JR14/JR14_gatk.vcf \
-V ${workdir}/JR18/JR18_gatk.vcf \
-V ${workdir}/JR19/JR19_gatk.vcf \
-V ${workdir}/JR1/JR1_gatk.vcf \
-V ${workdir}/JR2/JR2_gatk.vcf \
-V ${workdir}/JR3/JR3_gatk.vcf \
-V ${workdir}/JR4/JR4_gatk.vcf \
-V ${workdir}/JR9N/JR9N_gatk.vcf \
-V ${workdir}/JR9T/JR9T_gatk.vcf \
-V ${workdir}/LB1/LB1_gatk.vcf \
-V ${workdir}/LB2/LB2_gatk.vcf \
-V ${workdir}/MA10/MA10_gatk.vcf \
-V ${workdir}/MA11/MA11_gatk.vcf \
-V ${workdir}/MA12/MA12_gatk.vcf \
-V ${workdir}/MA13/MA13_gatk.vcf \
-V ${workdir}/MA14/MA14_gatk.vcf \
-V ${workdir}/MA15/MA15_gatk.vcf \
-V ${workdir}/MA18/MA18_gatk.vcf \
-V ${workdir}/MA1N/MA1N_gatk.vcf \
-V ${workdir}/MA1T/MA1T_gatk.vcf \
-V ${workdir}/MA5/MA5_gatk.vcf \
-V ${workdir}/MA7/MA7_gatk.vcf \
-V ${workdir}/MA8/MA8_gatk.vcf \
-V ${workdir}/MA9/MA9_gatk.vcf \
-V ${workdir}/SB10/SB10_gatk.vcf \
-V ${workdir}/SB12/SB12_gatk.vcf \
-V ${workdir}/SB13/SB13_gatk.vcf \
-V ${workdir}/SB14N/SB14N_gatk.vcf \
-V ${workdir}/SB14T/SB14T_gatk.vcf \
-V ${workdir}/SB15/SB15_gatk.vcf \
-V ${workdir}/SB16/SB16_gatk.vcf \
-V ${workdir}/SB1/SB1_gatk.vcf \
-V ${workdir}/SB2/SB2_gatk.vcf \
-V ${workdir}/SB3/SB3_gatk.vcf \
-V ${workdir}/SB4/SB4_gatk.vcf \
-V ${workdir}/SB5/SB5_gatk.vcf \
-V ${workdir}/SB9/SB9_gatk.vcf \
-V ${workdir}/SC11/SC11_gatk.vcf \
-V ${workdir}/SC12/SC12_gatk.vcf \
-V ${workdir}/SC13/SC13_gatk.vcf \
-V ${workdir}/SC16/SC16_gatk.vcf \
-V ${workdir}/SC17/SC17_gatk.vcf \
-V ${workdir}/SC1/SC1_gatk.vcf \
-V ${workdir}/SC3/SC3_gatk.vcf \
-V ${workdir}/SC4/SC4_gatk.vcf \
-V ${workdir}/SC5/SC5_gatk.vcf \
-V ${workdir}/SC6N/SC6N_gatk.vcf \
-V ${workdir}/SC6T/SC6T_gatk.vcf \
-V ${workdir}/SC7N/SC7N_gatk.vcf \
-V ${workdir}/SC7T/SC7T_gatk.vcf \
-V ${workdir}/SML11/SML11_gatk.vcf \
-V ${workdir}/SML12/SML12_gatk.vcf \
-V ${workdir}/SML13/SML13_gatk.vcf \
-V ${workdir}/SML14/SML14_gatk.vcf \
-V ${workdir}/SML15/SML15_gatk.vcf \
-V ${workdir}/SML16/SML16_gatk.vcf \
-V ${workdir}/SML1N/SML1N_gatk.vcf \
-V ${workdir}/SML1T/SML1T_gatk.vcf \
-V ${workdir}/SML2/SML2_gatk.vcf \
-V ${workdir}/SML3/SML3_gatk.vcf \
-V ${workdir}/SML4N/SML4N_gatk.vcf \
-V ${workdir}/SML4T/SML4T_gatk.vcf \
-V ${workdir}/SML5/SML5_gatk.vcf \
-V ${workdir}/SML7/SML7_gatk.vcf \
-V ${workdir}/SP2/SP2_gatk.vcf \
-V ${workdir}/SP3/SP3_gatk.vcf \
-V ${workdir}/vt1_02/vt1_02_gatk.vcf \
-V ${workdir}/vt1_03/vt1_03_gatk.vcf \
-V ${workdir}/vt1_05/vt1_05_gatk.vcf \
-V ${workdir}/vt1_07/vt1_07_gatk.vcf \
-V ${workdir}/vt2_02/vt2_02_gatk.vcf \
-V ${workdir}/vt2_05/vt2_05_gatk.vcf \
-V ${workdir}/vt2_06/vt2_06_gatk.vcf \
-V ${workdir}/vt2_10/vt2_10_gatk.vcf \
-V ${workdir}/vt2_14/vt2_14_gatk.vcf \
-V ${workdir}/vt2_16/vt2_16_gatk.vcf \
-V ${workdir}/vt2_18/vt2_18_gatk.vcf \
-V ${workdir}/vt2_23/vt2_23_gatk.vcf \
-V ${workdir}/vt2_28/vt2_28_gatk.vcf \
-V ${workdir}/vt2_29R/vt2_29R_gatk.vcf \
-V ${workdir}/vt2_30/vt2_30_gatk.vcf \
-V ${workdir}/vt2_32/vt2_32_gatk.vcf \
-V ${workdir}/vt2_37/vt2_37_gatk.vcf \
-V ${workdir}/vt2_40/vt2_40_gatk.vcf \
--genomicsdb-workspace-path ${workdir}/comb_and_genotype_VFC_short_rna_and_dna/ \
-L HL4_annotated # name of the mitochondrial fasta header


#######
# Genotype samples with GATK's GenotypeGVCFs
#######

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" GenotypeGVCFs \
   -R ${ref_genome} \
   -V gendb:///${workdir}/comb_and_genotype_VFC_short_rna_and_dna/\
   -O ${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito.vcf.gz

#######
# Select SNPs from samples with GATK's SelectVariants
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
 SelectVariants \
    -V /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito.vcf.gz \
    -select-type SNP \
    -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps.vcf.gz

#######
# Select INDELS from samples with GATK's SelectVariants
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
    SelectVariants \
        -V /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito.vcf.gz \
        -select-type INDEL \
        -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_indel.vcf.gz

#######
# Hard Filter SNPs with GATK's VariantFiltration
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      VariantFiltration \
          -V /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2.5" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps_filtered.vcf.gz --verbosity ERROR

#######
# Hard Filter INDELS with GATK's VariantFiltration
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      VariantFiltration \
      -V /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_indel.vcf.gz \
      -filter "QD < 2.0" --filter-name "QD2" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "FS > 200.0" --filter-name "FS200" \
      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
      -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_indel_filtered.vcf.gz \
      --verbosity ERROR

#######
# Merge filtered SNP and INDEL vcfs back together
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      MergeVcfs -I /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps_filtered.vcf.gz \
      -I /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_indel_filtered.vcf.gz \
      -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_filtered.vcf.gz

#######
# Extract the PASS variants from the merged SNP+INDEL vcf
#######

gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=2" \
      SelectVariants -R ${ref_genome} -V /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_filtered.vcf.gz \
      -O /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean.vcf.gz  --exclude-filtered


######
# Use vcftools to Filter out SNVs that are not found in at least 90% of samples and with at least 10 x coverage - then change name, compress, and index
######

vcftools --vcf /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean.vcf.gz --minDP 10 --max-missing 0.9  --recode --out /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x

mv /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.recode.vcf /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.vcf

bgzip /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.vcf
tabix /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.vcf.gz

######
# Use BCFtools to extract allele frequency
#####

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.vcf.gz > /${workdir}/comb_and_genotype_VFC_short_rna_and_dna/short_TN_mito_snps-and-indels_clean_0.9_10x.AD.txt


