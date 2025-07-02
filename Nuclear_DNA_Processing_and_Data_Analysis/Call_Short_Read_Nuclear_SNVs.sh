###########################
# Call short reads Nuclear SNVs
###########################

########
# set paths and variables
########

sample_name=<name or sampls>
ref_genome=/path_to/HL4.scaffolds.fa
sample_path=/path_to/cleaned_short_reads/
path=/path_to_general_output/
workdir=${path}/short_read_genome_mapping/${sample_name}/
map_name=${path}/short_read_genome_mapping/${sample_name}/${sample_name}_mapped
tempdir=${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_temp

mkdir -p ${tempdir}
mkdir -p ${workdir}
cd ${workdir}

########
# map to HL4 genome
########

# mapping
minimap2 -t 1 -ax sr ${ref_genome} ${sample_path}/${sample_name}_clean_R1_001.fastq.gz ${sample_path}/${sample_name}_clean_R2_001.fastq.gz > ${map_name}.sam

# convert to bam
samtools view -S -b ${map_name}.sam -o ${map_name}.bam

# sort bam file
samtools sort ${map_name}.bam -o ${map_name}.sorted.bam

# modify identifiers
samtools addreplacerg -r 'ID:UVM' -r 'PL:Singularity' -r 'SM:'${sample_name} -o ${map_name}.sorted.un.bam ${map_name}.sorted.bam


# delete superfluous files
rm ${map_name}.sam
rm ${map_name}.bam
rm ${map_name}.sorted.bam

# sort bam 
samtools sort ${map_name}.sorted.un.bam -o ${map_name}.sorted.bam


# mark duplicates with gatk
gatk --java-options "-Xms256G -Xmx256G -XX:ParallelGCThreads=16" \
    MarkDuplicatesSpark \
          -I ${map_name}.sorted.bam \
          -O ${map_name}.sorted.clean.bam --remove-sequencing-duplicates
            
# get mapping stats
samtools flagstat ${map_name}.sorted.clean.bam > ${map_name}.stats_after_mito_mapping.txt


# get mapping depth
mosdepth -t 16 ${sample_name} ${map_name}.sorted.clean.bam

  
# select only mapped bams
samtools view -F 0x04  -b ${map_name}.sorted.clean.bam > ${map_name}.genome_only.bam
samtools index ${map_name}.genome_only.bam


# select unmapped ams for virus / microbe hunting
samtools view -f 4  -b ${map_name}.sorted_for_virus_hunting.bam > ${map_name}_not_genome.bam
samtools index ${map_name}_not_genome.bam

# delete superfluous files
rm ${map_name}.sorted.un.bam
rm ${map_name}.sorted.bam
rm ${map_name}.sorted.clean.ba*


# Haplotype
gatk --java-options "-Xms200G -Xmx200G -XX:ParallelGCThreads=16" \
    HaplotypeCaller \
          -R ${ref_genome} \
          -I ${map_name}.genome_only.bam  \
          -O ./${sample_name}_gatk_2N.vcf \
          -ERC GVCF -ploidy 2


# make a genomics database for genotyping
# Run for each chunk of the genome!!!!!
########    - I ran the first 100 Scaffolds independently and merged them at the end!!!!!!!!!!!!!!!!!!

contig=<genome scaffold> #e.g. HL4_0001

gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" \
GenomicsDBImport \
  -R ${ref_genome} \
  -V ${workdir}BS1/BS1_gatk_2N.vcf \
  -V ${workdir}BS4/BS4_gatk_2N.vcf \
  -V ${workdir}FB11/FB11_gatk_2N.vcf \
  -V ${workdir}FB12/FB12_gatk_2N.vcf \
  -V ${workdir}FB13/FB13_gatk_2N.vcf \
  -V ${workdir}FB15/FB15_gatk_2N.vcf \
  -V ${workdir}FB16/FB16_gatk_2N.vcf \
  -V ${workdir}FB18/FB18_gatk_2N.vcf \
  -V ${workdir}FB1/FB1_gatk_2N.vcf \
  -V ${workdir}FB21/FB21_gatk_2N.vcf \
  -V ${workdir}FB2/FB2_gatk_2N.vcf \
  .......
  .......
  -V ${workdir}SML5/SML5_gatk_2N.vcf \
  -V ${workdir}SML7/SML7_gatk_2N.vcf \
  -V ${workdir}SP2/SP2_gatk_2N.vcf \
  -V ${workdir}SP3/SP3_gatk_2N.vcf \
  --genomicsdb-workspace-path ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig} \
  --tmp-dir ${tempdir} \
  -L ${contig}

# joint genotype all samples by contig
gatk --java-options "-Xmx50g -XX:ParallelGCThreads=2" \
  GenotypeGVCFs \
   -R ${ref_genome} \
   -V gendb:///${workdir}/comb_and_geneotype_VFC_Feb2025/${contig} \
   -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}.vcf.gz

# select SNPs
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
  SelectVariants \
    -V ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}.vcf.gz \
    -select-type SNP \
    -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps.vcf.gz


# select INDELSs
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
  SelectVariants \
        -V ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}.vcf.gz  \
        -select-type INDEL \
        -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_indel.vcf.gz \

# Hard filter SNPs
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
   VariantFiltration \
        -V {workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2.5" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps_filtered.vcf.gz --verbosity ERROR


# Hard filter INDELS
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
  VariantFiltration \
      -V ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_indel.vcf.gz \
      -filter "QD < 2.0" --filter-name "QD2" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "FS > 200.0" --filter-name "FS200" \
      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
      -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_indel_filtered.vcf.gz \
      --verbosity ERROR


# Merge filtered SNP and INDEL vcfs back together
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
  MergeVcfs \
      -I ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps_filtered.vcf.gz \
      -I ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_indel_filtered.vcf.gz \
      -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps_indels_filtered.vcf.gz


# Extract PASS variants only
gatk --java-options "-Xms48G -Xmx48G -XX:ParallelGCThreads=5" \
    SelectVariants \
      -R ${ref_genome} -V ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps_indels_filtered.vcf.gz \
      -O ${workdir}/comb_and_geneotype_VFC_Feb2025/${contig}_snps_indels_clean.vcf.gz  --exclude-filtered
      



