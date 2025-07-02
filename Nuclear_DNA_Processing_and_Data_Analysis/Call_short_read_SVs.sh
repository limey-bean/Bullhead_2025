
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
#sample=HC4N

ref=/gpfs1/home/e/g/eguswa/scratch/bullhead/Funannotate/HL4/no_mitochondria/HL4_no_mitochondriaFun_out/annotate_results/Ictalurus_punctatus_HL4.scaffolds.fa

input=/gpfs1/home/e/g/eguswa/scratch/bullhead/is_it_transmissible_cancer/nuclear/short_read_genome_mapping/${sample}/${sample}_mapped.genome_only.bam

output=/gpfs1/home/e/g/eguswa/scratch/bullhead/is_it_transmissible_cancer/nuclear/delly/downsampled_130M/

mkdir -p ${output}
cd ${output}

echo ${sample} 
COUNT=$(${samtools} view -@ 10 -c ${input})
echo ${COUNT}
FRACTION=$(awk -v COUNT=${COUNT} "BEGIN {print 130000000 / COUNT}") # use awk to calculate fraction
echo ${FRACTION}


${samtools} view -b -h -@ 5 -s ${FRACTION} ${input} -o ${output}/${sample}_down.bam

${samtools} sort -o ${output}/${sample}_down_sorted.bam -@ 5 ${output}/${sample}_down.bam

${samtools} index -b -@ 10 ${output}/${sample}_down_sorted.bam

${samtools} view -@ 5 -c ${output}/${sample}_down_sorted.bam >> ${output}/final.counts.txt
echo ${sample} >> ${output}/final.counts.txt
rm ${output}/${sample}_down.bam

${delly} call -o ${output}/${sample}.delly.bcf -q 30 -g ${ref} ${output}/${sample}_down_sorted.bam


${gatkcon} bcftools view ${output}/${sample}.delly.bcf  > ${output}/${sample}.delly.vcf

awk ' $0 ~ /^#.*/ && $0 !~ /^##.*/ { split($10, name, "." ); print "CHR", $2, "TYPE", $4, $5, "END", "INS_LEN", ${sample}"_NR", ${sample}"_NV"};
    $0 !~ /^#.*/ { split($8, info, ";" );
        split(info[2], type, "=" );
        split(info[4], end, "=" );
        split(info[5], inslen, "=" );
        split($10, counts, ":" );
        total = counts[11] + counts[12];
        if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] ~ "INS"){print $1, $2, type[2], $4, $5, end[2], inslen[2], total, counts[12]};
        if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] !~ "INS"){ print $1, $2, type[2], $4, $5, end[2], 0, total, counts[12]}
    } ' OFS='\t' ${output}/${sample}.delly.vcf > ${output}/${sample}.delly.table.txt
