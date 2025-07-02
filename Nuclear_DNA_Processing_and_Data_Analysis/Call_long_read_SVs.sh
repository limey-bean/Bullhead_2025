


sample=<sample name>

out=/path_to_output_for/sniffles/

datadir=/path_to/map_long/
ref_genome=/path_to/HL4.scaffolds.fa

mkdir ${out}

sniffles --input ${datadir}/${sample}/${sample}_aligned_sorted_filtered.bam --snf ${out}/${sample}.snf

sniffles --input ${out}/${sample1}.snf ${out}/${sample2}.snf ... ${out}/${sampleN}.snf --vcf sniffles_all.vcf
