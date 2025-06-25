path=~/concatonated_short_reads/
sample_name=HL4

## Assemble mitogenome
get_organelle_from_reads.py -1 ${path}/${sample_name}_all_R1_001.fastq.gz \
			    -2 ${path}/${sample_name}_all_R2_001.fastq.gz \
           -R 10 -k 21,45,65,85,105 -F animal_mt -o ~/getorganelle/${sample_name}_mt_out --target-genome-size 17000  

## Standardize origin and annotate
mtgrasp_standardize.py -i ~/getorganelle/HL4_mt_out/HL4_complete_mitogenome.fasta -o ${sample_name}_stand  -c 2 -p ${sample_name}_ann -a -mp ~/miniforge3/bin
