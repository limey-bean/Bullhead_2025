# this is the run command

nextflow run nf-core/mag -r 4.0.0 -profile singularity -resume -params-file nf-params.json



# this is the nf-params.json file:

{
   "input": "\/path\/to\/output\/viral_mapping\/DNA_viruses\/samplesheet.csv",
    "outdir": "\/path\/to\/output\/\/viral_mapping\/DNA_viruses\/mag_results",
    "centrifuge_db": "\/path\/to\/Databases\/Centrifuge\/p+h+v.tar.gz",
    "gtdb_db": "\/path\/to\/Databases\/gtdbtk\/gtdbtk_r202_data.tar.gz",
    "kraken2_db": "\/path\/to\/Databases\/Kraken\/minikraken2_v2_8GB_201904.tar.gz",
    "run_virus_identification": true,
    "genomad_db": "\/path\/to\/Databases\/genomad\/genomad_db"
}

###############
# the format of the sample file is below, samplesheet.csv
###############

# you need to locate the unmapped reads from the nuclear genome run mapping and convert then from bam to fastq


# I had to repaire the fastq files

nextflow run nf-core/fastqrepair    -profile singularity    --input samplesheet.csv    --outdir ./ -resume  -r dev


sample	group	short_reads_1	short_reads_2	long_reads	
BS1	0	/path_to/DNA_unmapped_repaired/repaired/BS1_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/BS1_2.fastq.gz		
BS4	0	/path_to/DNA_unmapped_repaired/repaired/BS4_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/BS4_2.fastq.gz		
FB11	0	/path_to/DNA_unmapped_repaired/repaired/FB11_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/FB11_2.fastq.gz		
FB12	0	/path_to/DNA_unmapped_repaired/repaired/FB12_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/FB12_2.fastq.gz		
FB13	0	/path_to/DNA_unmapped_repaired/repaired/FB13_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/FB13_2.fastq.gz		
FB15	0	/path_to/DNA_unmapped_repaired/repaired/FB15_1.fastq.gz	/path_to/DNA_unmapped_repaired/repaired/FB15_2.fastq.gz		
