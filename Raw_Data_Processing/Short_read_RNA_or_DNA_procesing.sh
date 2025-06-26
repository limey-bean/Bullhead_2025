
# All samples with multiple sequence output were concatonated by forward or reverse reads using cat and stored in ~/raw_concatonated_short_reads/

in=~/raw_concatonated_short_reads/
fastqc_out=~/short_reads/fastqc_short/
fastp_out=~/short_reads/fastp/

sample_name=<sample name>

# Run fastqc on all samples
fastqc ${in}* -o ${fastqc_out}

in=~/raw_concatonated_short_reads/
out=~/short_reads/fastp/

fastp -i ${in}/${sample_name}_all_R1_001.fastq.gz -I ${in}/${sample_name}_all_R2_001.fastq.gz -o ${fastp_out}${sample_name}_R1_001.fastq.gz -O ${fast_pout}${sample_name}_R2_001.fastq.gz -h ${fastp_out}${sample_name}_fastp.html -j ${fastp_out}${sample_name}_fastp.json

