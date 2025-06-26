#######
# Set variables
#######

sample_name=<sample_name>

in=<path to raw concatonated short reads>
fastqc_out=<path to fastqc output directory>
fastp_out=<path to fastp output directory>

#######
# make directories
#######

mkdir -p ${in}
mkdir -p ${fastqc_out}
mkdir -p ${fastp_out}

######
# Concatonate raw fastq data (if necessary)
######

# If there are multiple sequence output per sample, concatonate by forward and reverse reads:

## compile a list of all raw reads for a sample, below I am calling the place where I store my raw reads "~/raw_data_dump_dir/":
find ~/raw_data_dump_dir/ -type f -name "${sample_name}_*R1_001*" > ${workdir}/${sample_name}_r_R1_filepaths.txt
find ~/raw_data_dump_dir/  -type f -name "${sample_name}_*R2_001*" > ${workdir}/${sample_name}_r_R2_filepaths.txt

## concatonate reads
{ xargs cat < ${workdir}/${sample_name}_R1_filepaths.txt ; } > ${workdir2}/${sample_name}_all_R1_001.fastq.gz
{ xargs cat < ${workdir}/${sample_name}_R2_filepaths.txt ; } > ${workdir2}/${sample_name}_all_R2_001.fastq.gz

########
# Run fastqc on all samples
########
fastqc ${in}* -o ${fastqc_out}

########
# Run fastp on all samples
########
fastp -i ${in}/${sample_name}_all_R1_001.fastq.gz -I ${in}/${sample_name}_all_R2_001.fastq.gz -o ${fastp_out}${sample_name}_R1_001.fastq.gz -O ${fast_pout}${sample_name}_R2_001.fastq.gz -h ${fastp_out}${sample_name}_fastp.html -j ${fastp_out}${sample_name}_fastp.json

