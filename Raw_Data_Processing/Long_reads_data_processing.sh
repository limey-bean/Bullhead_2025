########################################################
# Basecalling/quality filtering and Demuxing for ONT Longreads was done using dorado.
# Nanopore sequencers collect sequence and base modification signal data and store it in a pod5 file. 
# We re-basecalled all long reads with the same version of dorado, but only called one run at a time. Bam files were merged by sample after basecalling.

# all new directories were made prior to running software.

############
# Basecall with dorado version 0.7.2 - dorado-0.7.2-linux-x64
# used high accuracy and called methylation for a different project at the same time.
# be careful with the trim option… If you trim before you demux, you will remove the bases used to sort barcodes. If you have only a single sample use --trim all.
# filtered out reads with q score < 7
############
pod5_in=<path to the pod5 directory for a single run>
base_called_output=<path to the basecalled bam file>

dorado basecaller hac,5mCG_5hmCG \
  --min-qscore 7 --kit-name SQK-NBD114-24  --no-trim \
  ${pod5_in} > ${base_called_output}

########
# Demux using dorado version 0.7.2
#######

demuxed_samples=<path to the demuxed samples output directory>

dorado \
  demux --kit-name SQK-NBD114-24 --output-dir ${demuxed_samples} ${base_called_output}

#######
# Rename demuxed bam files
#######

renamed_demuxed_bam=<dir path to store the renamed demuxed bam files>

mkdir -p ${renamed_demuxed_bam}
sample=<name of sample>
barcode=<barcode assigned to a sample>

cp ${demuxed_samples}/SQK-NBD114-24_barcode${barcode}.bam ${renamed_demuxed_bam}/${sample}.bam

######
# Merge bam samples if sample was sequenced across multiple runs
#####

merged_bam=<path to directory to store merged bam files>
samtools merge ${merged_bam}/bam/All_${sample}.bam ${renamed_demuxed_bam1}/${sample}.bam ${renamed_demuxed_bam2}/${sample}.bam ....


########
# Run fastqc on all samples
########

fastqc_out=<path to fastqc output directory>
fastqc ${merged_bam}/* -o ${fastqc_out}

