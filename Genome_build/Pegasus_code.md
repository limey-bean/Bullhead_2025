module load singularity

~/PEGASUS/pegasus.sh \
            -n ~/HL4/Long_reads/HL4_no_mito.fastq.gz \
            -s1 ~/HL4/reads/short_reads1  \
            -s2 ~/HL4/reads/short_reads2 \
            -phv ~/PHVindexes\
            -b ~/actinopterygii_odb10\
            -qg ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.gff \
            -qf ~/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.fna \
            -r ~/HL4RefGen.fasta
