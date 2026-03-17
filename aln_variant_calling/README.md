# Aln variant calling

### 01.fastp_clean.sh

Perform default trimming of reads using FASTP

### 02.align_sort_reads.sh

Align reads, sort and filter bam files

### 03.var_calling_freebayes_auto.sh

Run variant calling. First, the command was run without the -@ parameter on 55 samples, and we subsequently called genotypes for only these variants.

### 04.var_calling_freebayes_auto.sh

As above, for the X-chromosome

### 05.vcf_filter_makebed.sh

Filter the vcf file and make bed files used in other analyses.

### 06.prune_bedfiles.sh

Prune the bed files for linkage disequilibrium.