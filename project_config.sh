PROJECT_NAME=RNASeq_project
PROJECT_PATH=${PWD}

# Parameters for download_genome_and_annotation.sh
GENOME_TYPE=mouse # or human
GENOME_VERSION=24 #or 33 for human


# Parameters for reference_based_alignment.sh
ALIGNER=star #tophat or star

# Parameters for quantification.sh
QUANT_COUNT=per_gene_bedtools #per_gene_bedtools #per_gene_htseq or per_transcript_cufflinks or per_exon_dexseq
