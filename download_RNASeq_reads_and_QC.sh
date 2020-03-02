#Download rawreads from SRA and perform quality control
#!/bin/bash
#================================================================
#% RUN 
#+    ./download_reference_alignment_RNAseq_reads.sh
#%
#% DESCRIPTION
#%    This is a script will download reads from SRA.
#%		And perform quality assessment & quality control.
#%		Choose "QC" and "ADAPTER" variables based of reads information
#%		QC(Default:yes) and ADAPTER(Default:NEXTERA_UNIVERSAL) 
#%    In order run the analysis in a specified folder and path,
#%		then set the name in "project_config.sh" file.
#%    (Default:RNAseq_project)and (Default:current working directory).
#%
#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================

# Prepare the reads for analysis
# Select adapter sequences
ILLUMINA_TRUSEQ_LT_HT_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ILLUMINA_TRUSEQ_LT_HT_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
NEB_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
NEB_R2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"
NEXTERA_UNIVERSAL="CTGTCTCTTATACACATCT"
TRUSEQ_SRNA="TGGAATTCTCGGGTGCCAAGG"
RIBO_PRO_NEB_R1="CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
SURESELECTQXT="CTGTCTCTTGATCACA"

# Are reads Quality controlled?
QC=yes #yes or no 
# Choose the adapter from above list!
ADAPTER=$NEXTERA_UNIVERSAL

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101753
SRA_IDS="SRR5858228
         SRR5858229
         SRR5858230
         SRR5858231
         SRR5858232
         SRR5858233
         SRR5858234
         SRR5858235
         SRR5858236"


main(){
  import_project_config
  set_path
  create_folders
  download_SRA_data
  quality_control
}

import_project_config(){
. ${PWD}/project_config.sh
}

set_path(){
  FOLDER_NAME=${PROJECT_NAME}
  FOLDER_PATH=${PROJECT_PATH}
  INPUT_FOLDER=${FOLDER_PATH}/${FOLDER_NAME}
  SRA_BIN=${INPUT_FOLDER}/tools/sratoolkit/bin/
  FASTQC_BIN=${INPUT_FOLDER}/tools/FastQC/fastqc
  CUTADAPT_BIN=${INPUT_FOLDER}/tools/cutadapt/bin/cutadapt
}


create_folders(){
  echo "Creating Folder Structure"
  mkdir -p ${INPUT_FOLDER}
  mkdir -p ${INPUT_FOLDER}/data
  mkdir -p ${INPUT_FOLDER}/analysis
  mkdir -p ${INPUT_FOLDER}/tools
  mkdir -p ${INPUT_FOLDER}/reads
}


download_SRA_data(){
  # Download SRA data
  echo "Downloading SRA data"
  mkdir -p ${INPUT_FOLDER}/reads/rawreads
  cd ${INPUT_FOLDER}/reads/rawreads
  for SRA_ID in ${SRA_IDS}
   do
	  echo ""
     echo "#### Running for SRA ID $SRA_ID ####"
	  echo ""
     ${SRA_BIN}/fastq-dump \
     ${SRA_ID} \
     --split-files \
     --origfmt \
	  --gzip
   done
}


run_preQC_stats(){
  echo "#### Running pre QC stats using fastqc ####"
  mkdir -p ${INPUT_FOLDER}/reads/rawreads_QA_stats

  for FILE in $(ls ${INPUT_FOLDER}/reads/rawreads/*.fastq.gz)
    do
		echo ""
      echo "Running preQC on FILE $FILE"
		echo ""
      ${FASTQC_BIN} \
      ${FILE} \
      -o ${INPUT_FOLDER}/reads/rawreads_QA_stats
    done
}


run_cutadapt(){
  echo "#### Running cutadapt ####"
  mkdir -p ${INPUT_FOLDER}/reads/filtered_reads

  for FILE in $(ls ${INPUT_FOLDER}/reads/rawreads/*.fastq.gz)
    do
      FILENAME=$(basename ${FILE} | sed -e "s/.fastq.gz//")
      echo ""
		echo "Running on FILE $FILENAME"
      echo ""
		${CUTADAPT_BIN} \
      -q 20 \
      --minimum-length 35 \
      -a ${ADAPTER} \
      -o ${INPUT_FOLDER}/reads/filtered_reads/${FILENAME}_trimmed.gz \
      ${FILE} \
      > ${INPUT_FOLDER}/reads/filtered_reads/${FILENAME}_cutadapt_stats.txt
    done
}

main

