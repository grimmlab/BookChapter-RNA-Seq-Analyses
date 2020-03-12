#!/bin/bash
#================================================================
#% RUN 
#+    ./quantification.sh
#%
#% DESCRIPTION
#%    This is a script will perform quantification. Mutilple options
#%    are provided to perform quantification. It is possible to 
#%    perform different quantifications (gene, exon, isomforms)
#%    using this script. Default is gene quantification from 
#%    bedtools. In order to change theher quantification methods then 
#%    set it in "project_config.sh" file.
#%    In order run the analysis in a specified folder and path,
#%    then set the name in "project_config.sh" file.
#%    (Default:RNAseq_project) and (Default:current working directory).
#%
#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================

main(){
  import_project_config
  set_path
  create_folders
  run_counting_reads_per_gene_quantification
  run_counting_reads_per_transcript_quantification  
  run_counting_reads_per_exon_quantification
}

import_project_config(){
. ${PWD}/project_config.sh
	GENOME=${GENOME_TYPE}
	SOURCE=gencode
	VERSION=${GENOME_VERSION}
	GENOME_ANNO_FOLDER=${GENOME}-${SOURCE}-version-${VERSION}
	MAPPER=${ALIGNER}
	QUANT_COUNT_METHOD=${QUANT_COUNT}
	THREADS=8
}



set_path(){
  FOLDER_NAME=${PROJECT_NAME}
  FOLDER_PATH=${PROJECT_PATH}
  INPUT_FOLDER=${FOLDER_PATH}/${FOLDER_NAME}
  GENOME_FASTA=${FOLDER_PATH}/${FOLDER_NAME}/data/${GENOME_ANNO_FOLDER}/GRCm38.p6.genome.fa
  ANNOTATION_GTF=${FOLDER_PATH}/${FOLDER_NAME}/data/${GENOME_ANNO_FOLDER}/gencode.vM24.annotation.gtf
  ANNOTATION_GFF=${FOLDER_PATH}/${FOLDER_NAME}/data/${GENOME_ANNO_FOLDER}/gencode.vM24.annotation.gff3
  READS_FOLDER=${INPUT_FOLDER}/reads/filtered_reads
  MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/${MAPPER}/
  HTSEQ_BIN=${INPUT_FOLDER}/tools/htseq/bin/
  BEDTOOLS_BIN=${INPUT_FOLDER}/tools/bedtools2/bin/bedtools
  CUFFLINKS_BIN=${INPUT_FOLDER}/tools/cufflinks/
  SAMTOOLS_BIN=${INPUT_FOLDER}/tools/samtools/samtools
  DEXSEQ_BIN=${INPUT_FOLDER}/tools/DEXSeq/

}


create_folders(){
  mkdir -p ${INPUT_FOLDER}
  mkdir -p ${INPUT_FOLDER}/data
  mkdir -p ${INPUT_FOLDER}/analysis
  mkdir -p ${INPUT_FOLDER}/tools
  mkdir -p ${INPUT_FOLDER}/reads
}


run_counting_reads_per_gene_quantification(){
  if [ "${QUANT_COUNT_METHOD}" == "per_gene_bedtools" ]
  then
    echo "#### Running reads per gene quantification using BEDTOOLS ####"
    run_gene_quantification_with_bedtools
  elif [ "${QUANT_COUNT_METHOD}" == "per_gene_htseq" ]
  then
    echo "#### Running reads per gene quantification using HTSEQ ####"
    run_gene_quantification_with_htseq
  else
    echo "#### ERROR! Please choose BEDTOOLS or HTSEQ as quantification tool ####"
  fi
}


run_counting_reads_per_transcript_quantification(){
  if [ "${QUANT_COUNT_METHOD}" == "per_transcript_cufflinks" ]
    then
      echo ""
      echo "#### Running reads per transcript quantification using CUFFLINKS ####"
      echo ""
      run_transcript_quantification_with_cufflinks
   else
      echo ""
      echo "#### Skip running reads per transcript quantification using CUFFLINKS ####"
      echo ""
   fi 
}


run_counting_reads_per_exon_quantification(){
  if [ "${QUANT_COUNT_METHOD}" == "per_exon_dexseq" ]
    then
		echo ""
      echo "#### Running reads per exon quantification using DEXSEQ ####"
      echo ""
      run_exon_quantification_with_dexseq
  else
      echo ""
      echo "#### Skip running reads per exon quantification using DEXSEQ ####"
      echo ""
   fi
}

run_transcript_quantification_with_cufflinks(){
  CUFFLINKS_QUANT_FOLDER=${INPUT_FOLDER}/analysis/quantification/cufflinks-count  
  if [ ! -d "${CUFFLINKS_QUANT_FOLDER}" ]; then
    for BAM in $(ls ${MAPPING_FOLDER}/*bam)
      do
        mkdir -p ${CUFFLINKS_QUANT_FOLDER}
        echo ""
        echo "Running CUFFLINKS on $BAM file"
        echo ""
        OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam//")
        ${CUFFLINKS_BIN}/cufflinks \
        -G ${ANNOTATION_GTF} \
        -b ${GENOME_FASTA} \
        -p ${THREADS} \
        $BAM \
        -o ${CUFFLINKS_QUANT_FOLDER}

        mv ${CUFFLINKS_QUANT_FOLDER}/isoforms.fpkm_tracking ${CUFFLINKS_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_isoforms.fpkm_tracking)
        mv ${CUFFLINKS_QUANT_FOLDER}/genes.fpkm_tracking ${CUFFLINKS_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_genes.fpkm_tracking)
        mv ${CUFFLINKS_QUANT_FOLDER}/transcripts.gtf ${CUFFLINKS_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_transcripts.gtf)
        mv ${CUFFLINKS_QUANT_FOLDER}/skipped.gtf ${CUFFLINKS_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_skipped.gtf)
        echo "DONE!"
    done
  else
    echo ""
    echo "CUFFLINKS folder exists! Skipping CUFFLINKS quantification"
    echo ""
  fi
}


run_exon_quantification_with_dexseq(){
  DEXSEQ_QUANT_FOLDER=${INPUT_FOLDER}/analysis/quantification/dexseq-count
  ANNOTATION_GTF_DEX=${DEXSEQ_QUANT_FOLDER}/DEXSEQ_GTF_annotation.gff
  if [ ! -d "${DEXSEQ_QUANT_FOLDER}" ]; then
    mkdir -p ${DEXSEQ_QUANT_FOLDER}
    echo ""
    echo "Preparing annotation for DEXSEQ"
	 echo ""
    python2 ${DEXSEQ_BIN}/inst/python_scripts/dexseq_prepare_annotation.py \
	 ${ANNOTATION_GTF} \
    ${ANNOTATION_GTF_DEX}
    
    for BAM in $(ls $MAPPING_FOLDER/*bam)
      do
        mkdir -p ${DEXSEQ_QUANT_FOLDER}
        echo ""
        echo "Running DEXSEQ on $BAM file"
		  echo ""
        OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam//")
        #${SAMTOOLS_BIN} view \
        #-h \
        #-o ${OUTPUT_FILE_PREFIX}.sam \
        #${BAM}

        python2 ${DEXSEQ_BIN}/inst/python_scripts/dexseq_count.py \
        -p no \
        -s no \
        -r name \
        ${ANNOTATION_GTF_DEX} \
        -f bam \
	     -a 10 \
        ${BAM} \
        $DEXSEQ_QUANT_FOLDER/$(basename ${OUTPUT_FILE_PREFIX}_exon_counts.csv)
      done
  else
    echo ""
    echo "DEXSEQ folder exists! Skipping DEXSEQ quantification"
    echo ""
  fi
}


run_gene_quantification_with_htseq(){
  echo ""
  echo "Running per gene using HTSEQ"
  echo ""
  HTSEQ_QUANT_FOLDER=${INPUT_FOLDER}/analysis/quantification/htseq-count
  if [ ! -d "${HTSEQ_QUANT_FOLDER}" ]; then
    for BAM in $(ls ${MAPPING_FOLDER}/*bam)
      do
        mkdir -p ${HTSEQ_QUANT_FOLDER}
        OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam//")
		  echo $OUTPUT_FILE_PREFIX
        ${HTSEQ_BIN}/htseq-count \
        -f bam \
        -a 10 \
        -m intersect-strict \
        -s no \
        -t exon \
        -i gene_id \
		  ${BAM} \
        ${ANNOTATION_GTF} > ${HTSEQ_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_counts.csv)
     done
  else
    echo ""
    echo "HTSEQ folder exists! Skipping HTSEQ quantification"
    echo ""
  fi

  #merging the files
  rm ${HTSEQ_QUANT_FOLDER}/Read_per_features_combined.csv
  echo "Attributes" > tmp_combined
  tmp_file=$(ls -t ${HTSEQ_QUANT_FOLDER}/*counts.csv |head -n 1)
  cut -f -1 ${tmp_file} \
  >> tmp_combined
  
  for FILE in $(ls ${HTSEQ_QUANT_FOLDER})
    do
      CLEANED_NAME=$(echo $FILE | sed "s/_trimmed_counts.csv//")
      echo ${HTSEQ_QUANT_FOLDER}/$FILE
      echo $CLEANED_NAME > tmp1
      cut -f 2 ${HTSEQ_QUANT_FOLDER}/$FILE >> tmp1
      cp tmp_combined tmp_combined_curr
      paste tmp_combined_curr tmp1 > tmp_combined
      rm tmp1 tmp_combined_curr
    done
  cat tmp_combined | grep -vE  "__" > ${HTSEQ_QUANT_FOLDER}/Read_per_features_combined.csv
  rm tmp_combined
  echo "DONE!"
}



run_gene_quantification_with_bedtools(){
  echo ""
  echo "Running per gene using BEDTOOLS"
  echo ""
  mkdir -p ${INPUT_FOLDER}/analysis/quantification/tmp
  Q_TMP=${INPUT_FOLDER}/analysis/quantification/tmp
  BEDTOOLS_QUANT_FOLDER=${INPUT_FOLDER}/analysis/quantification/bedtools-count
  grep -P "\tgene\t.*protein_coding|^#" ${ANNOTATION_GFF} \
  > ${Q_TMP}/tmp_genes_only.gff3
  if [ ! -d "${BEDTOOLS_QUANT_FOLDER}" ]; then
    for BAM in $(ls ${MAPPING_FOLDER}/*bam)
      do
        mkdir -p ${BEDTOOLS_QUANT_FOLDER}
        OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam//")
        ${BEDTOOLS_BIN} intersect \
        -S \
        -wa \
        -c \
        -a ${Q_TMP}/tmp_genes_only.gff3 \
        -b ${BAM} \
        > ${BEDTOOLS_QUANT_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_counts.csv)
      done
      rm -r ${Q_TMP}

      echo -e "Chr\tSource\tFeature\tStart\tEnd\tFrame1\tStrand\tFrame2\tAttributes" \
       > tmp_combined
        tmp_file=$(ls -t ${BEDTOOLS_QUANT_FOLDER}/*counts.csv | head -n 1)
       cut -f -9 ${tmp_file} \
      >> tmp_combined

    rm ${BEDTOOLS_QUANT_FOLDER}/Read_per_features_combined.csv
    for FILE in $(ls ${BEDTOOLS_QUANT_FOLDER})
    do
      CLEANED_NAME=$(echo $FILE | sed "s/_trimmed.Aligned.sortedByCoord.out.csv//")
      echo ${BEDTOOLS_QUANT_FOLDER}/$FILE
      echo $CLEANED_NAME > tmp1
      cut -f 10 ${BEDTOOLS_QUANT_FOLDER}/$FILE >> tmp1
      cp tmp_combined tmp_combined_curr
      paste tmp_combined_curr tmp1 > tmp_combined
      rm tmp1 tmp_combined_curr
    done
    mv tmp_combined ${BEDTOOLS_QUANT_FOLDER}/Read_per_features_combined.csv
    echo "DONE!"
  else
    echo ""
    echo "BEDTOOLS folder exists! Skipping BEDTOOLS quantification"
    echo ""
  fi
}

main
