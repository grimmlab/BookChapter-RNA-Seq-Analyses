#!/bin/bash
#================================================================
#% RUN 
#+    ./reference_based_alignment.sh
#%
#% DESCRIPTION
#%    This is a script will perform reference based alignment.
#%    Alignment can be performed for Human or Mouse genomes.
#%		The script uses two alignment tools: STAR (Default) and TOPHAT
#%		and tool name can be changed in "project_config.sh".
#%    In order run the analysis in a specified folder and path,
#%    then set the name in "project_config.sh" file.
#%    (Default:RNAseq_project)and"(Default:current working directory).
#%
#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================
#Running reference based alignment using STAR and TopHat2

main(){
  import_project_config
  set_path
  create_folders
  run_alignment
}

import_project_config(){
. ${PWD}/project_config.sh
	ALIGNMENT=${ALIGNER}
	GENOME=${GENOME_TYPE}
	SOURCE=gencode
	VERSION=${GENOME_VERSION}
	GENOME_ANNO_FOLDER=${GENOME}-${SOURCE}-version-${VERSION}
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
  STAR_BIN=${INPUT_FOLDER}/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR
  SAMTOOLS_BIN=${INPUT_FOLDER}/tools/samtools/samtools
  TOPHAT_BIN=${INPUT_FOLDER}/tools/tophat/tophat #Note both Bowtie2 and SAMtools have to be on the path
  BOWTIE_BIN=${INPUT_FOLDER}/tools/bowtie2/
  RSEQC_BIN=${INPUT_FOLDER}/tools/rseqqc/bin
}

create_folders(){
  mkdir -p ${INPUT_FOLDER}
  mkdir -p ${INPUT_FOLDER}/data
  mkdir -p ${INPUT_FOLDER}/analysis
  mkdir -p ${INPUT_FOLDER}/tools
  mkdir -p ${INPUT_FOLDER}/reads
  }

run_alignment(){
  if [ "${ALIGNMENT}" == "star" ]
  then
	 echo "#### Running STAR alignment steps ####"
    run_genome_indexing_with_STAR
    run_alignment_with_STAR
	 run_index_STAR_bam_files
    run_mapping_stats_STAR
  elif [ "${ALIGNMENT}" == "tophat" ]
  then
	 echo "#### Running TOPHAT alignment steps ####"
    run_genome_indexing_with_BOWTIE2
	 run_alignment_with_TOPHAT
	 run_index_TOPHAT_bam_files
    run_mapping_stats_TOPHAT
  else
	 echo "#### ERROR! Please choose STAR or TOPHAT as alignment tool ####"
  fi	 
}


#################################################################################################################
                                                  #STAR#
#################################################################################################################

  run_genome_indexing_with_STAR(){
  #Running STAR Indexing
  echo ""
  echo "Running STAR genome indexing"
  echo ""
  #mkdir -p ${INPUT_FOLDER}/data/STAR_Genome_Index
  GENOME_INDEX_DIR_STAR=${INPUT_FOLDER}/data/STAR_Genome_Index
  if [ ! -d "${GENOME_INDEX_DIR_STAR}" ]; then
    mkdir -p ${GENOME_INDEX_DIR_STAR}
    $STAR_BIN \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir $GENOME_INDEX_DIR_STAR \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $ANNOTATION_GTF
    echo "DONE!"
  else
    echo "STAR genome index folder exists! Skipping Genome indexing"
  fi
}

  run_alignment_with_STAR(){
  #Running STAR alignment
  echo ""
  echo "Running STAR reference alignment"
  echo ""
  GENOME_INDEX_DIR_STAR=${INPUT_FOLDER}/data/STAR_Genome_Index
  STAR_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/star
  if [ ! -d "${STAR_MAPPING_FOLDER}" ]; then
    for LIB in $(ls $READS_FOLDER/*.gz)
      do
        mkdir -p ${STAR_MAPPING_FOLDER}
        OUTPUT_FILE_PREFIX=$(echo $LIB | sed -e "s/.gz/./")
        echo ""
        echo "Running alignment for $LIB"
        echo ""
        $STAR_BIN \
        --runMode alignReads \
        --runThreadN ${THREADS} \
        --genomeDir $GENOME_INDEX_DIR_STAR \
        --readFilesIn $LIB \
        --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION_GTF \
        --outFileNamePrefix $STAR_MAPPING_FOLDER/$(basename $OUTPUT_FILE_PREFIX) \
        --outSAMtype BAM SortedByCoordinate
    done
	 echo "DONE!"
  else
    echo ""
    echo "STAR mapping folder exists! Skipping STAR alignment"
    echo ""
  fi
}

  run_index_STAR_bam_files(){
  echo ""
  echo "Running samtools indexing for STAR BAM files"
  echo ""
  STAR_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/star/
  if test -n "$(find  ${STAR_MAPPING_FOLDER}/ -maxdepth 1 -name '*.bai' -print -quit)"
    then
		echo "" 
      echo "Indexed bam files from STAR already exists! Skipping STAR bam indexing"
      echo ""
  else
    for BAM in $(ls $STAR_MAPPING_FOLDER | grep bam$)
      do
        echo ""
        echo "Running samtools for indexing STAR generated $BAM file"
		  echo ""
        ${SAMTOOLS_BIN} \
        index $STAR_MAPPING_FOLDER/$BAM
        echo "DONE!"
     done
	  wait
  fi
}

  run_mapping_stats_STAR(){
  echo ""
  echo "Running STAR alignment stats using RSeQC"
  echo ""
  
  STAR_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/star/
  if test -a "$(find  ${STAR_MAPPING_FOLDER}/ -maxdepth 1 -name '*.bam' -print -quit)"
    then  
      for BAM in $(ls $STAR_MAPPING_FOLDER | grep bam$)
        do
          echo ""
          echo "Running samtools for STAR generated $BAM file"
          echo "" 
	       OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam/./")
 	       python3 $RSEQC_BIN/bam_stat.py \
			 -q 30 \
	       -i $STAR_MAPPING_FOLDER/$BAM > ${STAR_MAPPING_FOLDER}/${OUTPUT_FILE_PREFIX}_RSeQC_alignment_stats.txt
	   done
  grep "Uniquely mapped reads %" $STAR_MAPPING_FOLDER/*Log.final.out | sed -e "s/.Log.final.out//" -e s/:.*\|// \
  > $STAR_MAPPING_FOLDER/Percentage_uniquely_mapped_reads.csv
  else
    echo ""
    echo "Missing bam file.Cannot run alignment stats!"
    echo ""
  fi
}


#################################################################################################################
                                                  #TOPHAT#
#################################################################################################################

  #PART-2 Alignment with TOPHAT
  
  run_genome_indexing_with_BOWTIE2(){
  echo ""
  echo "Running genome indexing for TOPHAT"
  echo ""
  GENOME_INDEX_DIR_TOPHAT=${INPUT_FOLDER}/data/TOPHAT_Genome_Index
  GENOME_FILE_PREFIX=$(echo $(basename ${GENOME_FASTA}) | sed -e "s/.fa//")

  if [ ! -d "${GENOME_INDEX_DIR_TOPHAT}" ]; then
    mkdir -p ${GENOME_INDEX_DIR_TOPHAT}
    cd ${GENOME_INDEX_DIR_TOPHAT}
    ${BOWTIE_BIN}/bowtie2-build \
    -f $GENOME_FASTA \
    ${GENOME_FILE_PREFIX} \
    -p ${THREADS}
    echo "DONE!"
  else
    echo "TOPHAT genome index folder exists! Skipping Genome indexing"
  fi
}

  run_alignment_with_TOPHAT(){
  #Running TopHat2 alignment
  echo ""
  echo "Running TOPHAT reference alignment"
  echo ""
  GENOME_INDEX_DIR_TOPHAT=${INPUT_FOLDER}/data/TOPHAT_Genome_Index
  cp ${GENOME_FASTA} ${GENOME_INDEX_DIR_TOPHAT}
  TOPHAT_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/tophat/
  GENOME_FILE_PREFIX=$(echo $(basename ${GENOME_FASTA}) | sed -e "s/.fa//")

  if [ ! -d "${TOPHAT_MAPPING_FOLDER}" ]; then
    for LIB in $(ls $READS_FOLDER/*.gz)
      do
		  mkdir -p ${INPUT_FOLDER}/analysis/mapping/tophat/
		  OUTPUT_FILE_PREFIX=$(echo $LIB | sed -e "s/.gz//")
		  echo ""
		  echo "Running TOPHAT alignment for $LIB"
		  echo ""
		  $TOPHAT_BIN \
		  -o ${TOPHAT_MAPPING_FOLDER} \
		  -p ${THREADS} \
		  -G $ANNOTATION_GTF \
		  ${GENOME_INDEX_DIR_TOPHAT}/${GENOME_FILE_PREFIX} \
		  $LIB
		  mv ${TOPHAT_MAPPING_FOLDER}/accepted_hits.bam ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_accepted_hits.bam)
        mv ${TOPHAT_MAPPING_FOLDER}/unmapped.bam ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_unmapped.bam)
        mv ${TOPHAT_MAPPING_FOLDER}/align_summary.txt ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_align_summary.txt)
        mv ${TOPHAT_MAPPING_FOLDER}/deletions.bed ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_deletion.bed)
        mv ${TOPHAT_MAPPING_FOLDER}/junctions.bed ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_junctions.bed)
        mv ${TOPHAT_MAPPING_FOLDER}/prep_reads.info ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_prep_reads.info)
        mv ${TOPHAT_MAPPING_FOLDER}/logs ${TOPHAT_MAPPING_FOLDER}/$(basename ${OUTPUT_FILE_PREFIX}_logs)

    done
    mkdir -p ${TOPHAT_MAPPING_FOLDER}/unmapped
    mv ${TOPHAT_MAPPING_FOLDER}/*_unmapped.bam ${TOPHAT_MAPPING_FOLDER}/unmapped
    echo "DONE!"
   else
     echo ""
     echo "TOPHAT mapping folder exists! Skipping TOPHAT alignment"
     echo ""
  fi
}

  run_index_TOPHAT_bam_files(){
  echo ""
  echo "Running samtools indexing for TOPHAT BAM files"
  echo ""
  TOPHAT_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/tophat/
  if stat --printf='' ${TOPHAT_MAPPING_FOLDER}/*.bai 2>/dev/null
    then
		echo ""
      echo "Indexed bam files from TOPHAT already exists! Skipping TOPHAT bam indexing"
      echo ""
  else
    for BAM in $(ls ${TOPHAT_MAPPING_FOLDER} | grep bam$)
      do
        echo ""
        echo "Running samtools for indexing TOPHAT generated $BAM file"
        echo ""
        ${SAMTOOLS_BIN} \
        index ${TOPHAT_MAPPING_FOLDER}/${BAM}
        echo "DONE!"
     done
     wait
  fi
}

  run_mapping_stats_TOPHAT(){
  echo ""
  echo "Running TOPHAT alignment stats using RSeQC"
  echo ""
  TOPHAT_MAPPING_FOLDER=${INPUT_FOLDER}/analysis/mapping/tophat/
  if stat --printf='' ${TOPHAT_MAPPING_FOLDER}/*.bam 2>/dev/null
    then
      for BAM in $(ls ${TOPHAT_MAPPING_FOLDER} | grep bam$)
        do
        echo ""
        echo "Running alignment stats for TOPHAT generated $BAM file"
        echo "" 
        OUTPUT_FILE_PREFIX=$(echo $BAM | sed -e "s/.bam//")
        python3 ${RSEQC_BIN}/bam_stat.py \
        -i ${TOPHAT_MAPPING_FOLDER}/${BAM}> ${TOPHAT_MAPPING_FOLDER}/${OUTPUT_FILE_PREFIX}_RSeQC_alignment_stats.txt
      done
  else
    echo ""
    echo "Missing bam file.Cannot run alignment stats!"
    echo ""
  fi 
}

main
