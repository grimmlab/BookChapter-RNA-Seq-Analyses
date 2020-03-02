#Download and install required for the reference based alignment
#!/bin/bash
#================================================================
#% RUN 
#+    ./download_genome_and_annotation.sh
#%
#% DESCRIPTION
#%    This is a script will download human or mouse genome and its
#%    annotation as defined in "project_config.sh" file (Default:mouse)
#%    In order run the analysis in a specified folder and path,
#%    then set it in "project_config.sh file" 
#%		(default:RNASeq_project) and (default:current working directory)
#%
#================================================================
#- IMPLEMENTATION
#-    version         0.0.1
#-    author          Richa Bharti
#-    copyright       Copyright (c) 2020
#-    license         GNU General Public License
#================================================================

#Download annotation and genome
main(){
  import_project_config
  set_path # -> Donot comment this, has to be always active
	create_folders
	download_genome
	calculate_md5sum
	decompress_genome_annotation_files
	make_genome_annotation_unwritetable
}

import_project_config(){
. ${PWD}/project_config.sh
}


set_path(){
  FOLDER_NAME=${PROJECT_NAME}
  FOLDER_PATH=${PROJECT_PATH}
  INPUT_FOLDER=${FOLDER_PATH}/${FOLDER_NAME}
  VERSION=${GENOME_VERSION}
  GENOME=${GENOME_TYPE}
  SOURCE=gencode
  GENOME_ANNOTATION_FOLDER=${INPUT_FOLDER}/data/${GENOME}-${SOURCE}-version-${VERSION}
}

create_folders(){
  mkdir -p ${INPUT_FOLDER}
  mkdir -p ${INPUT_FOLDER}/data
  mkdir -p ${INPUT_FOLDER}/analysis
  mkdir -p ${INPUT_FOLDER}/tools
  mkdir -p ${INPUT_FOLDER}/reads
}


download_genome(){
  mkdir -p ${GENOME_ANNOTATION_FOLDER}
  if [ ${GENOME} == "mouse" ]
  then
    echo ''
    echo 'Downloading Mouse genome and annotation'
    echo ''
    SOURCE_MOUSE=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${VERSION}/
    if test ! -e ${GENOME_ANNOTATION_FOLDER}/GRCm38.p6.genome.fa.gz; then
		 echo "test NOT working"
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/GRCm38.p6.genome.fa.gz
	 fi
    if test ! -e  ${GENOME_ANNOTATION_FOLDER}/gencode.vM24.annotation.gff3.gz; then
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/gencode.vM24.annotation.gff3.gz
	 fi
    if test ! -e ${GENOME_ANNOTATION_FOLDER}/gencode.vM24.annotation.gtf.gz; then
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/gencode.vM24.annotation.gtf.gz 
	 fi
    if test ! -e ${GENOME_ANNOTATION_FOLDER}/gencode.vM24.2wayconspseudos.gff3.gz; then
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/gencode.vM24.2wayconspseudos.gff3.gz
	 fi
    if test ! -e ${GENOME_ANNOTATION_FOLDER}/gencode.vM24.long_noncoding_RNAs.gff3.gz; then
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/gencode.vM24.long_noncoding_RNAs.gff3.gz
	 fi
    if test ! -e ${GENOME_ANNOTATION_FOLDER}/gencode.vM24.tRNAs.gff3.gz; then
		 wget -cP ${GENOME_ANNOTATION_FOLDER} \
		 ${SOURCE_MOUSE}/gencode.vM24.tRNAs.gff3.gz
	 fi
  elif [ ${GENOME} == "human" ]
  then
    echo "Downloading Human genome and annotaion \n"
    SOURCE_HUMAN=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/
    wget -cP ${GENOME_ANNOTATION_FOLDER} \
    ${SOURCE_HUMAN}/GRCh38.p13.genome.fa.gz \
    ${SOURCE_HUMAN}/gencode.v33.annotation.gff3.gz \
    ${SOURCE_HUMAN}/gencode.v33.annotation.gtf.gz \
    ${SOURCE_HUMAN}/gencode.v33.2wayconspseudos.gff3.gz \
    ${SOURCE_HUMAN}/gencode.v33.long_noncoding_RNAs.gff3.gz \
    ${SOURCE_HUMAN}/gencode.v33.tRNAs.gff3.gz
  else
    echo "\n Error! Choose either human or mouse genome \n"
  fi
}


calculate_md5sum(){
  md5sum ${GENOME_ANNOTATION_FOLDER}/* > ${GENOME_ANNOTATION_FOLDER}/genome_annoatation_md5sum.txt
}


decompress_genome_annotation_files(){
  for FILE in $(ls $GENOME_ANNOTATION_FOLDER)
    do
      #NAME=$(echo $FILE | sed "s/.gz//")
      #zcat ${INPUT_FOLDER}/${FILE} > ${OUTPUT_FOLDER}/${NAME}
      echo $GENOME_ANNOTATION_FOLDER/${FILE}
      gunzip -k $GENOME_ANNOTATION_FOLDER/${FILE}
    done
}

make_genome_annotation_unwritetable(){
 chmod -R ugo-w $GENOME_ANNOTATION_FOLDER
}

main

