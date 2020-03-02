#Download rawreads from SRA and perform quality control
#!/bin/bash
#================================================================
#% RUN 
#+    ./download_install_tools.sh
#%
#% DESCRIPTION
#%    This is a script will download and install all the tools
#%    required for Reference based RNAseq analysis.
#%    In order run the analysis in a specified folder and path,
#%    then set the names in the project_config.sh
#%		(Default:RNAseq_project)and "FOLDER_PATH"(Default:current working directory).
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
	check_and_install
}

import_project_config(){
. ${PWD}/project_config.sh
}

check_and_install(){
  install_sratoolkit
  install_fastqc
  install_cutadapt
  install_star
  install_bowtie
  install_tophat
  install_samtools
  install_rseqqc
  install_bedtools
  install_htseq
  install_cufflinks
  install_dexseq
}

set_path(){
  FOLDER_NAME=${PROJECT_NAME}
  FOLDER_PATH=${PROJECT_PATH}
  INPUT_FOLDER=${FOLDER_PATH}/${FOLDER_NAME}/
}


create_folders(){
  mkdir -p ${INPUT_FOLDER}
  mkdir -p ${INPUT_FOLDER}/data
  mkdir -p ${INPUT_FOLDER}/analysis
  mkdir -p ${INPUT_FOLDER}/tools
  mkdir -p ${INPUT_FOLDER}/reads
}


install_sratoolkit(){
  cd ${INPUT_FOLDER}/tools
  #Download SRAtoolkit
  echo "Checking and downloading SRAtoolkit "
  if [ -d "${INPUT_FOLDER}/tools/sratoolkit/" ]; then
     echo "SRAtoolkit already installed "
  else
     echo "Installing SRAtoolkit"
     wget -cP ${INPUT_FOLDER}/tools/ \
     http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/sratoolkit.2.10.2-ubuntu64.tar.gz
     tar -xzf sratoolkit.2.10.2-ubuntu64.tar.gz
     rm sratoolkit.2.10.2-ubuntu64.tar.gz
     mv sratoolkit.2.10.2-ubuntu64 sratoolkit
  fi
  echo "DONE checking and downloading SRAtoolkit!"
}

install_fastqc(){
  cd ${INPUT_FOLDER}/tools
  echo "Checking and installing FastQC"
  if [ -d "${INPUT_FOLDER}/tools/FastQC" ]; then
     echo "FastQC already installed"
  else
     echo "Installing FastQC"
     wget -cP ${INPUT_FOLDER}/tools \
     https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
     unzip fastqc_v0.11.8.zip
     cd FastQC/
     chmod 777 fastqc
     cd ..
     rm fastqc_v0.11.8.zip
  fi
  echo "DONE checking and downloading FastQC!"
}

install_cutadapt(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing cutadapt \n"
  if [ -d "${INPUT_FOLDER}/tools/cutadapt" ]; then
     echo "\n cutadapt already installed \n"
  else
    echo "\n Installing cutadapt \n"
    pip3 install \
    --user \
    --install-option="--install-scripts=${INPUT_FOLDER}/tools/cutadapt/bin" \
    cutadapt==1.9.1
  fi
  echo "\n DONE checking and downloading cutadapt!\n"
}


install_star(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing STAR \n"
  if [ -d "${INPUT_FOLDER}/tools/STAR-2.7.3a" ]; then
     echo "\n STAR already installed\n"
  else
    echo "\n Installing STAR \n"
    wget -cP ${INPUT_FOLDER}/tools \
    https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
    tar -xvzf 2.7.3a.tar.gz
    rm 2.7.3a.tar.gz
  fi
  echo "\n DONE checking and downloading STAR!\n"
}


install_bowtie(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing Bowtie \n"
  if [ -d "${INPUT_FOLDER}/tools/bowtie2" ]; then
     echo "\n Bowtie already installed \n"
  else
    echo "\nInstalling Bowtie \n"
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
	 unzip bowtie2-2.3.5.1-linux-x86_64.zip
	 mv bowtie2-2.3.5.1-linux-x86_64 bowtie2
    rm bowtie2-2.3.5.1-linux-x86_64.zip
 fi
  export PATH=${INPUT_FOLDER}/tools/bowtie:$PATH
  echo "\n DONE checking and downloading Bowtie!\n"
}

#Download TopHat
install_tophat(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing TopHat \n"
  if [ -d "${INPUT_FOLDER}/tools/tophat" ]; then
     echo "\nTopHat already installed\n"
  else
    echo "\nInstalling TopHat \n"
    wget -cP ${INPUT_FOLDER}/tools \
	 https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz 
    tar -xvzf tophat-2.1.1.Linux_x86_64.tar.gz
	 mv tophat-2.1.1.Linux_x86_64 tophat
	 rm tophat-2.1.1.Linux_x86_64.tar.gz
  fi
   echo "\n DONE checking and downloading TopHat!\n"
}

install_samtools(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing SAMtools \n"
  if [ -d "${INPUT_FOLDER}/tools/samtools" ]; then
     echo "\n SAMtools already installed \n"
  else
    echo "\n Installing SAMtools \n"
    wget -cP ${INPUT_FOLDER}/tools \
	 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 
	 tar -xjf samtools-1.9.tar.bz2
	 mv samtools-1.9 samtools
	 cd samtools
	 ./configure --prefix=${INPUT_FOLDER}/tools/samtools/
	 make 
	 make install
	 cd ..
	 rm samtools-1.9.tar.bz2
  fi
  echo "DONE checking and downloading SAMtools!"
}

#export PATH=${TOOLS_FOLDER}/samtools-1.9:$PATH
#Download RseqQC
install_rseqqc(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing RseqQC \n"
  if [ -d "${INPUT_FOLDER}/tools/rseqqc" ]; then
     echo "\n RseqQC already installed\n"
   else
     echo "\n Installing RseqQC \n"
     pip3 install \
     --user \
     --install-option="--install-scripts=${INPUT_FOLDER}/tools/rseqqc/bin" \
     Cython


     pip3 install \
     --user \
     --install-option="--install-scripts=${INPUT_FOLDER}/tools/rseqqc/bin" \
     RSeQC
   fi
   echo "\n DONE checking and downloading RseqQC!\n"
}


install_bedtools(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing bedtools \n"
  if [ -d "${INPUT_FOLDER}/tools/bedtools2" ]; then
     echo "\n bedtools already installed\n"
   else
     echo "\n Installing bedtools \n"
     wget -cP ${INPUT_FOLDER}/tools \
     https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
     tar -zxvf bedtools-2.29.1.tar.gz
     cd bedtools2
     make
     cd ..
     rm bedtools-2.29.1.tar.gz
   fi
   echo "\n DONE checking and downloading bedtools!\n"
}


install_htseq(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing HTSeq \n"
  if [ -d "${INPUT_FOLDER}/tools/htseq" ]; then
     echo "\n HTSeq already installed\n"
   else
     echo "Installing HTSeq"
     pip3 install \
    --user \
    --install-option="--install-scripts=${INPUT_FOLDER}/tools/htseq/bin" \
     htseq
     #pip install 'matplotlib>=1.4'
     #pip install Cython
     #pip install 'pysam>=0.9'
     #git clone https://github.com/simon-anders/htseq.git
     #cd htseq
     #python3 setup.py --prefix=${INPUT_FOLDER}/tools/htseq build 
     #python3 setup.py install
     #cd ..
     #pip install HTSeq
   fi
   echo "\n DONE checking and downloading HTSeq!\n"
}


#Download cufflinks
install_cufflinks(){
  cd ${INPUT_FOLDER}/tools
  echo "\n Checking and installing cufflinks \n"
  if [ -d "${INPUT_FOLDER}/tools/cufflinks" ]; then
     echo "\n cufflinks already installed\n"
   else
     wget -cP ${INPUT_FOLDER}/tools \
     http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
     tar xzvf cufflinks-2.2.1.Linux_x86_64.tar.gz
     mv cufflinks-2.2.1.Linux_x86_64/ cufflinks
     rm cufflinks-2.2.1.Linux_x86_64.tar.gz     
  fi
  echo "\n DONE checking and downloading cufflinks!\n"
}


#Install dexseq
install_dexseq(){
  echo "\n Checking and installing R package DEXSeq\n"
  cd ${INPUT_FOLDER}/tools
  git clone https://github.com/RB786/DEXSeq.git
  echo "\nDONE checking and downloading DEXSeq!\n"
}


main

