![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg)[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.6.1-blue.svg)](https://cran.r-project.org/)[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# Design and Analysis of RNA Sequencing Data

This is the associated GitHub page for the book chapter *Design and Analysis of RNA Sequencing Data* *Bharti and Grimm et. al. 2020*. A comprehensive outline of RNASeq data analysis including quality assessment, pre-processing alignment/mapping of reads, quantification of genes and isoforms and differential expression (DE) analysis is provided. Step wise optimized codes and data formats alongside brief description and references to the corresponding sections of the book chapter are also provided.

### OS
Any Linux based distro should work. Out test OS is:

Distributor ID: Ubuntu <br/>
Description:    Ubuntu 18.04.2 LTS <br/>
Release: 	18.04 <br/>
Codename:       bionic <br/>

`lsb_release -a` on a Ubuntu based system.

### Hardware
<p style='text-align: justify;'> It is no secret that the hardware highly influences the speed of the workflows. And just one same can generate ~10-60 Million reads. The most time consuming tasks are the ones that involve reference based alignment. A modest configuration consists of 16+cores and 16 GB of RAM with TB external hard drive or 48 TB server and 3.6 GHz clock speed if one wants to run it on one’s own workstation. A majority of the diskspace is occupied by reference . Our HW configuration consists of 20 core CPU with 128 GB.</p>

### Software and packages
<p style='text-align: justify;'> All software should be installed by the user directly as the workflow depends on a lot of external software.
Without these the workflow will fail to run. </p>

- gcc, g++

- java

- python 2 python3: pip

- R version 3.6.2

- git



# Reference based analysis:

### Get the workflow from git using:  

```bash
git clone https://github.com/RB786/Chapter_RNASeq-analysis.git
```

After cloning the repository, you should have the following files in your repo directory:

```bash
$ls -1
Rscripts
differential_expression.sh
download_genome_and_annotation.sh
download_install_tools.sh
download_RNASeq_reads_and_QC.sh
metadata.txt
project_config.sh
quantification.sh
reference_based_alignment.sh
```

### Workflow organization:

The workflow centered around two variables which contain the **project name** and the **project path**. The **project name** is the root directory of the project under which the project relevant files are stored. The **project path** is the path on the system where the project root directory is created.  These two variables can be set in a configuration file called as `project_conf.sh`. The default values of the project_conf.sh are shown here:

```bash
PROJECT_NAME=RNASeq_project
PROJECT_PATH=${PWD}
```

Under the defined project name root directory, four sub-directories are created. They are:

```bash
|-- RNASeq_project
|   | analysis
|   | data
|   | reads
|   | tools
```

- *analysis*: This will contain the output of the analysis performed during each step of the workflow.
- *data*: This will contain the genome, annotation and index genome files.
- *reads*: This will contain the sequencing reads which is used as an input to perform the analysis.
- *tools*: This will contain all the tools and software required to run the workflow.



## Brief Description of the Each Step

1. #### **Download and install tools and python libraries (download_install_tools.sh):**

   The first step for running the workflow is to download and install all the required tools used by the workflow. This can be done by running the shell script:

   ```bash
   $./download_install_tools.sh
   ```

   This will first create a working folder structure as shown in the above section. This script will download all the tools required to perform the analysis inside the **tools** folder as listed below:

   ```bash
   ~/RNASeq_project/tools$ls -1
   DEXSeq
   FastQC
   STAR-2.7.3a
   bedtools2
   bowtie2
   cufflinks
   cutadapt
   htseq
   rseqqc
   samtools
   sratoolkit
   tophat
   ```



2. #### **Download Eukaryotic genome and its annotation (download_genome_and_annotation.sh)**
   Assumption here is that RNASeq data is from mouse or human. The genomes of these species are available and well annotated. The genome file is in fasta format and the annotation file is in gtf/gff format. This script will download the current [gencode genome](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/) (default: mouse) and its annotation into the **data** folder. Currently, the latest version of gencode is used for the mouse genome. The version can be changed by modifying `GENOME_VERSION` parameter and the genome can be changed to human by modifying `GENOME_TYPE` in the `project_conf.sh` file:

   ```bash
   GENOME_TYPE=mouse # or human
   GENOME_VERSION=24 #or 33 for human
   ```

   This can be done by running the shell script:

   ```bash
   $./download_genome_and_annotation.sh
   ```

   After the execution of the script, the data folder will contain:

   ```bash
   ~/RNASeq_project/data$ls -1
   mouse-gencode-version-24
   ```

   *Output:*

   ```bash
   ~/RNASeq_project/data/mouse-gencode-version-24$ls -1
   GRCm38.p6.genome.fa
   GRCm38.p6.genome.fa.gz
   gencode.vM24.2wayconspseudos.gff3
   gencode.vM24.2wayconspseudos.gff3.gz
   gencode.vM24.annotation.gff3
   gencode.vM24.annotation.gff3.gz
   gencode.vM24.annotation.gtf
   gencode.vM24.annotation.gtf.gz
   gencode.vM24.long_noncoding_RNAs.gff3
   gencode.vM24.long_noncoding_RNAs.gff3.gz
   gencode.vM24.tRNAs.gff3
   gencode.vM24.tRNAs.gff3.gz
   genome_annoatation_md5sum.txt
   ```

3. #### **Download the data and perform quality control (download_RNASeq_reads_and_QC.sh)**

   This step is used to download an example data to perform the RNASeq analysis. This can be done by running the shell script:

   ```bash
   $./download_RNASeq_reads_and_QC.sh
   ```

   This script contains three main sub-steps:
   
   __*a. Downloading the publically available data*__
   
   Publically available data from PMID: [28738885](https://www.ncbi.nlm.nih.gov/pubmed/28738885) is used for the analysis purpose. In the first step data is downloaded using `fastq-dump` binary from SRA toolkit in the **tools** folder.

   The `fastq-dump` tool uses the following parameters:

   ```bash
   fastq-dump \
      ${SRA_ID} \
      --split-files \
      --origfmt \
      --gzip
   ```
   
   *Description:*
   
   `fastq-dump` is the binary

   `${SRA_ID}` are list of fastq files. These SRA  files will be downloaded one at a time

   `--split-files` parameter will produce two reads if the data is paired end else single file

   `--origfmt`  parameter contains only original sequence name

   `--gzip` will provide the compressed output using gzip

   *Output:*
   The `fastq-dump` downloads the following files into the **reads/rawreads** folder:

   ```bash
   ~/RNASeq_project/reads/rawreads$ls -1
   SRR5858228_1.fastq.gz
   SRR5858229_1.fastq.gz
   SRR5858230_1.fastq.gz
   SRR5858231_1.fastq.gz
   SRR5858232_1.fastq.gz
   SRR5858233_1.fastq.gz
   SRR5858234_1.fastq.gz
   SRR5858235_1.fastq.gz
   SRR5858236_1.fastq.gz
   ```

   *NOTE! To read more about fastq-dump parameters click [here](https://ncbi.github.io/sra-tools/fastq-dump.html).*



   __*b. Quality assessment of downloaded data:*__

   The downloaded SRA fastq files are then assessed for sequence quality using `fastqc` program.

   The `fastqc` tool uses the following parameters:

   ```bash
   fastqc \
      SRR5858228_1.fastq.gz \
      -o rawreads_QA_stats
   ```

   *Description:*
   
   `fastqc` is the binary.
   
   `SRR5858229_1.fasta.gz`  is one example file as an input for quality assessment. This is run in a for loop which will run on all downloaded samples.
   
   `-o` is path to output directory, here it is `reads_QA_stats` .

   *Output:*
   The `fastqc` will generate the output for all downloaded files into the **reads/rawreads_QA_stats** folder:

   ```bash
   ~/RNASeq_project/reads/rawreads_QA_stats$ls -1
   SRR5858228_1_fastqc.html
   SRR5858228_1_fastqc.zip
   SRR5858229_1_fastqc.html
   SRR5858229_1_fastqc.zip
   SRR5858230_1_fastqc.html
   SRR5858230_1_fastqc.zip
   SRR5858231_1_fastqc.html
   SRR5858231_1_fastqc.zip
   SRR5858232_1_fastqc.html
   SRR5858232_1_fastqc.zip
   SRR5858233_1_fastqc.html
   SRR5858233_1_fastqc.zip
   SRR5858234_1_fastqc.html
   SRR5858234_1_fastqc.zip
   SRR5858235_1_fastqc.html
   SRR5858235_1_fastqc.zip
   SRR5858236_1_fastqc.html
   SRR5858236_1_fastqc.zip
   ```

   *Note! For more details about fastqc, read section XXX in book and for parameters check [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*
      

   __*c. Quality control of downloaded data:*__

   Based on the assessment of the report quality control (QC) is performed using the  `cutadapt` tool.

   The `cutadapt` tool uses the following parameters:

   ```bash
   cutadapt \
      -q 20 \
      --minimum-length 35 \
      -a "CTGTCTCTTATACACATCT" \
      -o reads/filtered_reads/SRR5858228_1_trimmed.gz \
      SRR5858229_1.fastq.gz \ #Input file name
      > SRR5858229_1_cutadapt_stats.txt # Output stats for given file
   ```

   *Description:*

   `-a`  is the adapter sequence for trimming

   `--minimum-length`  is the minimum sequence length (35)

   `-q` phred score of 20

   `-o` is the output file name

   *Output:*
   The `cutadapt`  will generated trimmed files into the **reads/filtered_reads** folder:

   ```bash
   ~/RNASeq_project/reads/filtered_reads$ls -1
   SRR5858228_1_trimmed.gz
   SRR5858229_1_trimmed.gz
   SRR5858230_1_trimmed.gz
   SRR5858231_1_trimmed.gz
   SRR5858232_1_trimmed.gz
   SRR5858233_1_trimmed.gz
   SRR5858234_1_trimmed.gz
   SRR5858235_1_trimmed.gz
   SRR5858236_1_trimmed.gz
   ```

   *Note! Details of quality control and how to choose the parameter are chosen is explained in chapter (XXX). To understand parameter in detail check [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html).*



4. #### **Perform reference alignment using genome sequence (reference_based_alignment.sh)**

   Aligning reads to genome can be challenging due to millions of reads, small read size (2nd generation sequencing), sequencing errors and variation leads to mismatches and indels. In addition to that eukaryotic genome have introns , this means when genome is used a reference, aligner is aware of splice junction in the reference.

   This step will perform reference based RNASeq alignment. This can be done by running the shell script:

   ```bash
   $./reference_based_alignment.sh
   ```

   Currently, several splice aware tools are available. In this workflow we have used the two most commonly used tools TopHat2 and STAR. The choice of `ALIGNER` can be set in the `project_conf.sh` file:

   ```bash
   ALIGNER=star #tophat or star
   ```

   The table below shows the tools used for each alignment step for both the STAR and TopHat2 aligners.

   | Alignment Steps           | STAR Alignment                                               | TopHat2 Alignment                                            |
   | :------------------------ | :----------------------------------------------------------- | :----------------------------------------------------------- |
   | Genome indexing           | *[STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)* (a.) | [*bowtie2*](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (c.) |
   | Reference based alignment | *[STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)* (b.) | [*TopHat2*](http://ccb.jhu.edu/software/tophat/manual.shtml) (d.) |
   | BAM file indexing         | [*samtools*](http://samtools.sourceforge.net/) (e.)          | [*samtools*](http://samtools.sourceforge.net/)               |
   | Alignment stats           | [*RseQC*](http://rseqc.sourceforge.net/) (f.)                | [*RseQC*](http://rseqc.sourceforge.net/)                     |
   | Visualization             | *IGV* (g.)                                                   | *IGV*                                                     |
   
      __*a. Genome indexing using STAR:*__
   
      ```bash
     STAR \
      --runThreadN 8 \
      --runMode genomeGenerate \
      --genomeDir data/STAR_Genome_Index \
      --genomeFastaFiles data/mouse-gencode-version-24/GRCm38.p6.genome.fa \
      --sjdbGTFfile data/mouse-gencode-version-24/gencode.vM24.annotation.gtf
     ```
    
     *Description:*

     `STAR` is the binary

     `--runMode` is set to `genomeGenerate`  which generate genome files

     `--genomeDir` is the path to output directory name

     `--genomeFastaFiles` is the path to genome fasta file with file name

     `--sjdbGTFfile` is the path to the annotation file with file name

     `--runThreadN` is number of threads (default used is 8)

     *Output:*
     The `STAR`  genome indexing will generated indexed genome files into the **data/STAR_Genome_Index** folder:

     ```bash
     ~/RNASeq_project/data/STAR_Genome_Index$ls -1
     Genome
     SA
     SAindex
     chrLength.txt
     chrName.txt
     chrNameLength.txt
     chrStart.txt
     exonGeTrInfo.tab
     exonInfo.tab
     geneInfo.tab
     genomeParameters.txt
     sjdbInfo.txt
     sjdbList.fromGTF.out.tab
     sjdbList.out.tab
     transcriptInfo.tab
     ```
     *Note! For more detail read section XXXX in the book* and for parameters check *[STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)*
     
     __*b. Alignment using STAR:*__
  
     ```bash
     STAR \
      --runMode alignReads \
      --runThreadN 8 \
      --genomeDir data/STAR_Genome_Index \
      --readFilesIn reads/filtered_reads/SRR5858228_1_trimmed.gz \
      --readFilesCommand zcat \
      --sjdbGTFfile data/mouse-gencode-version-24/gencode.vM24.annotation.gtf \
      --outFileNamePrefix /analysis/mapping/star/SRR5858228_1_trimmed  \
      --outSAMtype BAM SortedByCoordinate
      ```
      
      *Description:*

      `STAR` is the binary

      `--runMode`  is `alignReads`  which means map reads to the reference genome

      `--runThreadN` is the number of threads (default used is 8)

      `--genomeDir` is the path to genome indexed directory

      `--readFilesIn` is the path to trimmed reads including file name

      `--readFilesCommand` for gzipped files (*.gz) use zcat

      `--sjdbGTFfile` is the path annotation file with file name

      `--outFileNamePrefix`  is output prefix name with its path

      `--outSAMtype` is the output sorted by coordinate, similar to samtools sort command

      *Output:*
      The `STAR`  genome alignment will generated in files into the **analysis/mapping/star** folder:

      This is the output for one sample, similar files will be generated for other samples as well

      ```bash
      ~/RNASeq_project/analysis/mapping/star$ls -1
      Percentage_uniquely_mapped_reads.csv
      SRR5858228_1_trimmed.Aligned.sortedByCoord.out.bam
      SRR5858228_1_trimmed.Aligned.sortedByCoord.out.bam.bai
      SRR5858228_1_trimmed.Log.final.out
      SRR5858228_1_trimmed.Log.out
      SRR5858228_1_trimmed.Log.progress.out
      SRR5858228_1_trimmed.SJ.out.tab
      SRR5858228_1_trimmed._STARgenome
      ```

      - **SRR5858228_1_trimmed.Aligned.sortedByCoord.out.bam** is aligned bam file
      - **SRR5858228_1_trimmed.SJ.out.tab**  is a tab-delimited file that provides information about alignments to splice junctions
      - **Percentage_uniquely_mapped_reads.csv**  file represents the uniquely mapped reads

      These three files are as names suggests provides the information about ongoing samples alignment. Out of these file with Log.final.out extension provides the complete mapping stats 	

      - **SRR5858228_1_trimmed.Log.final.out**
      - **SRR5858228_1_trimmed.Log.out**
      - **SRR5858228_1_trimmed.Log.progress.out**

      *Note! For more detail read section XXXX in the book* and for parameters check *[STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)*
      
      
      __*c. Genome indexing using bowtie for TopHat2:*__

      TopHat2 uses bowtie2 (which is gapped aligner) for genome indexing

      ```bash
      bowtie2-build \
         -f data/mouse-gencode-version-24/GRCm38.p6.genome.fa \
         /data/TOPHAT_Genome_Index/ \
         -p 8
      ```

      *Description:*

      `bowtie2-build` is the binary

      `-f` is the path to genome fasta file with file name

      `-p` is to launch a specified number of parallel search threads

      `TOPHAT_Genome_Index` is the output directory created in **data** folder

      *Output:*
      The `bowtie-build`  genome indexing will generated indexed genome files into the **data/TOPHAT_Genome_Index** folder:

      ```bash
      ~/RNASeq_project/data/TOPHAT_Genome_Index$ls -1
      GRCm38.p6.genome.1.bt2
      GRCm38.p6.genome.2.bt2
      GRCm38.p6.genome.3.bt2
      GRCm38.p6.genome.4.bt2
      GRCm38.p6.genome.fa
      GRCm38.p6.genome.rev.1.bt2
      GRCm38.p6.genome.rev.2.bt2
      ```

      *Note! Bowtie2 itself is not able to perform spliced alignments. Bowtie2 is known for its speed and small-memory efficiency because it index the reference genome using an FM index which is a [Burrows–Wheeler transform](https://en.wikipedia.org/wiki/Burrows–Wheeler_transform) method. For more detail read section XXXX in the book.*
      
      __*d. Alignment using TopHat2:*__

      ```bash
      tophat \
         -o /analysis/mapping/tophat \
         -p 8 \
         -G data/mouse-gencode-version-24/gencode.vM24.annotation.gtf \
         /data/TOPHAT_Genome_Index/GRCm38.p6.genome \
         reads/filtered_reads/SRR5858228_1_trimmed.gz
      ```

      *Description:*

      `tophat` is the binary

      `-o`  is the name of the directory in which TopHat will write all of its outputs

      `-p` is the number of threads  to align reads (default used is 8)

      `-G` is the path to the annotation file with file name

      `/data/TOPHAT_Genome_Index/GRCm38.p6.genome` is the path and prefix name of genome fasta file

      `reads/filtered_reads/SRR5858228_1_trimmed.gz` is the path to trimmed reads including file name


      *Output:*
      The `tophat`  genome alignment will generated in files into the **analysis/mapping/tophat** folder

      TopHat generated many files as listed below. Here is example output for one SRA ID SRR5858228:

      - **SRR5858228_1_trimmed_accepted_hits.bam** represents the alignment file in bam format.
      - **SRR5858228_1_trimmed_junctions.bed** consists of a list of exon junctions in BED formation. An exon junction comprises two blocks where either block is as long as the longest overhang of any read present in the junction sequence. The score is the number of alignments identified for a junction sequence. 
      - **SRR5858228_1_trimmed_insertions.bed** consists of a list of discovered insertions. In each case chromLeft represents the last genomic base before individual insertions.
      - **SRR5858228_1_trimmed_deletions.bed** consists of a list of discovered deletions. In each case chromLeft represents the last genomic base before individual deletions.
      - **SRR5858228_1_trimmed_align_summary.txt** contains a summary of alignment rates along with the number of  reads and pairs showing multiple alignments.

      *Note! For more detail read section XXXX in the book* and for parameters check [*TopHat2*](http://ccb.jhu.edu/software/tophat/manual.shtml).
      
      
      __*e. BAM file indexing for both STAR or TopHat2:*__

      ```bash
      samtools \
         index analysis/mapping/star/SRR5858228_1_trimmed.bam
      ```

      *Description:*

      `samtools ` is the binary

      `index`  will index a coordinate-sorted BAM file for fast random access


      *Output:*
      This command will index all the bam files and create `.bai` files in the respective **mapping** folder.

      *Note! For more detail read section XXXX in the book* and for parameters check [*samtools*](http://samtools.sourceforge.net/).
      
      
      __*f. Generate stats for alignment:*__

      ```bash
      python3 bam_stat.py \
         -q 30 \
         -i analysis/mapping/star/SRR5858228_1_trimmed.bam >  \
         analysis/mapping/star/SRR5858228_1_trimmed_RSeQC_alignment_stats.txt
      ```

      *Description:*

      `bam_stat.py ` is the python script to calculate bam stats

      `-i`  input bam file with its path

      `-q`  is the mapping quality to determine uniquely mapped read

      RseQC produces table in the respective **mapping** folder where unique reads are considered if their mapping quality is more than 30.

      *Output:* 
      **SRR5858228_1_trimmed.Aligned.sortedByCoord.out._RSeQC_alignment_stats.txt** file will be created for SRA sample ID SRR5858228.

      |   #============== # All numbers are READ count  #============== |                            |
      | :--------------------------------------------- | :------------------------------------------ |
      |Total records:                                  |             	              72432644       |
      |                                                |                                             |
      | QC failed:                                     |                                          0  |
      | Optical/PCR duplicate:                         |                                          0  |
      | Non primary hits                               |                                   11660623  |
      | Unmapped reads:                                |                                           0 |
      |mapq < mapq_cut (non-unique):                   |                                    5935125  |
      |                                                |                                             |
      | mapq >= mapq_cut (unique):                     |          					          54836896  |
      | Read-1:                                        |                            				  0  |
      |Read-2:                                         |                            				  0  |
      | Reads map to '+':                              |                						  27245688 |
      | Reads map to '-':                              |              							 27591208  |
      |  Non-splice reads:                             |           								  42688777 |
      | Splice reads:                                  |        									  12148119 |
      |  Reads mapped in proper pairs:				       |			                                   0  |
      |  Proper-paired reads map to different chrom:	 |	                                         0  |


      *Note! For more detail read section XXXX in the book* and for parameters check [*RseQC*](http://rseqc.sourceforge.net/).
      
      
      __*g. Visualization of alignment (BAM) files: *__
      
      The indexed bam files along with genome and annotation can be loaded into several genome browsers, including the Integrative Genomics Viewer IGV `IGV` . Explaining this is beyond the scope of this chapter and recommend to go through https://software.broadinstitute.org/software/igv/UserGuide.



5. #### **Quantification (quantification.sh)**

   This step will perform quantification on alignment files generated form reference based RNASeq alignment. This can be done by running the shell script:

   ```bash
   $./quantification.sh
   ```

   The three ways to perform quantification is XXXXXX aware tools are available. In this workflow we have used the two most commonly used tools XXXXX. The choice of `QUANT_COUNT` can be set in the `project_conf.sh` file:

   ```bash
   QUANT_COUNT=per_gene_bedtools # per_gene_htseq or per_transcript_cufflinks or per_exon_dexseq
   ```

   | Quantification methods (Counting reads) | Tools                                                        |
   | :-------------------------------------- | :----------------------------------------------------------- |
   | per gene                                | [bedtools](https://bedtools.readthedocs.io/en/latest/) (a.)  |
   | per gene                                | [HTSeq](https://htseq.readthedocs.io/) (b.)                  |
   | per transcript                          | *[Cufflinks](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjBgZnNsPnnAhWJ5KQKHX-WDUEQFjAAegQIBBAB&url=http%3A%2F%2Fcole-trapnell-lab.github.io%2Fcufflinks%2Fcufflinks%2F&usg=AOvVaw1JN_NkATlDD-FIYoGiXJFf)* (c.) |
   | per exon                                | *[DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html)* (d.) |
   
   
   __*a. Counting reads per gene using *bedtools*:*__

   ```bash
   bedtools intersect \
      -S \
      -wa \
      -c \
      -a genes_only.gff3 \
      -b analysis/mapping/star/SRR5858228_1_trimmed.bam \
      > analysis/quantification/bedtools-count/SRR5858228_1_trimmed_counts.csv
   ```

   *Description:*
   
   `bedtools` is the binary
   
   `intersect` allows one to screen for overlaps between two sets of genomic features
   
   `-S`  is require for different strandedness
   
   `-c`  is for each entry in genomic features A, report the number of hits in genomic features B while restricting to `-f` (-f is  for minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp))
   
   `-a`  is gff annotation file where each genomic features in A is compared to genomic features B in search of overlaps
   
   `-b` is the input bam file with its path

   *Output:*
   Output will be generated in **analysis/quantification/bedtools-count** folder for each file. In the end all the  bedtools generated output files are merged to created final file **Read_per_features_combined.csv** in same folder which will be the input for differential expression analysis.

   *Note! For more details about bedtools read section XXX in the book and for parameters check [bedtools](https://bedtools.readthedocs.io/en/latest/)*


   __*b. Counting reads per gene using *HTSeq*:*__

   ```bash
   htseq-count \
      -f bam \
      -a 10 \
      -m intersect-strict \
      -s no \
      -t exon \
      -i gene_id \
      analysis/mapping/star/SRR5858228_1_trimmed.bam \
      data/mouse-gencode-version-24/gencode.vM24.annotation.gtf > \
      analysis/quantification/htseq-count/SRR5858228_1_trimmed_counts.csv
   ```

   *Description:*
   
   `htseq-count` is the binary
   
   `-f`  is the format of the input data.
   
   `-a `  will skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10).
   
   `-m` is mode to handle reads overlapping more than one feature. In this case its `intersect-strict`.
   
   `-s` is set whether the data is from a strand-specific assay. For `stranded=no`, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature
   
   `-t ` is the feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an [GENCODE GTF](https://www.gencodegenes.org/) file: `exon`)
   
   `-i`  is the attribute to be used as feature ID. The default, suitable for RNA-Seq analysis using an GENCODE GTF file, is `gene_id`.

   *Output:*
   The outputs will be generated in **analysis/quantification/htseq-count folder**. In the end all the  HTSeq generated output files are combined to created final file **Read_per_features_combined.csv** in same folder which will be used as an input for differential expression analysis.

   *Note! For more details about HTSeq read section XXX in the book and for parameters check [HTSeq](https://htseq.readthedocs.io/en/release_0.9.1/count.html#count)*


   __*c. Counting reads per transcript using *cufflinks*:*__

   ```bash
   cufflinks \
      -G data/mouse-gencode-version-24/gencode.vM24.annotation.gtf  \
      -b data/mouse-gencode-version-24/GRCm38.p6.genome.fa \
      -p 8 \
      analysis/mapping/star/SRR5858228_1_trimmed.bam \
      -o analysis/quantification/cufflinks-count
   ```

   *Description:*
   
   `cufflinks` is the binary
   
   `-G` is the path to the annotation file with file name
   
   `-b` is the path to genome fasta file with file name
   
   `-p` is the number threads
   
   `-o` is the name of the directory in which Cuffdiff will write all of its output

   *Output:*
   The outputs will be generated in **analysis/quantification/cufflinks-count folder**. The output consists of transcript and gene-level FPKM-tracking files, which contain FPKM values and their confidence intervals. FPKM tracking files are also produced when a set of samples is tested for differential expression using Cuffdiff.

   *Note! For more details about Cufflinks read section XXX in the book and for parameters check [Cufflinks](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjBgZnNsPnnAhWJ5KQKHX-WDUEQFjAAegQIBBAB&url=http%3A%2F%2Fcole-trapnell-lab.github.io%2Fcufflinks%2Fcufflinks%2F&usg=AOvVaw1JN_NkATlDD-FIYoGiXJFf).*
   
   
   __*d. Counting reads per  exon using *DEXSeq*:*__

   The initial steps of a *[DEXSeq](https://bioconductor.org/packages/3.11/DEXSeq)* analysis are done using two Python scripts. Importantly, these preprocessing steps can also be done using tools equivalent to these Python scripts, for example, using *[GenomicRanges](https://bioconductor.org/packages/3.11/GenomicRanges)* infrastructure ([10.2](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#preprocessing-within-r)) or *[Rsubread](https://bioconductor.org/packages/3.11/Rsubread)* ([10.3](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#preprocessing-using-featurecounts)). The following two steps describe how to do this steps using the Python scripts that are provided within *[DEXSeq](https://bioconductor.org/packages/3.11/DEXSeq)*.

   <u>*Preparing the annotation*</u>

   ```bash
   python2 dexseq_prepare_annotation.py \
      data/mouse-gencode-version-24/gencode.vM24.annotation.gtf \
      analysis/quantification/dexseq-count/DEXSEQ_GTF_annotation.gff
   ```

   The Python script *dexseq_prepare_annotation.py* takes an GTF file and translates it into a GFF file with collapsed exon counting bins. For more details see Figure 1 of the *[DEXSeq](https://bioconductor.org/packages/3.11/DEXSeq)* paper (Anders, Reyes, and Huber (2012)) for an illustration.

   *<u>Counting reads</u>*

   ```bash
   python2 dexseq_count.py \
      -p no \
      -s no \
      -r name \
      analysis/quantification/dexseq-count/DEXSEQ_GTF_annotation.gff \
      -f bam \
      -a 10 \
      analysis/mapping/star/SRR5858228_1_trimmed.bam \
      analysis/quantification/dexseq-count/SRR5858228_1_trimmed_exon_counts.csv
   ```

   *Description:*
   
   `-p` is to set the paired end information. If the data is from a paired-end sequencing run, you need to add the option `-p yes`
   
   `-s` If you have used a library preparation protocol that does not preserve  strand information (i.e., reads from a given gene can appear equally  likely on either strand), you need to inform the script by specifying  the option `-s no`
   
   `-r` is to indicate whether your data is sorted by alignment position or by read name
   
   `-f` input reads format
   
   `-a` to specify the minimum alignment quality. All reads with a lower quality than specified (with default `-a 10`) are skipped.

   *Output:*
   The outputs for both the steps will be generated in **analysis/quantification/dexseq-count folder** with files for all the samples.


   *Note! For more details about DEXSeq read section XXX in the book and for parameters check [DEXSeq](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html).*



6. #### **Differential Expression (differential_expression.sh)**

   Differential expression (DE) analysis refers to the identification of genes (or other types of genomic features, such as, transcripts or exons) that are expressed in significantly different quantities in distinct groups of samples, be it biological conditions (drug-treated vs. controls), diseased vs. healthy individuals, different tissues, different stages of development, or something else.

   This step will perform quantification on alignment files generated form reference based RNASeq alignment. This can be done by running the shell script:

   ```bash
   $./differential_expression.sh
   ```

   The choice of `QUANT_COUNT` set in the `project_conf.sh` file determines the differential expression tool. The table below shows the differential expression tool used for the generated count data:

   | Generated count data                | Differential expression tools                                |
   | :---------------------------------- | :----------------------------------------------------------- |
   | per gene for bedtools output        | *[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)* (a.) -> Differentially expressed genes |
   | per gene for HTSeq output           | *[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)* -> Differentially expressed genes |
   | per exon for DEXSeq output          | *[DEXSeq](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html)* (b.) -> Differentially expressed exon |
   | per transcript for cufflinks output | *[cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/)* (c.) -> Differentially expressed isoforms |


   __*a. Differential expression for count data generated using *bedtools/HTSeq*:*__

   This analysis of differential expression (DE) will identify of genes that are expressed in significantly different quantities in distinct groups biological conditions (vehicle_treated vs. drug_treated).

   ```bash
   Rscript DE_deseq_genes_bedtools.R \
      --genecount analysis/quantification/bedtools-count/Read_per_features_combined.csv \
      --metadata metadata.txt \
      --condition condition \
      --outputpath analysis/DE/deseq2/
   ```

   *Description:*
   
   `--genecount` is combined generated from bedtools or htseq along with its path
   
   `--metadata` is the sample information file
   
   `--condition` is the column in the metadata data file that will be used for pairwise comparision using DESeq2.
   
   `--outputpath` is the path to the output directory where the result table and plots will be generated.
   
   *Output:*
   Differentially expressed outfiles for bedtools and HTSeq are generated in **analysis/DE/deseq2/**

   ```bash
   ~/RNASeq_project/analysis/DE/deseq2$ls -1
   PCA-samples_bedtools.pdf
   diffexpr-maplot_bedtools.pdf
   pvalue_histogram_50_bedtools.pdf
   qc-dispersions_bedtools.pdf
   qc-heatmap-samples_bedtools.pdf
   vehicle_treated_vs_drug_treated_diffexpr-results_bedtools.csv
   ```

   - **vehicle_treated_vs_drug_treated_diffexpr-results_bedtools.csv** is the final table from deseq2 for the pairwise comparison between the two conditions with stats. The counts are normalized to rlog normalization.

   - **PCA-samples_bedtools.pdf** XXXXXXXXXXXXXXXXXX
   - **diffexpr-maplot_bedtools.pdf**  XXXXXXXXXXXXXXXXXX
   - **pvalue_histogram_50_bedtools.pdf** XXXXXXXXXXXXXXXXXX
   - **qc-dispersions_bedtools.pdf** XXXXXXXXXXXXXXXXXX
   - **qc-heatmap-samples_bedtools.pdf** XXXXXXXXXXXXXXXXXX

   *Note! For more details about differential expression read section XXX in the book and for analysis details check* *[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).*


   __*b. Differential expression for count data generated using *DEXSeq*:*__

   ```bash
   Rscript DE_dexseq_exon.R \
     --exoncountpath analysis/quantification/dexseq-count/ \
     --gffFile analysis/quantification/dexseq-count/DEXSEQ_GTF_annotation.gff \
     --metadata metadata.txt \
     --outputpath analysis/DE/dexseq/
   ```

   *Description:* 
   
   `--exoncountpath` is the dexseq generated exon count file path
   
   `--gffFile` is the dexseq generated gff file with its path
   
   `--metadata` is the sample information file
   
   `--outputpath` is the path to the output directory where the result table and plots will be generated.

   *Output:*
   Differentially expressed output files for bedtools and HTSeq are generated in **analysis/DE/dexseq/**

   ```
   ~/RNASeq_project/analysis/DE/dexseq2$ls -1
   diffexpr-results_dexseq.csv
   dispersion_plot_dexseq.pdf
   ma_plot_dexseq.pdf
   ```

   - **diffexpr-results_dexseq.csv** is the final table from deseq2 for the pairwise comparison between the two conditions with stats.
   - **dispersion_plot_dexseq.pdf** 
   - **ma_plot_dexseq.pdf** XXXXX

   *Note! For more details about differential expression read section XXX in the book and for analysis details check* *[DEXSeq](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html)*


  __*c. Differential expression for count data generated using *cufflinks*:*__

   ```bash
   cuffdiff \
    -o analysis/DE/cuffdiff \
    -L vt,dt \
    --FDR 0.01 \
    -u data/mouse-gencode-version-24/gencode.vM24.annotation.gtf \
    -p 8 \
    BAMLIST
   ```

   *Description:*
   
   `-o` path to the output folder
   
   `-L` lists the labels to be used as “conditions”
   
   `-FDR` cutoff for false discovery rate for the DE analysis
   
   `-u` path to annotation file
   
   `-p` is the number threads
   
   `BAMLIST` is list of all bam file in comma format


**In the end these script will create a folders and sub-folders based on the choice of analysis on the current data as shown below:**

```bash
|-- RNASeq_project
|   | analysis
|   |   |-- mapping
|   |   |   |-- tophat/star
|   |   |-- quantification
|   |   |   |-- bedtools-count/cufflinks-count/dexseq-count/htseq-count
|   |   |-- DE
|   |   |   |-- deseq2/dexseq/cuffdiff
|   | data
|   |   |-- mouse-gencode-version-24/
|   |   |-- STAR_Genome_Index/TOPHAT_Genome_Index/
|   | reads
|   |   |-- rawreads
|   |   |-- rawreads_QA_stats
|   |   |-- filtered_reads
|   | tools
|   |   |-- DEXSeq/FastQC/STAR-2.7.3a/bedtools2/bowtie2/cufflinks/cutadapt/htseq/rseqqc/samtools/sratoolkit/tophat
```
