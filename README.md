# **NGS Pipeline README**

This pipeline is designed to process Next-Generation Sequencing (NGS) data, starting from raw FASTQ files and performing steps such as read trimming, alignment, variant calling, and annotation.

## **Table of Contents**

* [Dependencies](#bookmark=id.bnsyvnndaehz)  
* [Installation](#bookmark=id.2a4ohkwmb0y7)  
* [Pipeline Overview](#bookmark=id.fb9jxiqw92qy)  
* [Directory Structure](#bookmark=id.aq0dhkwmve5r)  
  * [data/](#bookmark=id.2ia070g06hig)  
  * [data/untrimmed\_fastq/](#bookmark=id.ie2lyhyci7en)  
  * [data/trimmed\_fastq/](#bookmark=id.ewdx444luy45)  
  * [data/aligned\_data/](#bookmark=id.d0mgjz1jgxpr)  
  * [data/reference/](#bookmark=id.jxop93sla9pm)  
  * [results/](#bookmark=id.3yknbjqkgpzw)  
  * [logs/](#bookmark=id.g60hpn8bidgg)  
* [Usage](#bookmark=id.9oygho10z0df)  
* [Pipeline Steps](#bookmark=id.dgyyme2otlcs)  
* [Input Data](#bookmark=id.aeyzq5wfb47r)  
* [Output Data](#bookmark=id.6d8a7h6qjak6)  
* [Error Handling](#bookmark=id.ofgexj45vy98)  
* [Disk Space Management](#bookmark=id.itfunzawdqlb)  
* [Cleanup](#bookmark=id.ucqnij4xevqx)  
* [Notes](#bookmark=id.77677525dg5f)  
* [Version Control](#bookmark=id.tksk3gb8z2bn)  
* [Author](#bookmark=id.yh0nwjtj6op)  
* [License](#bookmark=id.t2q0lopaaeuv)

## **Dependencies**

The pipeline relies on the following software:

* **dos2unix:** Used to convert file line endings (if necessary).  
  * Installation: sudo apt-get install dos2unix  
* **wget:** Used to download data files.  
  * Installation: sudo apt-get install wget  
* **zcat:** Used to decompress .qz files. Often part of gzip package.  
  * Installation: Usually pre-installed, or sudo apt-get install gzip  
* **Trimmomatic:** Used to trim low-quality reads and adapter sequences.  
  * Installation: Download from [https://github.com/usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic) and follow the instructions. The pipeline assumes Trimmomatic is installed and the adapters file is located at /home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af\_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa. You may need to adjust the path in the script.  
* **Bowtie 2:** Used to align reads to the reference genome.  
  * Installation: sudo apt-get install bowtie2  
* **samtools:** Used to process SAM/BAM files (convert, sort, index, filter).  
  * Installation: sudo apt-get install samtools  
* **picard:** Used to mark duplicate reads.  
  * Installation: Download the latest .jar file from [https://github.com/broadinstitute/picard](https://github.com/broadinstitute/picard). The pipeline assumes picard is in the users PATH, or you need to specify the full path to the jar file in the script.  
* **freebayes:** Used for variant calling.  
  * Installation: sudo apt-get install freebayes  
* **bcftools:** Used as a fallback variant caller.  
  * Installation: sudo apt-get install bcftools  
* **vcffilter:** Part of the VCFtools package, used to filter variants.  
  * Installation: sudo apt-get install vcftools  
* **ANNOVAR:** Used to annotate variants.  
  * Installation: Download from [https://annovar.openbioinformatics.org/](https://annovar.openbioinformatics.org/) and follow the instructions. You will also need to download the required databases. The pipeline assumes ANNOVAR is installed at \~/annovar.  
* **snpEff:** Used to annotate variants.  
  * Installation: Download from [http://snpeff.sourceforge.net/](http://snpeff.sourceforge.net/) and follow the instructions. You will also need to download the hg19 database.

## **Installation**

1. **Clone the repository:**  
   git clone \<repository\_url\>  
   cd \<repository\_name\>

2. **Create the necessary directories:** The pipeline script will create these, but you can create them manually if needed.  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/data/untrimmed\_fastq  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/data/trimmed\_fastq  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/data/aligned\_data  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/results  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/logs  
   mkdir \-p \~/ngs\_course/assessment\_pipeline/data/reference  
   mkdir \-p \~/annovar  \# If you are using the default path

3. **Download the reference genome:** The pipeline downloads hg19, but you can place your own reference genome in the data/reference/ directory.  
4. **Download ANNOVAR databases:** The pipeline attempts to download the required databases, but you may need to do this manually depending on your network and ANNOVAR setup.  
5. **Set permissions:** Ensure the script has execute permissions:  
   chmod \+x pipeline.sh

## **Pipeline Overview**

The pipeline performs the following steps:

1. Downloads raw FASTQ data.  
2. Decompresses the FASTQ files.  
3. Trims the reads using Trimmomatic.  
4. Prepares the reference genome (downloads if necessary).  
5. Indexes the reference genome using Bowtie2 and samtools.  
6. Aligns the trimmed reads to the reference genome using Bowtie2.  
7. Processes the alignment files (converts SAM to BAM, sorts, indexes).  
8. Marks duplicate reads using Picard.  
9. Filters the BAM file.  
10. Calls variants using FreeBayes (with bcftools fallback).  
11. Filters variants.  
12. Annotates variants using ANNOVAR.  
13. Annotates variants using SnpEff.  
14. Performs basic variant prioritization.

## **Directory Structure**

ngs\_course/assessment\_pipeline/  
├── pipeline.sh  
├── data/  
│   ├── untrimmed\_fastq/  
│   ├── trimmed\_fastq/  
│   ├── aligned\_data/  
│   ├── reference/  
│   └── annotation.bed  
├── results/  
└── logs/

### **data/**

This directory is the main data directory. It contains the following subdirectories and files:

* **Description:** Stores all input and intermediate data files.  
* **Contents:**  
  * untrimmed\_fastq/: Stores the raw, untrimmed FASTQ files.  
  * trimmed\_fastq/: Stores the trimmed FASTQ files.  
  * aligned\_data/: Stores the aligned reads in SAM and BAM formats.  
  * reference/: Stores the reference genome files.  
  * annotation.bed: Stores the annotation BED file.  
* **Relevance:** Essential for organizing the input and output of each processing step.

### **data/untrimmed\_fastq/**

* **Description:** Stores the raw, untrimmed FASTQ files.  
* **Contents:**  
  * NGS0001.R1.fastq.qz: Raw reads 1 in compressed format.  
  * NGS0001.R2.fastq.qz: Raw reads 2 in compressed format.  
  * NGS0001.R1.fastq: Decompressed reads 1\.  
  * NGS0001.R2.fastq: Decompressed reads 2\.  
* **Relevance:** Input for the trimming step.  
* **Snippet:**  
  zcat $FILE\_R1 \> $UNTRIMMED/R1.fastq  
  zcat $FILE\_R2 \> $UNTRIMMED/R2.fastq

### **data/trimmed\_fastq/**

* **Description:** Stores the trimmed FASTQ files after Trimmomatic processing.  
* **Contents:**  
  * trimmed\_1P.fastq.gz: Trimmed paired-end reads 1 (compressed).  
  * trimmed\_1U.fastq.gz: Trimmed unpaired reads 1 (compressed).  
  * trimmed\_2P.fastq.gz: Trimmed paired-end reads 2 (compressed).  
  * trimmed\_2U.fastq.gz: Trimmed unpaired reads 2 (compressed).  
* **Relevance:** Output from the trimming step and input for the alignment step.  
* **Snippet:**  
  trimmomatic PE \-threads 4 \-phred33 $R1 $R2 \\  
    $TRIMMED/trimmed\_1P.fastq.gz $TRIMMED/trimmed\_1U.fastq.gz \\  
    $TRIMMED/trimmed\_2P.fastq.gz $TRIMMED/trimmed\_2U.fastq.gz ...

### **data/aligned\_data/**

* **Description:** Stores the aligned reads in SAM and BAM formats.  
* **Contents:**  
  * NGS01.sam: Aligned reads in SAM format (intermediate).  
  * NGS01.bam: Aligned reads in BAM format (intermediate).  
  * NGS01\_sorted.bam: Sorted BAM file.  
  * NGS01\_marked.bam: BAM file with marked duplicates.  
  * NGS01\_filtered.bam: Filtered BAM file.  
* **Relevance:** Output from the alignment and duplicate marking steps; input for variant calling.  
* **Snippet:**  
  bowtie2 \-p 4 \-x "$REF\_FA" \-1 $TRIMMED/trimmed\_1P.fastq.gz \-2 $TRIMMED/trimmed\_2P.fastq.gz \-S $SAM  
  samtools view \-b $SAM \> $BAM  
  samtools sort $BAM \> $SORT\_BAM  
  picard MarkDuplicates I=$SORT\_BAM O=$MARK\_BAM M=$ALIGNED/marked\_dup\_metrics.txt  
  samtools view \-F 0x4 \-q 20 \-b $MARK\_BAM \> $FILTER\_BAM

### **data/reference/**

* **Description:** Stores the reference genome files.  
* **Contents:**  
  * hg19.fa.gz: Compressed reference genome (downloaded).  
  * hg19.fa: Decompressed reference genome.  
  * hg19.fa.fai: Samtools index file.  
  * hg19.1.bt2 etc: Bowtie2 index files.  
* **Relevance:** Required for alignment and variant calling.  
* **Snippet:**  
  bowtie2-build "$REF\_FA" "$REF\_FA" \#Bowtie2 Indexing  
  samtools faidx $REF\_FA \#Samtools Indexing

### **results/**

* **Description:** Stores the final results of the pipeline.  
* **Contents:**  
  * NGS01.vcf: Raw variant calls.  
  * NGS01\_filtered.vcf: Filtered variant calls.  
  * NGS01\_annovar.\*: ANNOVAR annotation files.  
  * NGS01\_snpEff.vcf: snpEff annotation file.  
  * NGS01\_prioritized\_variants.txt: Prioritized variants.  
* **Relevance:** Contains the final, processed data for analysis.  
* **Snippet: Variant Calling**  
  freebayes \-f "$REF\_FA" \-v "$VCF" "$FILTER\_BAM"

* **Snippet: ANNOVAR Annotation**  
  $ANNOVAR\_DIR/table\_annovar.pl $FILTER\_VCF $HUMANDB\_DIR ...

### **logs/**

* **Description**: Although the script does not explicitly create or write to a log file, this directory is present, and the user may add logging.  
* **Contents**: Ideally, this directory should contain log files of the pipeline's progress, including any errors or warnings.  
* **Relevance**: Important for debugging and monitoring the pipeline's execution.  
* **Potential Snippet**:  
  \#  Example of redirecting standard error to a log file (add to relevant commands)  
  freebayes \-f "$REF\_FA" \-v "$VCF" "$FILTER\_BAM" 2\> $LOGS/freebayes.log

## **Usage**

1. **Navigate to the pipeline directory:**  
   cd \~/ngs\_course/assessment\_pipeline

2. **Run the pipeline script:**  
   ./pipeline.sh

## **Pipeline Steps**

The pipeline consists of the following steps:

1. **Downloading raw data files:** Downloads the raw FASTQ files and the annotation BED file using wget.  
2. **Decompressing FASTQ files:** Decompresses the downloaded .qz files using zcat.  
3. **Trimming reads:** Trims the reads using Trimmomatic to remove low-quality bases and adapter sequences.  
4. **Preparing reference genome:** Downloads the hg19 reference genome if it doesn't exist and decompresses it.  
5. **Indexing reference genome:** Indexes the reference genome using bowtie2-build and samtools faidx.  
6. **Aligning reads:** Aligns the trimmed reads to the reference genome using Bowtie2.  
7. **Processing alignment files:** Converts the SAM file to BAM format, sorts the BAM file by coordinate, and indexes the sorted BAM file using samtools.  
8. **Marking duplicates:** Marks duplicate reads using Picard MarkDuplicates to prepare for variant calling.  
9. **Filtering BAM file:** Filters the BAM file to keep only properly paired reads with a mapping quality greater than 20\.  
10. **Variant calling:** Calls variants using FreeBayes, with a fallback to bcftools if FreeBayes fails.  
11. **Filtering variants:** Filters the variants using vcffilter.  
12. **Variant annotation (ANNOVAR):** Annotates the variants using ANNOVAR.  
13. **Variant annotation (snpEff):** Annotates the variants using snpEff.  
14. **Variant prioritization:** Performs a basic variant prioritization step.

## **Input Data**

* Raw FASTQ files (NGS0001.R1.fastq.qz, NGS0001.R2.fastq.qz)  
* Annotation BED file (annotation.bed)  
* Reference genome (hg19.fa)

## **Output Data**

The pipeline generates the following output data:

* Trimmed FASTQ files  
* Aligned BAM files  
* Variant calls in VCF format  
* Annotated variants (ANNOVAR and snpEff)  
* Prioritized variants

## **Error Handling**

The pipeline includes error handling at each step. It checks for:

* Successful download of files.  
* Successful execution of tools.  
* Existence and size of output files.  
* Sufficient disk space.

If an error occurs, the pipeline will print an error message and exit.

## **Disk Space Management**

The pipeline includes checks for available disk space before performing disk-intensive operations. It also attempts to clean up intermediate files to free up space.

## **Cleanup**

The pipeline attempts to remove unnecessary intermediate files, such as:

* Uncompressed raw FASTQ files.  
* SAM files after BAM conversion.  
* Unsorted BAM files after sorting.  
* Compressed reference genome after decompression.  
* Temporary files in /tmp.

## **Notes**

* The pipeline assumes that all required software is installed and configured correctly. You may need to adjust paths to executables and adapter files.  
* The pipeline downloads the hg19 reference genome. You can replace this with your own reference genome if needed.  
* The ANNOVAR databases are downloaded dynamically by the script. Ensure you have a working internet connection and that ANNOVAR is configured correctly.  
* The variant prioritization step is basic. You may need to modify it for your specific needs.  
* Consider adding a logging mechanism to record the pipeline's progress and any errors.  
* This script is designed to be run on a Linux system. You may need to modify it for other operating systems.

## **Version Control**

It is recommended to use version control (e.g., Git) to track changes to the pipeline script.

## 
