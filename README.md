#   NGS Pipeline README

This pipeline is designed to process Next-Generation Sequencing (NGS) data, starting from raw FASTQ files and performing steps such as read trimming, alignment, variant calling, and annotation.

##   Table of Contents

* [Dependencies](#dependencies)
* [Installation](#installation)
* [Pipeline Overview](#pipeline-overview)
* [Directory Structure](#directory-structure)
    * [data/](#data)
    * [data/untrimmed_fastq/](#datauntrimmed_fastq)
    * [data/trimmed_fastq/](#data-trimmed_fastq)
    * [data/aligned_data/](#dataaligneddata)
    * [data/reference/](#datareference)
    * [analysis/](#analysis)
    * [results/](#results)
    * [docs/](#docs)
* [Usage](#usage)
* [Pipeline Steps](#pipeline-steps)
* [Input Data](#input-data)
* [Output Data](#output-data)
* [Error Handling](#error-handling)
* [Disk Space Management](#disk-space-management)
* [Cleanup](#cleanup)
* [Notes](#notes)
* [Prerequisites](#prerequisites)
* [Version Control](#version-control)
* [Author](#author)


##   Dependencies

The pipeline relies on the following software:

* **dos2unix:** Used to convert file line endings (if necessary).
    * Installation:

        ```bash
        sudo apt-get install dos2unix
        ```
* **wget:** Used to download data files.
    * Installation:

        ```bash
        sudo apt-get install wget
        ```
* **zcat:** Used to decompress .qz files. Often part of gzip package.
    * Installation: Usually pre-installed, or

        ```bash
        sudo apt-get install gzip
        ```
* **Trimmomatic:** Used to trim low-quality reads and adapter sequences.
    * Installation: Download from [https://github.com/usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic) and follow the instructions. The pipeline assumes Trimmomatic is installed and the adapters file is located at /home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af\_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa. You may need to adjust the path in the script.
* **BWA:** Used to align reads to the reference genome.
    * Installation:

        ```bash
        sudo apt-get install bwa
        ```
* **samtools:** Used to process SAM/BAM files (convert, sort, index, filter).
    * Installation:

        ```bash
        sudo apt-get install samtools
        ```
* **picard:** Used to mark duplicate reads.
    * Installation: Download the latest .jar file from [https://github.com/broadinstitute/picard](https://github.com/broadinstitute/picard). The pipeline assumes picard is in the users PATH, or you need to specify the full path to the jar file in the script.
* **freebayes:** Used for variant calling.
    * Installation:

        ```bash
        sudo apt-get install freebayes
        ```
* **bcftools:** Used as a fallback variant caller.
    * Installation:

        ```bash
        sudo apt-get install bcftools
        ```
* **vcffilter:** Part of the VCFtools package, used to filter variants.
    * Installation:

        ```bash
        sudo apt-get install vcftools
        ```
* **ANNOVAR:** Used to annotate variants.
    * Installation: Download from [https://annovar.openbioinformatics.org/](https://annovar.openbioinformatics.org/) and follow the instructions. You will also need to download the required databases. The pipeline assumes ANNOVAR is installed at \~/annovar.
* **snpEff:** Used to annotate variants.
    * Installation: Download from [http://snpeff.sourceforge.net/](http://snpeff.sourceforge.net/) and follow the instructions. You will also need to download the hg19 database.

##   Installation

1.  **Clone the repository:**

    ```bash
    git clone <repository_url>
    cd <repository_name>
    ```

2.  **Create the necessary directories:** The pipeline script will create these, but you can create them manually if needed.

    ```bash
    mkdir -p ~/ngs_course/assessment_pipeline/data/untrimmed_fastq
    mkdir -p ~/ngs_course/assessment_pipeline/data/trimmed_fastq
    mkdir -p ~/ngs_course/assessment_pipeline/data/aligned_data
    mkdir -p ~/ngs_course/assessment_pipeline/analysis
    mkdir -p ~/ngs_course/assessment_pipeline/results
    mkdir -p ~/ngs_course/assessment_pipeline/docs
    mkdir -p ~/ngs_course/assessment_pipeline/data/reference
    mkdir -p ~/annovar # If you are using the default path
    ```

3.  **Download the reference genome:** The pipeline downloads hg19, but you can place your own reference genome in the data/reference/ directory.
4.  **Download ANNOVAR databases:** The pipeline attempts to download the required databases, but you may need to do this manually depending on your network and ANNOVAR setup.
5.  **Set permissions:** Ensure the script has execute permissions:

    ```bash
    chmod +x pipeline.sh
    ```

##   Pipeline Overview

The pipeline performs the following steps:

1.  Downloads raw FASTQ data.
2.  Decompresses the FASTQ files.
3.  Trims the reads using Trimmomatic.
4.  Prepares the reference genome (downloads if necessary).
5.  Indexes the reference genome using BWA and samtools.
6.  Aligns the trimmed reads to the reference genome using BWA.
7.  Processes the alignment files (converts SAM to BAM, sorts, indexes).
8.  Marks duplicate reads using Picard.
9.  Filters the BAM file.
10. Calls variants using FreeBayes (with bcftools fallback).
11. Filters variants.
12. Annotates variants using ANNOVAR.
13. Annotates variants using snpEff.
14. Performs basic variant prioritization.

##   Directory Structure

ngs_course/assessment_pipeline/├── pipeline.sh├── data/│   ├── untrimmed_fastq/│   ├── trimmed_fastq/│   ├── aligned_data/│   ├── reference/│   └── annotation.bed├── analysis/├── results/└── docs/
###   data/

This directory is the main data directory. It contains the following subdirectories and files:

* **Description:** Stores all input and intermediate data files.
* **Contents:**
    * `untrimmed_fastq/`: Stores the raw, untrimmed FASTQ files.
    * [data/trimmed_fastq/](#data-trimmed_fastq): Stores the trimmed FASTQ files.
    * [data/aligned_data/](#dataaligneddata): Stores the aligned reads in SAM and BAM formats.
    * [data/reference/](#datareference): Stores the reference genome files.
    * `annotation.bed`: Stores the annotation BED file.
* **Relevance:** Essential for organizing the input and output of each processing step.

###   data/untrimmed\_fastq/

* **Description:** Stores the raw, untrimmed FASTQ files.
* **Contents:**
    * `NGS0001.R1.fastq.qz`: Raw reads 1 in compressed format.
    * `NGS0001.R2.fastq.qz`: Raw reads 2 in compressed format.
    * `NGS0001.R1.fastq`: Decompressed reads 1.
    * `NGS0001.R2.fastq`: Decompressed reads 2.
* **Relevance:** Input for the trimming step.
* **Snippet:**

    ```bash
    zcat $FILE_R1 > $UNTRIMMED/R1.fastq
    zcat $FILE_R2 > $UNTRIMMED/R2.fastq
    ```

###   data/trimmed\_fastq/

* **Description:** Stores the trimmed FASTQ files after Trimmomatic processing.
* **Contents:**
    * `trimmed_1P.fastq.gz`: Trimmed paired-end reads 1 (compressed).
    * `trimmed_1U.fastq.gz`: Trimmed unpaired reads 1 (compressed).
    * `trimmed_2P.fastq.gz`: Trimmed paired-end reads 2 (compressed).
    * `trimmed_2U.fastq.gz`: Trimmed unpaired reads 2 (compressed).
* **Relevance:** Output from the trimming step and input for the alignment step.
* **Snippet:**

    ```bash
    trimmomatic PE -threads 4 -phred33 $R1 $R2 \
        $TRIMMED/trimmed_1P.fastq.gz $TRIMMED/trimmed_1U.fastq.gz \
        $TRIMMED/trimmed_2P.fastq.gz $TRIMMED/trimmed_2U.fastq.gz ...
    ```

###   data/aligned\_data/

* **Description:** Stores the aligned reads in SAM and BAM formats.
* **Contents:**
    * `NGS01.sam`: Aligned reads in SAM format (intermediate).
    * `NGS01.bam`: Aligned reads in BAM format (intermediate).
    * `NGS01_sorted.bam`: Sorted BAM file.
    * `NGS01_marked.bam`: BAM file with marked duplicates.
    * `NGS01_filtered.bam`: Filtered BAM file.
* **Relevance:** Output from the alignment and duplicate marking steps; input for variant calling.
* **Snippet:**

    ```bash
    bwa mem "$REF_FA" $TRIMMED/trimmed_1P.fastq.gz  $TRIMMED/trimmed_2P.fastq.gz -o $SAM
    samtools view -b $SAM > $BAM
    samtools sort $BAM > $SORT_BAM
    picard MarkDuplicates I=$SORT_BAM O=$MARK_BAM M=$ALIGNED/marked_dup_metrics.txt
    samtools view -F 0x4 -q 20 -b $MARK_BAM > $FILTER_BAM
    ```

###   data/reference/

* **Description:** Stores the reference genome files.
* **Contents:**
    * `hg19.fa.gz`: Compressed reference genome (downloaded).
    * `hg19.fa`: Decompressed reference genome.
    * `hg19.fa.fai`: Samtools index file.
    * `hg19.1.bt2` etc: BWA index files.
* **Relevance:** Required for alignment and variant calling.
* **Snippet:**

    ```bash
    bwa index "$REF_FA" # BWA Indexing
    samtools faidx $REF_FA #Samtools Indexing
    ```

###   analysis/

* **Description:** Stores the primary variant calling and annotation results.
* **Contents:**
    * `NGS01.vcf`: Raw variant calls.
    * `NGS01_filtered.vcf`: Filtered variant calls.
    * `NGS01_annovar.*`: ANNOVAR annotation files.
    * `NGS01_prioritized_variants.txt`: Prioritized variants.
* **Relevance:** Contains the main processed data for analysis.
* **Snippet: Variant Calling**

    ```bash
    freebayes -f "$REF_FA" -v "$VCF" "$FILTER_BAM"
    ```
* **Snippet: ANNOVAR Annotation**

    ```bash
    $ANNOVAR_DIR/table_annovar.pl $FILTER_VCF $HUMANDB_DIR ...
    ```

###   results/

* **Description:** Stores additional results, specifically snpEff annotations.
* **Contents:**
    * `NGS01_snpEff.vcf`: snpEff annotation file.
* **Relevance:** Contains results from snpEff variant annotation.
* **Snippet: snpEff Annotation**

    ```bash
    java -Xmx4g -jar "$snpeff_dir/snpEff.jar" "$genome_version" "$filter_vcf" > "$results_dir/NGS01_snpEff.vcf"
    ```

###   docs/

* **Description:** This directory stores log files and other documentation produced by the pipeline.  It can also be used to track download progress.
* **Contents:**
    * `/home/ubuntu/bioinformatics_course/scripts/docs/available_dbs.txt`:  A file produced by the pipeline.
    * Other log files.
* **Relevance:** Important for tracking pipeline execution, debugging, and recording download status.

##   Usage

To use the pipeline, follow these steps:

1.  **Clone the repository to your local machine:**

    ```bash
    git clone https://github.com/LDolanLDolan/bioinformatics_course
    ```

2.  **Navigate to the pipeline directory:**

    ```bash
    cd bioinformatics_course
    cd scripts
    ```

3.  **Ensure the pipeline script has execute permissions:**

    ```bash
    chmod +x pipeline_one.sh
    ```

4.  **Run the pipeline script:**

    ```bash
    ./pipeline_one.sh
    ```

##   Pipeline Steps

The pipeline consists of the following steps:

1.  [Downloading raw data files](#downloading-raw-data-files): Downloads the raw FASTQ files and the annotation BED file using wget.
2.  [Decompressing FASTQ files](#decompressing-fastq-files): Decompresses the downloaded .qz files using zcat.
3.  [Trimming reads](#trimming-reads): Trims the reads using Trimmomatic to remove low-quality bases and adapter sequences.
4.  [Preparing reference genome](#preparing-reference-genome): Downloads the hg19 reference genome if it doesn't exist and decompresses it.
5.  [Indexing reference genome](#indexing-reference-genome): Indexes the reference genome using BWA and samtools faidx.
6.  [Aligning reads](#aligning-reads): Aligns the trimmed reads to the reference genome using BWA.
7.  [Processing alignment files](#processing-alignment-files): Converts the SAM file to BAM format, sorts the BAM file by coordinate, and indexes the sorted BAM file using samtools.
8.  [Marking duplicates](#marking-duplicates): Marks duplicate reads using Picard MarkDuplicates to prepare for variant calling.
9.  [Filtering BAM file](#filtering-bam-file): Filters the BAM file to keep only properly paired reads with a mapping quality greater than 20.
10. [Variant calling](#variant-calling): Calls variants using FreeBayes, with a fallback to bcftools if FreeBayes fails.
11. [Filtering variants](#filtering-variants): Filters the variants using vcffilter.
12. [Variant annotation (ANNOVAR)](#variant-annotation-annovar): Annotates the variants using ANNOVAR.
13.  [Variant annotation (snpEff)](#variant-annotation-snpeff): Annotates the variants using snpEff.
14. [Variant prioritization](#variant-prioritization): Performs a basic variant prioritization step.

##   Input Data

* Raw FASTQ files (NGS0001.R1.fastq.qz, NGS0001.R2.fastq.qz)
* Annotation BED file (annotation.bed)
* Reference genome (hg19.fa)

##   Output Data

The pipeline generates the following output data:

* Trimmed FASTQ files
* Aligned BAM files
* Variant calls in VCF format
* Annotated variants (ANNOVAR and snpEff)
* Prioritized variants

##   Error Handling

The pipeline includes error handling at each step. It checks for:

* Successful download of files.
* Successful execution of tools.
* Existence and size of output files.
* Sufficient disk space.

If an error occurs, the pipeline will print an error message and exit.

##   Disk Space Management

The pipeline includes checks for available disk space before performing disk-intensive operations. It also attempts to clean up intermediate files to free up space.

##   Cleanup

The pipeline attempts to remove unnecessary intermediate files, such as:

* Uncompressed raw FASTQ files.
* SAM files after BAM conversion.
* Unsorted BAM files after sorting.
* Compressed reference genome after decompression.
* Temporary files in /tmp.

##   Notes

* The pipeline assumes that all required software is installed and configured correctly. You may need to adjust paths to executables and adapter files.
* The pipeline downloads the hg19 reference genome. You can replace this with your own reference genome if needed.
* The ANNOVAR databases are downloaded dynamically by the script. Ensure you have a working internet connection and that ANNOVAR is configured correctly.
* The variant prioritization step is basic. You may need to modify it for your specific needs.
* Consider adding a logging mechanism to record the pipeline's progress and any errors.
* This script is designed to be run on a Linux system. You may need to modify it for other operating systems.

##   Prerequisites

* **Bash**: Ensure you have Bash installed (version 4.0 or higher).
* **BWA**: For sequence alignment.
* **SAMtools**: For processing SAM/BAM files.
* **Picard**: For marking duplicates.
* **FreeBayes**: For variant calling.
* **ANNOVAR**: For variant annotation.

##   Version Control

It is recommended to use version control (e.g., Git) to track changes to the pipeline script.

##   Author

[Lita Doolan]
[litadoolan.net]

