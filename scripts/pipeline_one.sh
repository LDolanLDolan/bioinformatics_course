 #!/bin/bash

# Check and install dos2unix if not present
if ! command -v dos2unix &> /dev/null; then
    echo "dos2unix not found. Installing..."
    sudo apt-get update
    sudo apt-get install -y dos2unix
else
    echo "dos2unix already installed."
fi

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed."
    echo "Please install Conda before running this pipeline."
    echo "You can download Miniconda or Anaconda from:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Array of tools to check
tools=("freebayes" "samtools" "bwa" "bowtie2" "fastqc" "trimmomatic" "bedtools" "vcftools" "bcftools" "picard")

echo -e "${YELLOW}Checking for NGS tools installation...${NC}"
echo

missing_tools=()
for tool in "${tools[@]}"; do
    echo -n "Checking for $tool: "
    if command -v "$tool" &> /dev/null; then
        version=$($tool -version 2>&1 | head -n 1)
        echo -e "${GREEN}INSTALLED${NC} - $version"
    else
        echo -e "${RED}NOT INSTALLED${NC}"
        missing_tools+=("$tool")
    fi
done

# Check for ANNOVAR
echo -n "Checking for ANNOVAR: "
if command -v table_annovar.pl &> /dev/null; then
    echo -e "${GREEN}INSTALLED${NC}"
else
    echo -e "${RED}NOT INSTALLED${NC}"
    missing_tools+=("ANNOVAR")
fi

# Check Java for snpEff
echo
echo -n "Checking for Java (required for snpEff): "
if command -v java &> /dev/null; then
    version=$(java -version 2>&1 | head -n 1)
    echo -e "${GREEN}INSTALLED${NC} - $version"
else
    echo -e "${RED}NOT INSTALLED${NC}"
    echo "You need to install Java before installing snpEff:"
    echo "sudo apt install default-jre"
    missing_tools+=("Java")  # Add Java to the missing_tools array if it's not installed
fi

# If any tools are missing, exit with a message
if [ ${#missing_tools[@]} -gt 0 ]; then
    echo -e "${RED}The following tools are missing: ${missing_tools[@]}${NC}"
    exit 1
fi

echo -e "${GREEN}All required tools are installed!${NC}"


# Array of old file paths to delete that may be present from previous runs
echo "Deleting old files..."
files_to_delete=(
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_sorted.bam.bai"
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_sorted.bam"
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_marked.bam.bai"
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_marked.bam"
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_filtered.bam"
    "/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_filtered.bam.bai"
    "/home/ubuntu/bioinformatics_course/scripts/analyis/NGS01_annovar.avinput"
    "/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01.vcf"
    "/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01_snpEff.vcf"
    "/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01_filtered.vcf"
    "/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01_annovar.hg19_multianno.vcf"
)

# Loop through the array and delete each file if it exists
for file in "${files_to_delete[@]}"; do
    if [ -f "$file" ]; then
        rm "$file"
        echo "Deleted: $file"
    fi
done



# Directories
BASE=~/bioinformatics_course/scripts
DATA=$BASE/data
UNTRIMMED=$DATA/untrimmed_fastq
REF=$DATA/reference
TRIMMED=$DATA/trimmed_fastq
ALIGNED=$DATA/aligned_data
RESULTS=$BASE/analysis
LOGS=$BASE/docs
ANNOVAR_DIR=~/annovar
FILE_ANNOTATION=$DATA/annotation.bed  # Full path for annotation.bed file

# Create directories (with -p to avoid errors if they already exist)
mkdir -p $UNTRIMMED $TRIMMED $ALIGNED $RESULTS $LOGS $DATA $ANNOVAR_DIR $REF

# Output file names
SAM=$ALIGNED/NGS01.sam
BAM=$ALIGNED/NGS01.bam
SORT_BAM=$ALIGNED/NGS01_sorted.bam
MARK_BAM=$ALIGNED/NGS01_marked.bam
FILTER_BAM=$ALIGNED/NGS01_filtered.bam
VCF=$RESULTS/NGS01.vcf
FILTER_VCF=$RESULTS/NGS01_filtered.vcf



# Reference genome
REF_GZ=$REF/hg19.fa.gz
REF_FA=$REF/hg19.fa

# URLs of the input files
URL_R1="https://emea01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fs3-eu-west-1.amazonaws.com%2Fworkshopdata2017%2FNGS0001.R1.fastq.qz&data=05%7C02%7C%7C840fe0668565488531dc08dd6a499765%7C84df9e7fe9f640afb435aaaaaaaaaaaa%7C1%7C0%7C638783586482395346%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=%2F13acujjahhi3aUJRgtOm84LQoxYBzl5kIRWVgKewhY%3D&reserved=0"
URL_R2="https://emea01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fs3-eu-west-1.amazonaws.com%2Fworkshopdata2017%2FNGS0001.R2.fastq.qz&data=05%7C02%7C%7C840fe0668565488531dc08dd6a499765%7C84df9e7fe9f640afb435aaaaaaaaaaaa%7C1%7C0%7C638783586482416855%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=Ga%2Fu9HmiW13mLq9vcJB8rzsu24GlUUNFLmcEabKz7%2FQ%3D&reserved=0"
URL_ANNOTATION="https://emea01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fs3-eu-west-1.amazonaws.com%2Fworkshopdata2017%2Fannotation.bed&data=05%7C02%7C%7C840fe0668565488531dc08dd6a499765%7C84df9e7fe9f640afb435aaaaaaaaaaaa%7C1%7C0%7C638783586482426312%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=1bUbknTG3Ahgnh2TYqzIq1y54F%2BMnQj8j8BQTMVrNDs%3D&reserved=0"

# Output file names for raw data
FILE_R1=$UNTRIMMED/NGS0001.R1.fastq.qz
FILE_R2=$UNTRIMMED/NGS0001.R2.fastq.qz

FILE_ANNOTATION=$DATA/annotation.bed  # Full path for annotation.bed file

# Check if the annotation file exists, if not, download it

if [ ! -f "$FILE_ANNOTATION" ]; then
    echo "Annotation file not found. Downloading..."
    wget -q $URL_ANNOTATION -O $FILE_ANNOTATION

    # Check if the file was successfully downloaded
    if [ ! -f "$FILE_ANNOTATION" ]; then
        echo "Error: Annotation file could not be downloaded."
        exit 1
    fi
fi

if [ ! -r "$FILE_ANNOTATION" ]; then
    echo "Error: Annotation file $FILE_ANNOTATION is not readable."
    exit 1
fi



# Function to check if a bwa index is empty
bwa_index_valid() {
    local ref="$1"
    local required_exts=(".amb" ".ann" ".bwt" ".pac" ".sa")
    for ext in "${required_exts[@]}"; do
        if [[ ! -s "${ref}${ext}" ]]; then
            return 1  # At least one index file is missing or empty
        fi
    done
    return 0  # All required files exist and are non-empty
}

# Function to check if a file is empty
check_file_empty() {
    if [ ! -s "$1" ]; then
        echo "Error: $1 is empty or failed to download."
        return 1
    else
        echo "$1 downloaded successfully."
        return 0
    fi
}

# Function to check disk space and cleanup if needed
check_disk_space() {
    local required_mb=$1
    local required_kb=$((required_mb * 1024))
    local available_kb=$(df -k . | awk 'NR==2 {print $4}')

    if [ $available_kb -lt $required_kb ]; then
        echo "WARNING: Low disk space detected. Available: $((available_kb / 1024))MB, Required: ${required_mb}MB"
        echo "Attempting to free up space..."
        return 1
    else
        echo "Sufficient disk space available: $((available_kb / 1024))MB"
        return 0
    fi
}

# Function to clean up unnecessary files to free space
cleanup_files() {
    echo "Cleaning up files to free space..."

    # Remove core dump files
    find $BASE -name "core.*" -delete

    # Remove temporary and unnecessary files
    # Only delete if processing has moved past that stage
    if [ -f "$TRIMMED/trimmed_1P.fastq.gz" ] && [ -f "$TRIMMED/trimmed_2P.fastq.gz" ]; then
        echo "Removing uncompressed raw fastq files..."
        rm -f $UNTRIMMED/R1.fastq $UNTRIMMED/R2.fastq
    fi

    # If raw data has already been processed, remove original downloaded files
    if [ -f "$SAM" ] && [ -s "$SAM" ]; then
        echo "Removing original .qz files..."
        rm -f $FILE_R1 $FILE_R2
    fi

    # If BAM file exists, remove SAM file
    if [ -f "$BAM" ] && [ -s "$BAM" ]; then
        echo "Removing SAM file as BAM exists..."
        rm -f $SAM
    fi

    # Clean up any temporary files in /tmp
    if [ -d "/tmp" ]; then
        echo "Cleaning up temporary files in /tmp..."
        find /tmp -user $(whoami) -type f -mtime +1 -delete 2>/dev/null || true
    fi

    echo "Cleanup completed."
}

echo "=== STEP 1: Downloading raw data files ==="
# Check disk space before downloading
if ! check_disk_space 1000; then
    cleanup_files
    if ! check_disk_space 1000; then
        echo "ERROR: Not enough disk space even after cleanup. Cannot proceed."
        exit 1
    fi
fi

# Download the raw data files using wget
wget -q $URL_R1 -O $FILE_R1
if ! check_file_empty $FILE_R1; then
    echo "Retrying download from alternate source..."
    # an alternate source URL can be added here if available
    # If no alternate source, exit
    exit 1
fi

wget -q $URL_R2 -O $FILE_R2
if ! check_file_empty $FILE_R2; then
    echo "Retrying download from alternate source..."
    # an alternate source URL here if available
    # If no alternate source, exit
    exit 1
fi

wget -q $URL_ANNOTATION -O $FILE_ANNOTATION
if ! check_file_empty $FILE_ANNOTATION; then
    echo "Retrying download from alternate source..."
    # You can add an alternate source URL here if available
    # If no alternate source, exit
    exit 1
fi

echo "=== STEP 2: Decompressing fastq files ==="
# Check disk space before decompression
if ! check_disk_space 5000; then
    cleanup_files
    if ! check_disk_space 5000; then
        echo "ERROR: Not enough disk space for decompression even after cleanup."
        exit 1
    fi
fi

# Decompress .qz files
echo "Decompressing .qz files..."
zcat $FILE_R1 > $UNTRIMMED/R1.fastq
zcat $FILE_R2 > $UNTRIMMED/R2.fastq

# Define input variables for Trimmomatic using the decompressed files
R1=$UNTRIMMED/R1.fastq
R2=$UNTRIMMED/R2.fastq

# Check if decompression worked
if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
    echo "Error: Decompression failed. Check if the downloaded files are in the correct format."
    exit 1
fi

echo "=== STEP 3: Trimming reads ==="
# Check disk space before trimming
if ! check_disk_space 10000; then
    cleanup_files
    if ! check_disk_space 10000; then
        echo "ERROR: Not enough disk space for trimming even after cleanup."
        exit 1
    fi
fi

# Run Trimmomatic
echo "Trimming reads using Trimmomatic..."
trimmomatic PE -threads 4 -phred33 $R1 $R2 \
    $TRIMMED/trimmed_1P.fastq.gz $TRIMMED/trimmed_1U.fastq.gz \
    $TRIMMED/trimmed_2P.fastq.gz $TRIMMED/trimmed_2U.fastq.gz \
    ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
    TRAILING:25 MINLEN:50

# Check if Trimmomatic ran successfully
if [ $? -ne 0 ]; then
    echo "Error: Trimmomatic failed."
    exit 1
fi

# Count and report how many reads were retained
echo "Counting reads in trimmed output..."
FORWARD_READS=$(zcat $TRIMMED/trimmed_1P.fastq.gz | wc -l)
REVERSE_READS=$(zcat $TRIMMED/trimmed_2P.fastq.gz | wc -l)
echo "Trimmed paired forward reads: $FORWARD_READS lines ($(($FORWARD_READS/4)) reads)"
echo "Trimmed paired reverse reads: $REVERSE_READS lines ($(($REVERSE_READS/4)) reads)"

# Check if trimmed files are valid
if [ $FORWARD_READS -eq 0 ] || [ $REVERSE_READS -eq 0 ]; then
    echo "Error: Trimming produced empty output files."
    exit 1
fi

# After successful trimming, we can remove the raw uncompressed files
echo "Removing uncompressed raw fastq files to save space..."
rm -f $UNTRIMMED/R1.fastq $UNTRIMMED/R2.fastq

echo "=== STEP 4: Preparing reference genome ==="
# Check disk space before downloading reference genome
if ! check_disk_space 10000; then
    cleanup_files
    if ! check_disk_space 10000; then
        echo "ERROR: Not enough disk space for reference genome even after cleanup."
        exit 1
    fi
fi

# Reference genome URL (UCSC hg19)
HG19_URL="https://emea01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fhgdownload.soe.ucsc.edu%2FgoldenPath%2Fhg19%2FbigZips%2Fhg19.fa.gz&data=05%7C02%7C%7C840fe0668565488531dc08dd6a499765%7C84df9e7fe9f640afb435aaaaaaaaaaaa%7C1%7C0%7C638783586482435611%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=YjX%2Fg%2BaWuK7h85qWfxEzDi04e2ohQEHg3qxFUDiyZBg%3D&reserved=0"
# Alternative URL in case the primary one fails
ALT_HG19_URL="https://emea01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fstorage.googleapis.com%2Fgenomics-public-data%2Fresources%2Fbroad%2Fhg19%2Fv0%2FHomo_sapiens_assembly19.fasta&data=05%7C02%7C%7C840fe0668565488531dc08dd6a499765%7C84df9e7fe9f640afb435aaaaaaaaaaaa%7C1%7C0%7C638783586482444585%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=T30%2F4QFHWVMwprktFdcPkyKcOGydJLjuxHX8C5SHUvU%3D&reserved=0"

# Check if the reference genome exists and is not empty
if [ ! -s "$REF_FA" ]; then
    echo "Reference genome hg19.fa not found or is empty. Downloading..."
    wget -q $HG19_URL -O $REF_GZ

    # Check if download was successful
    if [ ! -s "$REF_GZ" ]; then
        echo "Error: Failed to download hg19.fa.gz or file is empty."
        echo "Trying alternative source..."
        wget -q $ALT_HG19_URL -O $REF_FA

        if [ ! -s "$REF_FA" ]; then
            echo "Error: Failed to download reference genome from alternative source."
            exit 1
        fi
    else
        # Decompress the genome
        echo "Decompressing hg19.fa.gz..."
        gunzip -c $REF_GZ > $REF_FA

        # Check if the decompressed file exists and is not empty
        if [ ! -s "$REF_FA" ]; then
            echo "Error: Failed to decompress hg19.fa.gz or file is empty."
            exit 1
        fi

        # After successful decompression, remove the compressed file to save space
        echo "Removing compressed reference file to save space..."
        rm -f $REF_GZ
    fi
else
    echo "Reference genome hg19.fa exists and is not empty."
fi

echo "=== STEP 5: Indexing reference genome ==="

# Function to check if all BWA index files exist and are valid
bwa_index_valid() {
    local REF=$1
    local INDEX_FILES=(".amb" ".ann" ".bwt" ".pac" ".sa")

    for ext in "${INDEX_FILES[@]}"; do
        if [ ! -s "${REF}${ext}" ]; then
            return 1  # At least one index file is missing or empty
        fi
    done

    return 0  # All index files exist and are not empty
}


# Only run BWA index if the index files are missing or invalid
if bwa_index_valid "$REF_FA"; then
    echo "BWA index files already exist and are valid. Skipping indexing."
else
    echo "BWA index files missing or incomplete. Running bwa index..."
    bwa index $REF_FA
    if [ $? -ne 0 ]; then
        echo "Error: bwa index failed."
        exit 1
    fi
fi



# Check disk space before indexing
if ! check_disk_space 15000; then
    cleanup_files
    if ! check_disk_space 15000; then
        echo "ERROR: Not enough disk space for indexing even after cleanup."
        exit 1
    fi
fi



# Function to check if samtools index file exists and is valid
samtools_index_valid() {
    local REF=$1
    if [ ! -s "${REF}.fai" ]; then
        return 1  # Index file is missing or empty
    fi
    return 0  # Index file exists and is not empty
}

# Always attempt to verify index integrity
echo "Checking integrity of BWA index files..."
if ! bwa_index_valid "$REF_FA"; then
    echo "BWA index files missing, incomplete, or potentially corrupted."
    echo "Removing any existing BWA index files and re-indexing..."
    # Remove any existing index files
    rm -f ${REF_FA}.amb ${REF_FA}.ann ${REF_FA}.bwt ${REF_FA}.pac ${REF_FA}.sa

    # Re-index the reference
    echo "Indexing reference genome with BWA..."
    bwa index $REF_FA

    if [ $? -ne 0 ] || ! bwa_index_valid "$REF_FA"; then
        echo "Error: BWA indexing failed or produced incomplete index files."
        exit 1
    fi
else
    echo "BWA index files appear to be complete."
fi

echo "Checking integrity of samtools index file..."
if ! samtools_index_valid "$REF_FA"; then
    echo "Samtools index file missing or potentially corrupted."
    echo "Removing any existing samtools index file and re-indexing..."
    # Remove any existing index file
    rm -f ${REF_FA}.fai

    # Re-index the reference
    echo "Indexing reference genome with samtools..."
    samtools faidx $REF_FA

    if [ $? -ne 0 ] || ! samtools_index_valid "$REF_FA"; then
        echo "Error: Samtools indexing failed or produced an incomplete index file."
        exit 1
    fi
else
    echo "Samtools index file appears to be complete."
fi

# Add a verification step to test if the index is functional
echo "Verifying BWA index functionality..."
echo -e ">test_seq\nACGTACGTACGT" > $REF/test.fa
if ! bwa mem $REF_FA $REF/test.fa > /dev/null 2>&1; then
    echo "BWA index verification failed. Re-indexing reference genome..."
    rm -f ${REF_FA}.amb ${REF_FA}.ann ${REF_FA}.bwt ${REF_FA}.pac ${REF_FA}.sa
    bwa index $REF_FA

    if [ $? -ne 0 ]; then
        echo "Error: BWA re-indexing failed."
        exit 1
    fi
else
    echo "BWA index verification successful."
fi
rm -f $REF/test.fa

touch "$REF_FA.indexed"


echo "=== STEP 6: Aligning reads ==="
# Check disk space before alignment (BWA mem needs significant space)
if ! check_disk_space 20000; then
    echo "WARNING: Low disk space detected before alignment. This step requires significant space."
    # Aggressive cleanup before BWA mem
    cleanup_files

    # Remove original downloaded files if trimmed files exist
    if [ -f "$TRIMMED/trimmed_1P.fastq.gz" ] && [ -f "$TRIMMED/trimmed_2P.fastq.gz" ]; then
        echo "Removing original downloaded files to save space..."
        rm -f $FILE_R1 $FILE_R2
    fi

    # Check space again after aggressive cleanup
    if ! check_disk_space 20000; then
        echo "ERROR: Not enough disk space for alignment even after aggressive cleanup."
        echo "Current system disk space:"
        df -h
        exit 1
    fi
fi

# Align reads using the trimmed FASTQ files
echo "Aligning reads using BWA..."

# Try running BWA with increased file size limit and reduced thread count if needed
echo "Setting resource limits for BWA..."
# Increase file size limit if possible (may require root)
ulimit -f unlimited 2>/dev/null || echo "Could not increase file size limit (not running as root)"

# Run BWA with careful error handling
set +e  # Don't exit on error yet, we want to handle it ourselves
bwa mem -t 4 -R "@RG\tID:NGS01\tSM:NGS01\tPL:ILLUMINA" $REF_FA $TRIMMED/trimmed_1P.fastq.gz $TRIMMED/trimmed_2P.fastq.gz > $SAM
BWA_EXIT=$?
set -e  # Turn exit on error back on

# Check if BWA alignment was successful
if [ $BWA_EXIT -ne 0 ] || [ ! -s "$SAM" ]; then
    echo "Error: BWA alignment failed with exit code $BWA_EXIT or produced an empty SAM file."
    echo "Checking system resources..."
    free -h
    df -h

    # Try again with reduced thread count
    echo "Retrying BWA alignment with reduced thread count..."
    cleanup_files  # Clean up first to maximize available space

    # Try with 2 threads instead of 4
    bwa mem -t 2 -R "@RG\tID:NGS01\tSM:NGS01\tPL:ILLUMINA" $REF_FA $TRIMMED/trimmed_1P.fastq.gz $TRIMMED/trimmed_2P.fastq.gz > $SAM

    if [ $? -ne 0 ] || [ ! -s "$SAM" ]; then
        echo "Error: BWA alignment retry failed. Cannot proceed."
        exit 1
    fi
else
    echo "BWA alignment completed successfully!"
fi

# After successful alignment, remove trimmed fastq files to save space
echo "Removing unpaired trimmed fastq files to save space..."
rm -f $TRIMMED/trimmed_1U.fastq.gz $TRIMMED/trimmed_2U.fastq.gz


echo "=== Cleanup: Removing trimmed FASTQ files to save space before BAM conversion ==="

echo "TRIMMED directory: $TRIMMED"
echo "Checking files for deletion..."

if [ -f "$TRIMMED/trimmed_1P.fastq.gz" ]; then
    echo "Deleting $TRIMMED/trimmed_1P.fastq.gz"
    rm -f "$TRIMMED/trimmed_1P.fastq.gz"
else
    echo "trimmed_1P.fastq.gz not found."
fi

if [ -f "$TRIMMED/trimmed_2P.fastq.gz" ]; then
    echo "Deleting $TRIMMED/trimmed_2P.fastq.gz"
    rm -f "$TRIMMED/trimmed_2P.fastq.gz"
else
    echo "trimmed_2P.fastq.gz not found."
fi

echo "Files in $TRIMMED after cleanup:"
ls -lh "$TRIMMED"






echo "=== STEP 7: Processing alignment files ==="
# Check disk space before BAM conversion
if ! check_disk_space 10000; then
    cleanup_files
    if ! check_disk_space 10000; then
        echo "ERROR: Not enough disk space for BAM conversion."
        exit 1
    fi
fi

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -b $SAM > $BAM

if [ ! -s "$BAM" ]; then
    echo "Error: SAM to BAM conversion failed or produced an empty file."
    exit 1
fi

# After successful conversion to BAM, remove SAM file
echo "Removing SAM file to save space..."
rm -f $SAM

echo "Sorting BAM file..."
samtools sort $BAM > $SORT_BAM

if [ ! -s "$SORT_BAM" ]; then
    echo "Error: BAM sorting failed or produced an empty file."
    exit 1
fi

# After successful sorting, remove unsorted BAM
echo "Removing unsorted BAM file to save space..."
rm -f $BAM

echo "Indexing sorted BAM file..."
samtools index $SORT_BAM

echo "=== STEP 8: Marking duplicates ==="
# Check disk space before marking duplicates
if ! check_disk_space 10000; then
    cleanup_files
    if ! check_disk_space 10000; then
        echo "ERROR: Not enough disk space for marking duplicates."
        exit 1
    fi
fi

# Mark duplicates
echo "Marking duplicate reads..."
picard MarkDuplicates I=$SORT_BAM O=$MARK_BAM M=$ALIGNED/marked_dup_metrics.txt

if [ $? -ne 0 ] || [ ! -s "$MARK_BAM" ]; then
    echo "Error: Marking duplicates failed or produced an empty file."
    exit 1
fi

echo "Indexing marked BAM file..."
samtools index $MARK_BAM

# disk space check
echo "Current disk space usage:"
df -h

echo "Output files:"
echo "- Marked duplicates BAM file: $MARK_BAM"


echo "=== STEP 9: Filtering BAM file ==="
# Check disk space before filtering
if ! check_disk_space 5000; then
    cleanup_files
    if ! check_disk_space 5000; then
        echo "ERROR: Not enough disk space for BAM filtering even after cleanup."
        exit 1
    fi
fi

# Filter BAM file to keep only properly paired reads with mapping quality > 20
echo "Filtering BAM to keep only high-quality reads..."
samtools view -F 0x4 -q 20 -b $MARK_BAM > $FILTER_BAM

if [ $? -ne 0 ] || [ ! -s "$FILTER_BAM" ]; then
    echo "Error: BAM filtering failed or produced an empty file."
    exit 1
fi



FILTERED_BAM="/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_filtered.bam"

# Check the number of reads in the filtered BAM file
echo "Number of reads in the filtered BAM file:"
samtools view -c "$FILTERED_BAM" > "$RESULTS/NGS01_filtered_read_count.txt"
echo "Number of reads outputted to $RESULTS/NGS01_filtered_read_count.txt"

# Check the alignment statistics of the filtered BAM file
echo "Alignment statistics for filtered BAM file:"
samtools flagstat "$FILTERED_BAM" > "$RESULTS/NGS01_filtered_flagstat.txt"
echo "Alignment statistics outputted to $RESULTS/NGS01_filtered_flagstat.txt"

# Check the basic coverage information (first 10 lines)
echo "Basic coverage information for filtered BAM file (first 10 lines):"
samtools depth "$FILTERED_BAM" | head -n 10 > "$RESULTS/NGS01_filtered_depth_head.txt"
echo "Basic coverage information outputted to $RESULTS/NGS01_filtered_depth_head.txt"

# Check index statistics
echo "Index statistics for filtered BAM file:"
samtools idxstats "$FILTERED_BAM" > "$RESULTS/NGS01_filtered_idxstats.txt"
echo "Index statistics outputted to $RESULTS/NGS01_filtered_idxstats.txt"


# Index the filtered BAM file
echo "Indexing filtered BAM file..."
samtools index $FILTER_BAM

# After successful filtering, remove intermediate files to save space
#echo "Removing intermediate files to save space..."
#rm files here


echo "=== STEP 10: Variant Calling using FreeBayes ==="

REF_FA="/home/ubuntu/bioinformatics_course/scripts/data/reference/hg19.fa"
VCF="/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01.vcf"
FILTER_BAM="/home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_filtered.bam"


# Check if reference genome is indexed
echo "Checking if reference genome is indexed:"
if [ ! -f /home/ubuntu/bioinformatics_course/scripts/data/reference/hg19.fa.fai ]; then
  echo "Error: Reference genome is not indexed. Please index the reference."
  exit 1
fi

# Check if the input BAM file is sorted and indexed
echo "Checking if input BAM file is sorted and indexed:"
if [ ! -f /home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_sorted.bam ]; then
  echo "Error: Sorted BAM file not found."
  exit 1
fi
if [ ! -f /home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_sorted.bam.bai ]; then
  echo "Error: BAM file is not indexed. Please index the BAM file."
  exit 1
fi

# Print BAM file header to confirm the data looks fine
echo "Printing BAM file header:"
samtools view -H /home/ubuntu/bioinformatics_course/scripts/data/aligned_data/NGS01_sorted.bam | head -n 10




if [ -z "$REF_FA" ] || [ -z "$VCF" ] || [ -z "$FILTER_BAM" ]; then
    echo "Error: One or more required variables are not set:"
    echo "  REF_FA = $REF_FA"
    echo "  VCF = $VCF"
    echo "  FILTER_BAM = $FILTER_BAM"
    exit 1
fi


# Validate critical files
if [ ! -f "$REF_FA" ]; then
    echo "Error: Reference genome file ($REF_FA) does not exist."
    exit 1
fi

if [ ! -f "$FILE_ANNOTATION" ]; then
    echo "Error: Annotation file ($FILE_ANNOTATION) does not exist."
    exit 1
fi

if [ ! -f "$FILTER_BAM" ]; then
    echo "Error: Filtered BAM file ($FILTER_BAM) does not exist."
    exit 1
fi

# Validate BED file contents
echo "Validating BED file contents..."
head -n 5 "$FILE_ANNOTATION"

# First, try a simple FreeBayes call without region constraints
echo "Attempting basic FreeBayes variant calling..."
freebayes -f "$REF_FA" -v "$VCF" "$FILTER_BAM"

# If basic call fails, try with target regions
if [ $? -ne 0 ]; then
    echo "Basic FreeBayes call failed. Attempting with target regions..."

    # Try with -t instead of -L
    freebayes -f "$REF_FA" -t "$FILE_ANNOTATION" -v "$VCF" "$FILTER_BAM"

    # If that fails, try alternative tools
    if [ $? -ne 0 ]; then
        echo "FreeBayes with target regions failed. Using alternative..."

        # Fallback to bcftools
        bcftools mpileup -f "$REF_FA" "$FILTER_BAM" | bcftools call -mv -Ov > "$VCF"
    fi
fi

# Final check on VCF file
if [ ! -s "$VCF" ]; then
    echo "Error: No variants were called. Please check input files and tool configurations."
    exit 1
fi

echo "Variant calling completed. Output VCF: $VCF"



echo "=== STEP 11: Filtering Variants ==="

# Perform quality filtering on the variants using your preferred filtering criteria
# Here, we filter out variants with a quality score below 20 and coverage less than 10
echo "Filtering variants based on quality and coverage..."
vcffilter -f "QUAL > 20 & DP > 10" $VCF > $FILTER_VCF

# Check if filtering was successful
if [ $? -ne 0 ] || [ ! -s "$FILTER_VCF" ]; then
    echo "Error: Variant filtering failed or filtered VCF file is empty."
    exit 1
fi

echo "Filtered variants saved to: $FILTER_VCF"



 echo "=== STEP 12: Variant Annotation using ANNOVAR ==="


ANNOVAR_DIR=~/annovar
HUMANDB_DIR=$ANNOVAR_DIR/humandb
mkdir -p $HUMANDB_DIR

# Function to check and optionally download an ANNOVAR database
check_or_download_db() {
    local db=$1
    local buildver="hg19"
    local db_file_txt="$HUMANDB_DIR/${buildver}_${db}.txt"
    local db_file_gz="$HUMANDB_DIR/${buildver}_${db}.txt.gz"

    if [ -f "$db_file_txt" ] || [ -f "$db_file_gz" ]; then
        echo "✔️  $db database already exists."
        echo "$db" >> $LOGS/available_dbs.txt
    else
        echo "⬇️  Downloading $db database..."
        $ANNOVAR_DIR/annotate_variation.pl -downdb -buildver $buildver -webfrom annovar $db $HUMANDB_DIR
        if [ $? -eq 0 ]; then
            echo "✅ $db downloaded successfully."
            echo "$db" >> $LOGS/available_dbs.txt

            # Cleanup immediately after successful download
            echo "🧹 Cleaning up temporary files for $db..."
            rm -f "$HUMANDB_DIR/${buildver}_${db}.*gz" # Remove any compressed files.
            rm -f "$HUMANDB_DIR/${buildver}_${db}.*.tmp" # remove any tmp files.
            echo "🧹 Cleanup for $db completed."

        else
            echo "❌ Failed to download $db. Skipping..."
        fi
    fi
}

# Clean previous db log and prepare list of requested databases
rm -f $LOGS/available_dbs.txt
DATABASES=("refGene" "exac03" "avsnp150")

# Check disk space before download
check_disk_space 5000

# Check/download all requested databases
for db in "${DATABASES[@]}"; do
    check_or_download_db "$db"
done

# Read available databases from log
AVAILABLE_DBS=($(cat $LOGS/available_dbs.txt))
if [ ${#AVAILABLE_DBS[@]} -eq 0 ]; then
    echo "Error: No ANNOVAR databases are available for annotation."
    exit 1
fi

# Create protocol and operation strings based on available DBs
PROTOCOL=$(IFS=,; echo "${AVAILABLE_DBS[*]}")
OPERATION=""
for db in "${AVAILABLE_DBS[@]}"; do
    case $db in
        refGene) OPERATION+="g," ;;
        exac03|avsnp150) OPERATION+="f," ;;
        *) OPERATION+="f," ;;
    esac
done
OPERATION=${OPERATION%,}  # Remove trailing comma

# Run ANNOVAR annotation
echo "Running ANNOVAR annotation with: $PROTOCOL"
$ANNOVAR_DIR/table_annovar.pl $FILTER_VCF $HUMANDB_DIR -buildver hg19 -out $RESULTS/NGS01_annovar -remove -protocol $PROTOCOL -operation $OPERATION -nastring . -vcfinput

if [ $? -ne 0 ]; then
    echo "Error: ANNOVAR annotation failed."
    exit 1
fi

echo "Variants annotated using ANNOVAR. Output files located in: $RESULTS/"

# Delete ANNOVAR human database directory before snpEff run to free up space.
if [ -d "$HOME/annovar/humandb" ]; then
    echo "Deleting existing ANNOVAR human database to free up space."
    rm -rf "$HOME/annovar/humandb"
fi


echo "=== STEP 13: Variant Annotation using snpEff ==="

# Function to check and install Java
install_java() {
    if ! command -v java &> /dev/null; then
        echo "Java not found. Installing Java..."
        sudo apt-get update -y
        sudo apt-get install -y default-jre
    else
        echo "Java is already installed."
    fi
}

# Function to download and install snpEff
install_snpeff() {
    if [ ! -d "snpEff" ]; then
        echo "snpEff not found. Downloading and installing snpEff..."
        wget -q https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
        unzip -o snpEff_latest_core.zip
        rm snpEff_latest_core.zip
    else
        echo "snpEff is already installed."
    fi
}




# Function to annotate variants using snpEff
annotate_variants() {
    local genome_version="hg19"
    local filter_vcf="/home/ubuntu/bioinformatics_course/scripts/analysis/NGS01_filtered.vcf"
    local results_dir="$(pwd)/results"
    local snpeff_dir="$(pwd)/snpEff"

    mkdir -p "$results_dir"

    if [ ! -f "$filter_vcf" ]; then
        echo "Error: Input VCF file '$filter_vcf' not found."
        exit 1
    fi

    echo "Annotating variants with snpEff..."
    java -Xmx4g -jar "$snpeff_dir/snpEff.jar" "$genome_version" "$filter_vcf" > "$results_dir/NGS01_snpEff.vcf"

    if [ $? -ne 0 ]; then
        echo "Error: snpEff annotation failed."
        exit 1
    fi

    echo "Variants annotated using snpEff. Output: $results_dir/NGS01_snpEff.vcf"
}

# Function to clean up snpEff installation
cleanup() {
    echo "Removing snpEff directory after use."
    rm -rf "$(pwd)/snpEff"
    echo "snpEff directory removed."
}

# Main script execution
install_java
install_snpeff
annotate_variants
cleanup

echo "=== STEP 13: Variant Annotation using snpEff - Completed ==="




echo "=== STEP 14: Variant Prioritization ==="

# Basic variant prioritization: filter to exonic variants not present in dbSNP
echo "Prioritizing variants: filtering to exonic variants not seen in dbSNP..."
grep -v "rs" $RESULTS/NGS01_annovar.hg19_multianno.txt | awk '$6=="exonic"' > $RESULTS/NGS01_prioritized_variants.txt

if [ $? -ne 0 ] || [ ! -s "$RESULTS/NGS01_prioritized_variants.txt" ]; then
    echo "Error: Variant prioritization failed or output file is empty."
    exit 1
fi

echo "Prioritized variants saved to: $RESULTS/NGS01_prioritized_variants.txt"




