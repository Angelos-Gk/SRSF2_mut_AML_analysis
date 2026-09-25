#!/bin/bash
#SBATCH --job-name=bam2fastq_srsf2
#SBATCH --output=/path/to/project/nanoranger_srsf2/logs/bam2fastq_%A_%a.out
#SBATCH --error=/path/to/project/nanoranger_srsf2/logs/bam2fastq_%A_%a.err
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --array=17-24
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Load conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wes

# Define directories
# INPUT_DIR = demultiplexed BAMs from 02_dorado_demux.sh
INPUT_DIR="/path/to/project/nanoranger_srsf2/output/long_read_preprocessing/demux"
OUTPUT_DIR="/path/to/project/nanoranger_srsf2/output/nanoranger/fastq"

# Make sure output dir exists
mkdir -p "$OUTPUT_DIR"

# Format barcode (e.g. 17, 18, ..., 24)
BARCODE=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

# Find the BAM file for this barcode (dorado may prefix the kit name, e.g. SQK-NBD114-24_barcode17.bam)
BAM_FILE=$(ls "${INPUT_DIR}"/*barcode${BARCODE}.bam 2>/dev/null | head -n 1)

if [[ ! -f "$BAM_FILE" ]]; then
    echo "No BAM file found for barcode${BARCODE} in ${INPUT_DIR}"
    exit 1
fi

# Output FASTQ named by barcode
OUTPUT_FASTQ="${OUTPUT_DIR}/barcode${BARCODE}.fastq"

echo "Converting $BAM_FILE to ${OUTPUT_FASTQ}.gz"
bedtools bamtofastq -i "$BAM_FILE" -fq "$OUTPUT_FASTQ"
gzip "$OUTPUT_FASTQ"

echo "Task $SLURM_ARRAY_TASK_ID (barcode${BARCODE}) finished."
