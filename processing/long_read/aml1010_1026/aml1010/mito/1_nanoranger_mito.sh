#!/bin/bash
#SBATCH --job-name=AML1010_mito_nano
#SBATCH --output=/path/to/project/AML1010/logs/AML1010_mito_nano_%A_%a.out
#SBATCH --error=/path/to/project/AML1010/logs/AML1010_mito_nano_%A_%a.err
#SBATCH --partition=medium
#SBATCH --array=1-3
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Activate Python 3.7 environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python3_7

# Base paths for this project
BASE_DIR="/path/to/project/AML1010"
FASTQ_DIR="${BASE_DIR}/data/mito/fastq"
OUTPUT_BASE="${BASE_DIR}/data/mito/nanoranger"
SCRIPT_DIR="/path/to/nanoranger"   # local clone of the nanoranger GitHub repository
BARCODES_FILE="${SCRIPT_DIR}/data/737K-august-2016.txt.gz"
GENOME_REF="/path/to/reference/GRCh38.primary_assembly.genome.fa"
PANEL_FILE="/path/to/project/nanoranger/panel/panel_MT_trns_extended_1010.fa"

# Sample list (0-based indexing)
SAMPLES=( "AML1010.1_mito" "AML1010.3_mito" "AML1010.5_mito" )

# Convert SLURM array task ID (1-based) to bash array index (0-based)
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Select sample for this array task
SAMPLE_NAME="${SAMPLES[$INDEX]}"

# Paths
FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_combined.fastq.gz"
SAMPLE_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
PANEL_PATH="${PANEL_FILE}"

# Create output directory for sample
mkdir -p "$SAMPLE_DIR"
cd "$SAMPLE_DIR"

[[ -f "$FASTQ_FILE" ]] || { echo "Missing FASTQ: $FASTQ_FILE"; exit 1; }
[[ -f "$PANEL_PATH" ]] || { echo "Missing panel: $PANEL_PATH"; exit 1; }
[[ -f "$GENOME_REF" ]] || { echo "Missing genome: $GENOME_REF"; exit 1; }
[[ -f "$BARCODES_FILE" ]] || { echo "Missing barcodes: $BARCODES_FILE"; exit 1; }

echo "Running nanoranger pipeline for sample ${SAMPLE_NAME}"
echo "FASTQ: ${FASTQ_FILE}"
echo "Panel: ${PANEL_PATH}"

# Run nanoranger pipeline
python "${SCRIPT_DIR}/pipeline.py" \
    --c 8 \
    --i "$FASTQ_FILE" \
    --o "$SAMPLE_DIR" \
    --e "$SAMPLE_NAME" \
    --t "$PANEL_PATH" \
    --g "$GENOME_REF" \
    --m 5p10XGEX \
    --barcodes "$BARCODES_FILE"
