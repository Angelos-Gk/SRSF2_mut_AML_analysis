#!/bin/bash
#SBATCH --job-name=nanoranger_srsf2
#SBATCH --output=/path/to/project/nanoranger_srsf2/logs/nano/nanoranger_srsf2_%A_%a.out
#SBATCH --error=/path/to/project/nanoranger_srsf2/logs/nano/nanoranger_srsf2_%A_%a.err
#SBATCH --partition=medium
#SBATCH --array=1-14
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Activate Python 3.7 environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python3_7

# Base paths for this project
BASE_DIR="/path/to/project/nanoranger_srsf2"
FASTQ_DIR="${BASE_DIR}/output/nanoranger/fastq"
OUTPUT_BASE="${BASE_DIR}/output/nanoranger"
SCRIPT_DIR="/path/to/nanoranger"   # local clone of the nanoranger GitHub repository
GENOME_REF="/path/to/reference/GRCh38.primary_assembly.genome.fa"
BARCODES_FILE="${SCRIPT_DIR}/data/3M-5pgex-jan-2023.txt.gz"
PANEL_FILE="panel_srsf2.fa"

# Make sure some dirs exist
mkdir -p "${OUTPUT_BASE}"
mkdir -p "${BASE_DIR}/logs/nano"

# Mapping arrays (barcode -> sample ID -> mutation file)
declare -a BARCODE=("barcode11" "barcode12" "barcode13" "barcode14" "barcode15" "barcode16" "barcode17" "barcode18" "barcode19" "barcode20" "barcode21" "barcode22" "barcode23" "barcode24")
declare -a SAMPLE=("P1" "P2" "P3" "P4" "P6" "P7" "P10" "P11" "P12" "P13" "P14" "P15" "P16" "P17")
declare -a MUTATION=("P1.csv" "P2.csv" "P3.csv" "P4.csv" "P6.csv" "P7.csv" "P10.csv" "P11.csv" "P12.csv" "P13.csv" "P14.csv" "P15.csv" "P16.csv" "P17.csv")

# Get values for this task
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
BC=${BARCODE[$INDEX]}
SAMPLE_NAME=${SAMPLE[$INDEX]}
MUT_FILE=${MUTATION[$INDEX]}

# Paths
FASTQ_FILE="${FASTQ_DIR}/${BC}.fastq.gz"
SAMPLE_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
TAGGED_BAM="${SAMPLE_DIR}/${SAMPLE_NAME}_genome_tagged.bam"
PANEL_PATH="${BASE_DIR}/panels/${PANEL_FILE}"
MUT_PATH="${BASE_DIR}/mutations/${MUT_FILE}"

# Create output directory for sample
mkdir -p "$SAMPLE_DIR"
cd "$SAMPLE_DIR"

echo "Running nanoranger pipeline for sample ${SAMPLE_NAME} (BC=${BC})"
echo "FASTQ: ${FASTQ_FILE}"
echo "Panel: ${PANEL_PATH}"
echo "Mutations: ${MUT_PATH}"

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

# Create pileups output dir (relative!)
mkdir -p pileups

# Run perform_pileup (relative output path!)
echo "Running perform_pileup for ${SAMPLE_NAME}"

python "${SCRIPT_DIR}/scripts/downstream/perform_pileup.py" \
    --b "$TAGGED_BAM" \
    --s "$SAMPLE_NAME" \
    --o "pileups" \
    --m "$MUT_PATH"

echo "Done with ${SAMPLE_NAME}"
