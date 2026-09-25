#!/bin/bash
#SBATCH --job-name=nanoranger_AML1026
#SBATCH --output=/path/to/project/nanoranger/logs/nanoranger_AML1026_%A_%a.out
#SBATCH --error=/path/to/project/nanoranger/logs/nanoranger_AML1026_%A_%a.err
#SBATCH --partition=medium
#SBATCH --array=1-6
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Activate Nanoranger environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python3_7

# --- Paths ---
BASE_DIR="/path/to/project/nanoranger"
FASTQ_DIR="${BASE_DIR}/fastq"
OUTPUT_BASE="${BASE_DIR}/output/nanoranger"
SCRIPT_DIR="/path/to/nanoranger"   # local clone of the nanoranger GitHub repository
GENOME_REF="/path/to/reference/GRCh38.primary_assembly.genome.fa"
BARCODES_FILE="${SCRIPT_DIR}/data/737K-august-2016.txt.gz"
PANEL_PATH="${BASE_DIR}/panel/panel_MT_trns_extended_1026.fa"
MUT_PATH="${BASE_DIR}/mutations/mutations.1026.csv"

PILEUP_SAMPLE="AML1026"
mkdir -p "${BASE_DIR}/logs/nano"
mkdir -p "${OUTPUT_BASE}"

# Array mapping
declare -a SAMPLE=(
  "AML1026.1_rhPCR"
  "AML1026.1_taqman"
  "AML1026.1_taqman_biotin"
  "AML1026.2_rhPCR"
  "AML1026.3_rhPCR"
  "AML1026.4_rhPCR"
)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_NAME="${SAMPLE[$INDEX]}"

FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_NAME}.combined.fastq.gz"
SAMPLE_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"

mkdir -p "$SAMPLE_DIR"
cd "$SAMPLE_DIR"

echo "Running Nanoranger for ${SAMPLE_NAME}"
echo "FASTQ:   ${FASTQ_FILE}"
echo "PANEL:   ${PANEL_PATH}"
echo "GENOME:  ${GENOME_REF}"
echo "MUT CSV: ${MUT_PATH}"
echo "OUTDIR:  ${SAMPLE_DIR}"

# --- checks ---
[[ -s "$FASTQ_FILE" ]]     || { echo "[ERROR] FASTQ missing: $FASTQ_FILE"; exit 1; }
[[ -s "$PANEL_PATH" ]]     || { echo "[ERROR] Panel FASTA missing: $PANEL_PATH"; exit 1; }
[[ -s "$GENOME_REF" ]]     || { echo "[ERROR] Genome FASTA missing: $GENOME_REF"; exit 1; }
[[ -s "$BARCODES_FILE" ]]  || { echo "[ERROR] Whitelist missing: $BARCODES_FILE"; exit 1; }
[[ -s "$MUT_PATH" ]]       || { echo "[ERROR] Mutation CSV missing: $MUT_PATH"; exit 1; }

# --- Run Nanoranger ---
python "${SCRIPT_DIR}/pipeline.py" \
  --c "${SLURM_CPUS_PER_TASK}" \
  --i "$FASTQ_FILE" \
  --o "$SAMPLE_DIR" \
  --e "$SAMPLE_NAME" \
  --t "$PANEL_PATH" \
  --g "$GENOME_REF" \
  --m 5p10XGEX \
  --barcodes "$BARCODES_FILE"

# --- Pileup ---
mkdir -p pileups
python "${SCRIPT_DIR}/scripts/downstream/perform_pileup.py" \
  --b "${SAMPLE_DIR}/${SAMPLE_NAME}_genome_tagged.bam" \
  --s "$PILEUP_SAMPLE" \
  --o "pileups" \
  --m "$MUT_PATH"

echo "Done with ${SAMPLE_NAME}"
