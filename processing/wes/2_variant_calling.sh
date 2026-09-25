#!/bin/bash
#SBATCH --job-name=strelka2
#SBATCH --error=/path/to/project/logs/strelka2/strelka2_%A_%a.err
#SBATCH --output=/path/to/project/logs/strelka2/strelka2_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-13
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Load environment
source "$HOME/.bashrc"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate strelka2-env

# Paths and config
BASE_DIR="/path/to/project"
MD_BAM_DIR="${BASE_DIR}/output/mark_dups"
OUT_BASE="${BASE_DIR}/output/strelka2"
REF_FASTA="/path/to/reference/Homo_sapiens_assembly38.masked.fasta"
TARGET_BED="/path/to/target_regions/S33266340_hg38/clean.sorted.bed.gz"
STRELKA_BIN="/path/to/strelka-2.9.10.centos6_x86_64/bin"

PAIR_FILE="${BASE_DIR}/input/strelka_pairs.tsv"

mkdir -p "$OUT_BASE"
mkdir -p "${BASE_DIR}/logs/strelka2"

# Read tumor/normal from 2-column pair file (line = SLURM_ARRAY_TASK_ID + 1)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PAIR_FILE")
TUMOR=$(echo "$LINE" | awk '{print $1}')
NORMAL=$(echo "$LINE" | awk '{print $2}')

TUMOR_BAM="${MD_BAM_DIR}/${TUMOR}.dedup.bam"
NORMAL_BAM="${MD_BAM_DIR}/${NORMAL}.dedup.bam"

OUT_DIR="${OUT_BASE}/${TUMOR}_vs_${NORMAL}"
mkdir -p "$OUT_DIR"

echo "$(date): Running Strelka2"
echo "PAIR_FILE = $PAIR_FILE"
echo "TUMOR     = $TUMOR"
echo "NORMAL    = $NORMAL"
echo "TUMOR_BAM = $TUMOR_BAM"
echo "NORMAL_BAM= $NORMAL_BAM"
echo "OUT_DIR   = $OUT_DIR"

# Run strelka2 somatic workflow
"${STRELKA_BIN}/configureStrelkaSomaticWorkflow.py" \
  --normalBam "$NORMAL_BAM" \
  --tumorBam "$TUMOR_BAM" \
  --referenceFasta "$REF_FASTA" \
  --callRegions "$TARGET_BED" \
  --exome \
  --runDir "$OUT_DIR"

"${OUT_DIR}/runWorkflow.py" -m local -j 8

echo "$(date): Finished Strelka2 for $TUMOR vs $NORMAL"
