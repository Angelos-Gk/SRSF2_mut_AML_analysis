#!/bin/bash
#SBATCH --job-name=AML1010_maegatk
#SBATCH --output=/path/to/project/AML1010/logs/AML1010_maegatk_%A_%a.out
#SBATCH --error=/path/to/project/AML1010/logs/AML1010_maegatk_%A_%a.err
#SBATCH --partition=medium
#SBATCH --array=1-3
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Run AFTER 01_nanoranger_mito.sh (uses its *_genome_tagged.bam output)

# Activate maegatk environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mito

BASE_DIR="/path/to/project/AML1010"
NANO_BASE="${BASE_DIR}/data/mito/nanoranger"
MAEGATK_BASE="${BASE_DIR}/data/mito/maegatk"
SAMPLES=( "AML1010.1_mito" "AML1010.3_mito" "AML1010.5_mito" )
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_NAME="${SAMPLES[$INDEX]}"
SAMPLE_NANO_DIR="${NANO_BASE}/${SAMPLE_NAME}"
TAGGED_BAM="${SAMPLE_NANO_DIR}/${SAMPLE_NAME}_genome_tagged.bam"
OUTDIR="${MAEGATK_BASE}/${SAMPLE_NAME}"

mkdir -p "${OUTDIR}"

echo "Sample:  ${SAMPLE_NAME}"
echo "BAM:     ${TAGGED_BAM}"
echo "OUTDIR:  ${OUTDIR}"
echo "CPUs:    ${SLURM_CPUS_PER_TASK}"

# Run maegatk
maegatk bcall \
  -i "${TAGGED_BAM}" \
  -o "${OUTDIR}" \
  --ncores "${SLURM_CPUS_PER_TASK}"
