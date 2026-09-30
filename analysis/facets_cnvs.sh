#!/bin/bash
#SBATCH --job-name=cnv_facets_srsf2_mut_aml
#SBATCH --error=/path/to/project/logs/cnv_facets/cnv_facets_srsf2_mut_aml_%A_%a.err
#SBATCH --output=/path/to/project/logs/cnv_facets/cnv_facets_srsf2_mut_aml_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-13
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

conda activate facets

# dbSNP VCF (shipped with the Funcotator data sources) and exome capture targets (as in 02_strelka2_somatic.sh)
DBSNP_VCF="/path/to/funcotator_dataSources.v1.8.hg38.20230908s/dbsnp/hg38/hg38_All_20180418.vcf.gz"
TARGETS_BED="/path/to/target_regions/S33266340_hg38/clean.sorted.bed.gz"
# Tab-separated with a header line; columns: tumor_id, normal_id, tumor_bam, normal_bam
# (e.g. P1_tumor  P1_germ  .../mark_dups/P1_tumor.dedup.bam  .../mark_dups/P1_germ.dedup.bam)
PAIRS_FILE="/path/to/project/input/facets_pairs.tsv"
OUT_BASE="/path/to/project/cnv_facets/output"
mkdir -p "$OUT_BASE" "/path/to/project/logs/cnv_facets"

# Read pair for this array task (skip header, 0-indexed)
TUMOR_ID=$(awk  -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+2{print $1}' "$PAIRS_FILE")
NORMAL_ID=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+2{print $2}' "$PAIRS_FILE")
TUMOR_BAM=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+2{print $3}' "$PAIRS_FILE")
NORMAL_BAM=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+2{print $4}' "$PAIRS_FILE")

if [[ -z "$TUMOR_BAM" || -z "$NORMAL_BAM" ]]; then
    echo "ERROR: Could not read pair for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

if [[ ! -f "$TUMOR_BAM" ]]; then
    echo "ERROR: Tumor BAM not found: $TUMOR_BAM" >&2
    exit 1
fi
if [[ ! -f "$NORMAL_BAM" ]]; then
    echo "ERROR: Normal BAM not found: $NORMAL_BAM" >&2
    exit 1
fi

PAIR_TAG="${TUMOR_ID}__vs__${NORMAL_ID}"
OUT_DIR="${OUT_BASE}/${PAIR_TAG}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "$(date) Running FACETS on $PAIR_TAG"
echo "  Tumor BAM : $TUMOR_BAM"
echo "  Normal BAM: $NORMAL_BAM"

cnv_facets.R \
  -t "$TUMOR_BAM" \
  -n "$NORMAL_BAM" \
  -vcf "$DBSNP_VCF" \
  -o "$PAIR_TAG" \
  --snp-mapq 15 \
  --snp-baq 20 \
  --snp-nprocs "${SLURM_CPUS_PER_TASK}" \
  --depth 25 4000 \
  --cval 25 150 \
  --nbhd-snp auto \
  --targets "$TARGETS_BED" \
  --gbuild hg38 \
  --rnd-seed 42

EXIT_CODE=$?
if [[ $EXIT_CODE -ne 0 ]]; then
    echo "ERROR: cnv_facets.R failed for $PAIR_TAG with exit code $EXIT_CODE" >&2
    exit $EXIT_CODE
fi

echo "$(date) Done $PAIR_TAG"
