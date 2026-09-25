#!/bin/bash
#SBATCH --job-name=strelka2_filter
#SBATCH --error=/path/to/project/logs/strelka2_filtering/strelka2_filter_%A_%a.err
#SBATCH --output=/path/to/project/logs/strelka2_filtering/strelka2_filter_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-13
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Load environment
source "$HOME/.bashrc"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wes

# Paths and config
BASE="/path/to/project"
STRELKA_DIR="${BASE}/output/strelka2"
LOG_DIR="${BASE}/logs/strelka2_filtering"
OUT_CSV="${LOG_DIR}/filter_summary.csv"
FUNCOTATOR_DB="/path/to/funcotator_dataSources.v1.8.hg38.20230908s"
REFERENCE="/path/to/reference/Homo_sapiens_assembly38.masked.fasta"
SCRIPT="/path/to/repo/preprocessing/wes/strelka_vaf_parser.py"
PAIR_FILE="${BASE}/input/strelka_pairs.tsv"

mkdir -p "$LOG_DIR"

# Read tumor/normal for this array task (line = SLURM_ARRAY_TASK_ID + 1)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PAIR_FILE")
TUMOR=$(echo "$LINE" | awk '{print $1}')
NORMAL=$(echo "$LINE" | awk '{print $2}')

PAIR_NAME="${TUMOR}_vs_${NORMAL}"
PAIR_FOLDER="${STRELKA_DIR}/${PAIR_NAME}"
VAR_DIR="${PAIR_FOLDER}/results/variants"

cd "$VAR_DIR"

# === 1) Parse VAF/AD/DP from Strelka VCFs ===
python3 "$SCRIPT" --input somatic.snvs.vcf.gz --output somatic.snvs.parsed.vcf --variant snv
python3 "$SCRIPT" --input somatic.indels.vcf.gz --output somatic.indels.parsed.vcf --variant indel

bgzip -f somatic.snvs.parsed.vcf && tabix -f somatic.snvs.parsed.vcf.gz
bgzip -f somatic.indels.parsed.vcf && tabix -f somatic.indels.parsed.vcf.gz

# === 2) Merge SNVs + indels ===
bcftools concat -a -Oz somatic.snvs.parsed.vcf.gz somatic.indels.parsed.vcf.gz -o merged_parsed.vcf.gz
tabix -f merged_parsed.vcf.gz

# === 3) Sequential filters (tumor is sample index [1]) ===
bcftools view -f PASS -Oz merged_parsed.vcf.gz -o merged_parsed.pass.vcf.gz
tabix -f merged_parsed.pass.vcf.gz

bcftools view -i 'FORMAT/DPVAF[1] >= 50' merged_parsed.pass.vcf.gz -Oz -o merged_parsed.pass.dp50.vcf.gz
tabix -f merged_parsed.pass.dp50.vcf.gz

bcftools view -i 'FORMAT/AD[1:1] >= 5' merged_parsed.pass.dp50.vcf.gz -Oz -o merged_parsed.pass.dp50.ad5.vcf.gz
tabix -f merged_parsed.pass.dp50.ad5.vcf.gz

bcftools view -i 'FORMAT/VAF[1:0] >= 0.05' merged_parsed.pass.dp50.ad5.vcf.gz -Oz -o merged_parsed.pass.dp50.ad5.vaf5.vcf.gz
tabix -f merged_parsed.pass.dp50.ad5.vaf5.vcf.gz

# === 4) Count & log (single pair) ===
count_vars() { bcftools view -H "$1" 2>/dev/null | wc -l; }

echo "${PAIR_NAME},$(count_vars merged_parsed.vcf.gz),$(count_vars merged_parsed.pass.vcf.gz),\
$(count_vars merged_parsed.pass.dp50.vcf.gz),$(count_vars merged_parsed.pass.dp50.ad5.vcf.gz),\
$(count_vars merged_parsed.pass.dp50.ad5.vaf5.vcf.gz)" >> "${OUT_CSV}.${SLURM_ARRAY_TASK_ID}"

# === 5) Annotate with Funcotator (canonical) ===
gatk Funcotator \
    --variant merged_parsed.pass.dp50.ad5.vaf5.vcf.gz \
    --reference "$REFERENCE" \
    --ref-version hg38 \
    --data-sources-path "$FUNCOTATOR_DB" \
    --output merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf \
    --output-file-format VCF

bgzip -f merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf
tabix -f merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz

# === Cleanup ===
rm -f header.vcf body.vcf
echo "Successfully processed ${PAIR_NAME}"

# === 6) Annotate with Funcotator (best effect) ===
#gatk Funcotator \
#    --variant merged_parsed.pass.dp50.ad5.vaf5.vcf.gz \
#    --reference "$REFERENCE" \
#    --ref-version hg38 \
#    --data-sources-path "$FUNCOTATOR_DB" \
#    --output merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf \
#    --output-file-format VCF \
#    --transcript-selection-mode BEST_EFFECT

#bgzip -f merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf
#tabix -f merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
