#!/bin/bash
#SBATCH --job-name=wes_alignment
#SBATCH --error=/path/to/project/logs/preprocessing/wes_alignment_%A_%a.err
#SBATCH --output=/path/to/project/logs/preprocessing/wes_alignment_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-27
#SBATCH --mail-user=your.email@institution.org
#SBATCH --mail-type=END,FAIL

# Load environment
source $HOME/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wes

# Paths and config
REF="/path/to/reference/Homo_sapiens_assembly38.masked.fasta"
FASTQ_DIR="/path/to/project/input/fastq"
OUTPUT_DIR="/path/to/project/output"
SORTED_BAM_DIR="${OUTPUT_DIR}/sorted_bams"
MD_OUTPUT_DIR="${OUTPUT_DIR}/mark_dups"
LOGDIR="/path/to/project/logs/preprocessing"
THREADS=16

mkdir -p "$SORTED_BAM_DIR" "$MD_OUTPUT_DIR" "$LOGDIR"

# Use user-provided sample list
SAMPLE_LIST="/path/to/project/input/samplelist.txt"

# Sample name for this array task
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLE_LIST)

START_TIME=$(date +%s)
echo "$(date): Processing sample $SAMPLE"

# Collect all lanes for this sample (anchored match to avoid substring collisions)
mapfile -t R1_FILES < <(
  ls -1 ${FASTQ_DIR}/*_DNA_*_S*_L*_R1_001.fastq.gz | \
  grep -E "_DNA_[0-9]+_${SAMPLE}_S" | \
  sort
)

mapfile -t R2_FILES < <(
  ls -1 ${FASTQ_DIR}/*_DNA_*_S*_L*_R2_001.fastq.gz | \
  grep -E "_DNA_[0-9]+_${SAMPLE}_S" | \
  sort
)

echo "Resolved FASTQ files:"
printf "R1:\n%s\n" "${R1_FILES[@]}"
printf "R2:\n%s\n" "${R2_FILES[@]}"

# Align → Sort → Index
echo "$(date): Starting alignment and sorting"

bwa-mem2 mem -t $THREADS \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}" \
  $REF \
  <(zcat "${R1_FILES[@]}") \
  <(zcat "${R2_FILES[@]}") | \
samtools view -@ $THREADS -b | \
samtools sort -@ $THREADS -o ${SORTED_BAM_DIR}/${SAMPLE}.sorted.bam

samtools index ${SORTED_BAM_DIR}/${SAMPLE}.sorted.bam
echo "$(date): Sorting and indexing complete"

# Mark duplicates
echo "$(date): Marking duplicates"

picard MarkDuplicates \
  I=${SORTED_BAM_DIR}/${SAMPLE}.sorted.bam \
  O=${MD_OUTPUT_DIR}/${SAMPLE}.dedup.bam \
  M=${MD_OUTPUT_DIR}/${SAMPLE}.metrics.txt

samtools index ${MD_OUTPUT_DIR}/${SAMPLE}.dedup.bam

# Runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
SECONDS=$((RUNTIME % 60))

echo "$(date): Total runtime for $SAMPLE: ${HOURS}h ${MINUTES}m ${SECONDS}s"
