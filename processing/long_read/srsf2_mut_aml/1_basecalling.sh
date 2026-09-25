#!/bin/bash
#SBATCH --job-name=ont_basecalling    # Job name
#SBATCH --error=/path/to/project/nanoranger_srsf2/logs/basecalling_%j.err    # Error file (%j will be replaced with the job ID)
#SBATCH --output=/path/to/project/nanoranger_srsf2/logs/basecalling_%j.out   # Output file (%j will be replaced with the job ID)
#SBATCH --time=60:00:00               # Time limit (60 hours)
#SBATCH --partition=gpu               # GPU partition
#SBATCH --gres=gpu:tesla:2            # Request 2 GPUs (adjust the GPU type to your cluster)
#SBATCH --ntasks=1                    # Number of tasks (1 in this case)
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=128G                    # Memory per node
#SBATCH --mail-user=your.email@institution.org # Mail to send to
#SBATCH --mail-type=END,FAIL          # Mail reasons

CUDA_LAUNCH_BLOCKING=1

# Set the model and data paths (dorado v0.8.2)
model=/path/to/dorado-0.8.2-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
data=/path/to/project/nanoranger_srsf2/pod5
output=/path/to/project/nanoranger_srsf2/output/long_read_preprocessing/basecalled.bam

# Change to the directory containing the data
cd $data

# Run the basecaller command with the specified GPU
dorado basecaller --kit-name SQK-NBD114-24 $model $data > $output
dorado summary $output > summary.tsv
