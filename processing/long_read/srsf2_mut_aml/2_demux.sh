#!/bin/bash
#SBATCH --job-name=demux_SRSF2mut_AML # Job name
#SBATCH --error=/path/to/project/nanoranger_srsf2/logs/demux_%j.err
#SBATCH --output=/path/to/project/nanoranger_srsf2/logs/demux_%j.out
#SBATCH --time=10:00:00               # Time limit (10 hours)
#SBATCH --partition=medium
#SBATCH --ntasks=1                    # Number of tasks (1 in this case)
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=128G                    # Memory per node
#SBATCH --mail-user=your.email@institution.org # Mail to send to
#SBATCH --mail-type=END,FAIL          # Mail reasons

# Set the data path
output=/path/to/project/nanoranger_srsf2/output/long_read_preprocessing/basecalled.bam
demux_dir=/path/to/project/nanoranger_srsf2/output/long_read_preprocessing/demux

dorado demux --output-dir $demux_dir --no-classify $output
