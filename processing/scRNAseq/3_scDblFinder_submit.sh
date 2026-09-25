#!/bin/bash
#SBATCH --job-name=scDblFinder
#SBATCH --output=/path/to/project/logs/scDblFinder.out
#SBATCH --error=/path/to/project/logs/scDblFinder.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.org

# Load Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RScluster

echo $(R RHOME)  # Check if R is correctly loaded

# Run R script
Rscript /path/to/repo/preprocessing/scrna/03_scdblfinder.R > /path/to/project/logs/scDblFinder_res.log 2>&1
