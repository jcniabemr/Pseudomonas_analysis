#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

#Script to provide a fastQC analysis of raw and trimmed sequencing reads. 

InfileF=$1
InfileR=$2
OutDir=$3

WORK_DIR=/mnt/shared/scratch/jconnell/${SLURM_JOB_USER}_${SLURM_JOBID}_pseudomonasQC
mkdir -p $WORK_DIR

cp -r $InfileF $WORK_DIR
cp -r $InfileR $WORK_DIR
cd $WORK_DIR	

fastqc *.gz
unzip \*.zip

cp -r $WORK_DIR/*fastqc $OutDir
rm -r $WORK_DIR



