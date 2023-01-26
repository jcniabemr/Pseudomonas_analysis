#!/usr/bin/env bash
#SBATCH -J busco
#SBATCH --partition=short
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4


##########################################################################
#INPUT:
# 1st argument: Genome assembly
# 2nd argument: Database option  
#OUTPUT:
# Busco result


# NOTE to easily prepare a figure with a comparison of complete, fragmented, duplicated
# and missing BUSCO genes in each genome, use the script BUSCO_plot.py in the BUSCO folder.
# Instructions in the user guide BUSCO_v2.0_userguide.pdf
###



Assembly=$1
DatabaseOpt=$2
OutDir=$3

Filename=$(basename $Assembly .fna)
#outname="${Filename%.*}"
#extension="${filename##*.}"


WorkDir=/home/jconnell/pseudomonas/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $Assembly $WorkDir
cp -r $DatabaseOpt $WorkDir	
cd $WorkDir


busco=/home/jconnell/miniconda3/envs/busco/bin/busco
$busco \
 -i $Assembly \
 -l $DatabaseOpt \
 -m genome \
 -c 8 \
 -o $Filename


cp -r $Filename $OutDir
rm -r $WorkDir



