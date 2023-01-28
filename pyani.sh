#!/usr/bin/env bash
#SBATCH --mem=50G
#SBATCH -p long
#SBATCH -J pyani
#SBATCH --cpus-per-task=12

#Calculate ANI for a set of genome assemblies 


input_files=/home/jconnell/pseudomonas/sequenced_genome_analysis/assembled_genomes/
result_file=/home/jconnell/pseudomonas/results_ani_seq_genomes

average_nucleotide_identity.py \
-i $input_files \
-o $result_file \
-m ANIm \
-g 

