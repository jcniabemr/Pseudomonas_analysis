#!/usr/bin/env bash
#SBATCH -J SPAdes
#SBATCH --partition=short 
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4


########################################################################
#Input 
# 1st argument: Forward read 
# 2nd argument: Reverse read  
# 3rd argument: Illumina adapters
# 4rd argument: Output directory 
# 5th 
#Output
# Will filter poor quality reads, perform trimming and
# remove illumina adapters.
# rna_qc_fastq-mcf <RNASeq_F.fq> <RNASeq_R.fq> <illumina_adapters.fa> [DNA/RNA]


R1=$1
R2=$2
OutDir=$3
Correction=$4
Cutoff='auto'
if [ $5 ]; then
  Cutoff=$5
fi

WorkDir=/home/jconnell/pseudomonas/${SLURM_JOB_USER}_${SLURM_JOB_ID}
mkdir -p $WorkDir

SpadesDir=/home/jconnell/miniconda3/pkgs/spades-3.13.0-0/share/spades-3.13.0-0/bin

F_Read=$(basename $R1)
R_Read=$(basename $R2)

cp $R1 $WorkDir/$F_Read
cp $R2 $WorkDir/$R_Read

$SpadesDir/spades.py -k 21,33,55,77,99,127 -m 200 --phred-offset 33 --careful -1 $WorkDir/$F_Read -2 $WorkDir/$R_Read -t 30  -o $WorkDir/. --cov-cutoff "$Cutoff"

mkdir -p $WorkDir/filtered_contigs
FilterDir=/home/jconnell/johnc/git_repos/niab_repos/pipeline_canker_cherry/cherry_canker_pipeline
$FilterDir/filter_abyss_contigs.py $WorkDir/scaffolds.fasta 500 > $WorkDir/filtered_contigs/contigs_min_500bp.fasta

rm $WorkDir/$F_Read
rm $WorkDir/$R_Read
cp -r $WorkDir/* $OutDir
rm -r $WorkDir

