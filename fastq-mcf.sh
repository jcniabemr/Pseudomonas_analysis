#!/usr/bin/env bash
#SBATCH -J FastQ-mcf
#SBATCH --partition=medium
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10


########################################################################
#Input 
# 1st argument: Forward read 
# 2nd argument: Reverse read  
# 3rd argument: Illumina adapters
# 4rd argument: Output directory 
#Output
# Will filter poor quality reads, perform trimming and
# remove illumina adapters.
# rna_qc_fastq-mcf <RNASeq_F.fq> <RNASeq_R.fq> <illumina_adapters.fa> [DNA/RNA]


F_read=$1
R_read=$2
illumina_adapters=$3
F_outdir=$4
R_outdir=$5

WorkDir=/mnt/shared/scratch/jconnell/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
mkdir -p $WorkDir/F
mkdir -p $WorkDir/R

cp $F_read $WorkDir
cp $R_read $WorkDir
cp $illumina_adapters $WorkDir
cd $WorkDir


F_file=$(echo $F_read | rev | cut -d "/" -f1 | rev | sed 's/.gz//')
R_file=$(echo $R_read | rev | cut -d "/" -f1 | rev | sed 's/.gz//')

F_out=$(echo "$F_file" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')
R_out=$(echo "$R_file" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')

cat $F_read | gunzip > read1
cat $R_read | gunzip > read2

fastq-mcf \
$illumina_adapters \
read1 \
read2 \
-o $WorkDir/F/"$F_out" \
-o $WorkDir/R/"$R_out" \
-C 1000000 \
-u \
-k 20 \
-t 0.01 \
-q 30 \
-p 5 


gzip $WorkDir/F/"$F_out"
gzip $WorkDir/R/"$R_out"
cp -r $WorkDir/F/* $F_outdir
cp -r $WorkDir/R/* $R_outdir
rm -r $WorkDir
