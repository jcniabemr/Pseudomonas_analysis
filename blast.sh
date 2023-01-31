#!/usr/bin/env bash
#SBATCH -J blast
#SBATCH -p long 
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4


infile=$1
outdir=$2
db=$3


WorkDir=/mnt/shared/scratch/jconnell/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp -s $infile $WorkDir
cd $WorkDir
name=$(basename $infile .txt)

blastn -query $infile -db $db -out "$name"_result.txt #-outfmt 6

cp -r "$name"_result.txt $outdir
rm -r $WorkDir



#Script to setup and run blast for species identification 

#Creata conda env for blast and install
#conda create -n blast 
#conda install -c bioconda blast

#Show all available blast databases
#update_blastdb.pl --showall [*]

#Install and decompress required databases 
#update_blastdb.pl --decompress ref_prok_rep_genomes


#Usage for a nucleotide blast of nt.fsa against the database nt, results go to results.out, --remote also searches NCBA servers 
###blastn –db nt –query nt.fsa –out results.out -remote###

# #blastn -query test.txt -db ref_prok_rep_genomes -out myseqs.txt



# #i=$(ls -l /home/jconnell/niab/pseudomonas/genomes/new_sequencing/renamed_final_genomes | awk '{print $9}'| head -n 2)
# # i=$(ls /home/jconnell/pseudomonas/blast_test_databases | grep -w "test")
# # blastdb=/home/jconnell/niab/pseudomonas/blast/prokaryote/ref_prok_rep_genomes
# blastdb_16s=/home/jconnell/niab/pseudomonas/blast/16s

# for x in $(echo $i); do 
# 	cp /home/jconnell/projects/niab/pseudomonas/genomes/new_sequencing/renamed_final_genomes/$x $WorkDir
# 	cd $WorkDir
# 	blastn -query $x -db $blastdb_16s -out "$x"_blast_16s_result.txt -outfmt "6 qseqid sseqid evalue bitscore qstart qend slen length mismatch gapopen gaps sseq"  -word_size 5 -perc_identity 80
# 	cp "$x"_blast_16s_result.txt /home/jconnell/pseudomonas/blast_test_databases/16s_test
# done 

# for x in $(echo $i); do 
# 	cp /home/jconnell/projects/niab/pseudomonas/genomes/new_sequencing/renamed_final_genomes/$x $WorkDir
# 	cd $WorkDir
# 	blastn -query $x -db $blastdb -out "$x"_blast_prok_result.txt -outfmt "6 qseqid sseqid evalue bitscore qstart qend slen length mismatch gapopen gaps sseq"  -word_size 5 -perc_identity 80
# 	cp "$x"_blast_prok_result.txt /home/jconnell/pseudomonas/blast_test_databases/prokaryotic_test
# done 

#rm -r $WorkDir
#cd /home/jconnell/niab/pseudomonas/blast/prokaryote

# for x in $(echo $i); do 
# 	cp /home/jconnell/projects/niab/pseudomonas/genomes/new_sequencing/renamed_final_genomes/$x $WorkDir
# 	cd $WorkDir
# 	blastn -query $x -db $blastdb -out "$x"_blast_16s_result.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# 	cp "$x"_blast_16s2_result.txt /home/jconnell/pseudomonas/blast_test_databases/16s_test
# done 

# for x in $(echo $i); do 
# 	cp /home/jconnell/pseudomonas/blast_test_databases/$x $WorkDir
# 	cd $WorkDir
# 	blastn -query $x -db $blastdb_16s -out "$x"_blast_16s_result.txt 
# 	cp "$x"_blast_16s2_result.txt /home/jconnell/pseudomonas/blast_test_databases/16s_test
# done 