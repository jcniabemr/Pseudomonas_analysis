#!/usr/bin/env bash
#SBATCH -J Pseudomonas_project
#SBATCH -p medium
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

#Script to assemble, filter and annotate Ps genomes downloaded from "https://microbesng.com/portal/projects/"
#Written in 2022-23 by John Connell 

####Download data, the raw read data was downloaded as below

#screen -S download_data
#srun --partition=short --cpus-per-task=4 --mem=8G --pty bash
# for x in $(cat untrimmed_urls.txt); do
# 	wget ${x}
# done  

####Seperate reads into folders with readname as foldername 

# datadir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data
# for data in 241278_241370 241185_241277 241464_241556 241371_241463; do 
# i=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${data}/reads | awk '{print $9}' | grep -v "download_data\|run.txt" | cut -c 1-6)
# 	for x in ${i}; do 
# 		cd ${datadir}/${data}/reads 
# 		mkdir -p ${x}
# 		mv ${x}_* ${x}
# 	done 
# done 

####Run fastQC on untrimmed reads
 
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
# 	a=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads | awk '{print $9}')
# 	for i in $a; do
# 		StrainPathF=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads/${i}/ | grep "_1_")
# 		StrainPathR=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads/${i}/ | grep "_2_")
# 		infileF=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads/${i}/$StrainPathF	
# 		infileR=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads/${i}/$StrainPathR	
# 		OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${i}/pre_trim_fastQC
# 		mkdir -p $OutDir
# 		ProgDir=git_repos/niab_repos/pipeline_canker_cherry/cherry_canker_pipeline
# 		sbatch $ProgDir/fastqc_ud.sh $infileF $infileR $OutDir	
# 	done  
# done 

#Trim adapter sequences from raw reads 
IlluminaAdapters=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
	file=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}/untrimmed_reads)
	for x in $file; do
 	StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data
  	ReadsF=$(ls $StrainPath/${strain}/*/${x} | grep "_1.fastq.gz")
  	ReadsR=$(ls $StrainPath/${strain}/*/${x} | grep "_2.fastq.gz")
  	echo $ReadsR
  	echo $ReadsF
  done 
done 
#    	F_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/trimmed_reads/F
#   	R_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/trimmed_reads/R
#   	mkdir -p $F_outdir $R_outdir
#   	ProgDir=git_repos/niab_repos/pipeline_canker_cherry/cherry_canker_pipeline
#   	sbatch $ProgDir/DNA_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters $F_outdir $R_outdir
#  	done   
# done




