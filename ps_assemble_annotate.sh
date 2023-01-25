#!/usr/bin/env bash
#SBATCH -J Pseudomonas_project
#SBATCH -p short
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
# datadir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data
# for data in 241278_241370 241185_241277 241464_241556 241371_241463; do 
# i=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${data}/ | awk '{print $9}' | grep -v "untrimmed_urls.txt" | cut -c 1-6)
# 	for x in ${i}; do 
# 		cd ${datadir}/${data} 
# 		mkdir -p ${x}
# 		mv ${x}_* ${x}
# 	done 
# done 

###Run fastQC on untrimmed reads
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
# 	i=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain} | awk '{print $9}' |  grep -v "untrimmed_urls.txt")
# 	for a in $i; do 
# 		StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/
# 		ReadsF=$(ls $StrainPath/${strain}/${a}/* | grep "_1.fastq.gz")
#     	ReadsR=$(ls $StrainPath/${strain}/${a}/* | grep "_2.fastq.gz")
# 		OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${a}/pre_trim_fastQC
# 		mkdir -p $OutDir
# 		ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# 		sbatch $ProgDir/fasQC_analysis.sh $ReadsF $ReadsR $OutDir
# 	done  
# done 

####Trim adapter sequences from raw reads 
# IlluminaAdapters=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
# 	for file in $(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain} | awk '{print $9}' | grep -v "untrimmed_urls.txt"); do 
# 		StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data
# 		ReadsF=$(ls $StrainPath/${strain}/${file}/* | grep "_1.fastq.gz")
#    	    ReadsR=$(ls $StrainPath/${strain}/${file}/* | grep "_2.fastq.gz")
# 	   	F_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${file}/trimmed_reads/F
#   		R_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${file}/trimmed_reads/R
#   		mkdir -p $F_outdir $R_outdir
#   		ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#   		sbatch $ProgDir/fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters $F_outdir $R_outdir
#  	done   
# done

####Assemble genomes with SPAdes
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
#   a=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain} | awk '{print $9}')
#   for x in $a; do
#       StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data
#       F_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/F/*_1_*)
#       R_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/R/*_2_*)
#       OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/
#       mkdir -p $OutDir
#       Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
#       while [ $Jobs -gt 20 ]; do
#         sleep 30s
#         printf "."
#         Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
#       done
#       printf "\n"
#       ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#       sbatch $ProgDir/spades.sh $F_Read $R_Read $OutDir correct 10
#   done
# done 

####Run checkM on control Ps genomes 
# data=$(cat /home/jconnell/pseudomonas/2795_busco.txt | sed 's/%//g' | awk '($2>=99) {print $1}')
# files=/home/jconnell/pseudomonas/Ps_genomes_2795/ncbi_dataset/data
# temp_files=/mnt/shared/scratch/jconnell/temp_checkm_2795
# outdir=/home/jconnell/pseudomonas/Ps_genomes_2795/ncbi_dataset/checkM
# mkdir -p $temp_files $outdir
# for x in $(echo $data); do
#  cp --symbolic ${files}/${x}/*.fna $temp_files
# done 
# scriptdir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# sbatch $scriptdir/checkM.sh $temp_files $outdir



for strain in 241278_241370 241185_241277 241464_241556 241371_241463; do
  a=$(ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain} | awk '{print $9}')
  for x in $a; do
      StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data
      F_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/F/*_1_*)
      R_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/R/*_2_*)
      OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/
      mkdir -p $OutDir
      ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
      sbatch $ProgDir/spades.sh $F_Read $R_Read $OutDir correct 10
  done
done 


