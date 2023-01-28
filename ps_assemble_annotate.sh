#!/usr/bin/env bash
#SBATCH -J Pseudomonas_project
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

####Script to assemble, filter and annotate Ps genomes downloaded from "https://microbesng.com/portal/projects/"
####Written in 2022-23 by John Connell 

####Download data, the raw read data was downloaded as below
# screen -S download_data
# srun --partition=short --cpus-per-task=4 --mem=8G --pty bash
# for x in $(cat untrimmed_urls.txt); do
# 	wget ${x}
# done  

####Seperate reads into folders with readname as foldername 
# datadir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data
# for data in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do 
# i=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${data} | cut -c 1-6)
# 	for x in ${i}; do 
# 		cd ${datadir}/${data} 
# 		mkdir -p ${x}
# 		mv ${x}_* ${x}
# 	done 
# done 

####Run fastQC on untrimmed reads
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
# 	for a in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}); do
# 		StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/
# 		ReadsF=$(ls $StrainPath/${strain}/${a}/* | grep "_1.fastq.gz")
#     ReadsR=$(ls $StrainPath/${strain}/${a}/* | grep "_2.fastq.gz")
# 		OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${a}/pre_trim_fastQC
# 		mkdir -p $OutDir
# 		ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# 		sbatch $ProgDir/fasQC_analysis.sh $ReadsF $ReadsR $OutDir
# 	done  
# done 

####Collect fastQC data
# rm /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC/result_table_trim.txt

# for x in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for y in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${x}); do
#     ls -l /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${x} | awk '{print $9}' > filenames


#     ls -l /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC | awk '{print $9 "-F"}' > newnames_f
#     ls -l /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC | awk '{print $9 "-R"}' > newnames_r




# sed -i 1d newnames_f
# sed -i 1d newnames_r
# #Get read counts 
# for name in $(cat filenames); do 
#   for file in /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC/$name/*_1_*; do 
#     cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  f_read_counts_pre
#     cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  f_GC_counts_pre
#     cat $file/summary.txt | cut -f1 | paste -s >> result_f_pre
#   done 
# done 
# for name in $(cat filenames); do 
#   for file in /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC/$name/*_2_*; do 
#     cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  r_read_counts_pre
#     cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  r_GC_counts_pre
#     cat $file/summary.txt | cut -f1 | paste -s >> result_r_pre
#   done 
# done  

# paste newnames_f f_read_counts_pre f_GC_counts_pre result_f_pre >> final_result_f_pre
# paste newnames_r r_read_counts_pre r_GC_counts_pre result_r_pre >> final_result_r_pre
# cat final_result_r_pre >> final_result_f_pre
# echo -e "Read name"'\t'"Read count"'\t'"GC content" > headder1
# cat /home/jconnell/pseudomonas/sequenced_genome_analysis/pre_trim_fastQC/20220819/238496/*_2_*/summary.txt | cut -f2 | paste -s > headder2
# paste headder1  headder2 > headders
# sort -V final_result_f_pre > result_table
# cat result_table >> headders
# mv headders /home/jconnell/pseudomonas/sequenced_genome_analysis/trimmed_QC/result_table_trim.txt

# rm filenames newnames_f newnames_r f_read_counts_pre f_GC_counts_pre result_f_pre result_r_pre r_read_coun








####Trim adapter sequences from raw reads 
# IlluminaAdapters=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
# 	for file in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data/${strain}); do 
# 		StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/raw_data
# 		ReadsF=$(ls $StrainPath/${strain}/${file}/* | grep "_1.fastq.gz")
#    	ReadsR=$(ls $StrainPath/${strain}/${file}/* | grep "_2.fastq.gz")
# 	  F_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${file}/trimmed_reads/F
#   	R_outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${file}/trimmed_reads/R
#   	mkdir -p $F_outdir $R_outdir
#   	ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#   	sbatch $ProgDir/fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters $F_outdir $R_outdir
#  	done   
# done

####Assemble genomes with SPAdes
#  for strain in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for x in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}); do
#       StrainPath=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data
#       F_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/F/*_1_*)
#       R_Read=$(ls $StrainPath/${strain}/${x}/trimmed_reads/R/*_2_*)
#       OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/
#       mkdir -p $OutDir
#       Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
#       while [ $Jobs -gt 50 ]; do
#         sleep 30s
#         printf "."
#         Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
#       done
#       printf "\n"
#       ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#       sbatch $ProgDir/spades.sh $F_Read $R_Read $OutDir correct 10
#   done
# done 

####Rename genomes by stain name 
# for x in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for y in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}); do
#     file=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}/${y}/spades/filtered_contigs/contigs_min_500bp.fasta
#     renamedir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}/${y}/spades/filtered_contigs
#     mkdir -p ${renamedir}/renamed_files_contigs
#     cp ${file} ${renamedir}/renamed_files_contigs/"$y"_contigs_unmasked.fa
#   done
# done 

####Rename contigs as per NCBI
 # for strain in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
 #  for x in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}); do
 #    Assembly=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/filtered_contigs/renamed_files_contigs/*.fa
 #    touch tmp.csv
 #    OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs
 #    mkdir -p $OutDir
 #    name=$(basename $Assembly)
 #    ProgDir=/home/jconnell/johnc/git_repos/niab_repos/pipeline_canker_cherry/cherry_canker_pipeline
 #    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$name" --coord_file tmp.csv
 #    rm tmp.csv
 #   done
 #  done  
 

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

####Checkm on Ps control genomes split
# data=$(cat /home/jconnell/pseudomonas/2795_busco.txt | sed 's/%//g' | awk '($2>=99) {print $1}')
# for x in $(echo $data); do
#   echo $x >> data2
# done 
# cat data2 |  sed -n 1,500p > ps_ref_checkm_1_500.txt
# cat data2 |  sed -n 501,1000p > ps_ref_checkm_501_1000.txt
# cat data2 |  sed -n 1001,1500p > ps_ref_checkm_1001_1500.txt
# cat data2 |  sed -n 1501,2000p > ps_ref_checkm_1501_2000.txt
# cat data2 |  sed -n 2001,2529p > ps_ref_checkm_2001_2529.txt
# for x in ps_ref_checkm_1_500.txt ps_ref_checkm_501_1000.txt  ps_ref_checkm_1001_1500.txt ps_ref_checkm_1501_2000.txt ps_ref_checkm_2001_2529.txt; do 
#     files=/home/jconnell/pseudomonas/Ps_genomes_2795/ncbi_dataset/data
#     data=$(basename $x .txt)
#     temp_files=/mnt/shared/scratch/jconnell/temp_"$data"
#     outdir=/home/jconnell/pseudomonas/Ps_genomes_2795/batch_checkM/$data
#     mkdir -p $temp_files $outdir
#     for n in $(cat $x); do
#         cp --symbolic ${files}/${n}/*.fna $temp_files
#     done 
#     scriptdir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#     sbatch $scriptdir/checkM.sh $temp_files $outdir
# done 

####Filter checkM on control genomes 
# data=/home/jconnell/pseudomonas/Ps_genomes_2795/ncbi_dataset/checkM/storage/edited_data_for_filter.txt
# progdir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# outdir=/home/jconnell/pseudomonas/Ps_genomes_2795/ncbi_dataset/checkM/storage/filterd_data.txt
# python $progdir/subset_checkm.py -i $data -o $outdir


# ####Run Busco on assembled query genomes
# for strain in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for x in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}); do
#     Assembly=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/spades/filtered_contigs/renamed_files_contigs/*.fa
#     BuscoDB=/home/jconnell/pseudomonas/busco_db/bacteria_odb10
#     OutDir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${strain}/${x}/busco/
#     mkdir -p $OutDir
#     ProgDir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
#     sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
#   done
# done 

###Collect Busco data
# for x in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for y in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}); do
#      file=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}/${y}/busco/"$y"_contigs_unmasked
#      cat ${file}/*.txt | sed -n 9p | awk '{print $1}' | cut -c 3-8 >> completeness 
#      cat ${file}/*.txt | sed -n 10p | awk '{print $1}' >> complete_buscos 
#      cat ${file}/*.txt | sed -n 11p | awk '{print $1}' >> complete_single_copy_buscos
#      cat ${file}/*.txt | sed -n 12p | awk '{print $1}' >> complete_and_duplicated_buscos 
#      cat ${file}/*.txt | sed -n 13p | awk '{print $1}' >> fragmented_buscos 
#      cat ${file}/*.txt | sed -n 14p | awk '{print $1}' >> missing_buscos 
#      cat ${file}/*.txt | sed -n 15p | awk '{print $1}' >> total_buscos_searched 
#      echo $y >> file_names
#   done 
# done
# echo -e "Genome""\t""Busco % complete""\t""Complete BUSCOs""\t""Complete and single-copy BUSCOs""\t""Complete and duplicated BUSCOs""\t""Fragmented BUSCOs""\t""Missing BUSCOs""\t""Total BUSCO groups searched" > headder
# paste file_names completeness complete_buscos complete_single_copy_buscos complete_and_duplicated_buscos fragmented_buscos missing_buscos total_buscos_searched > res1
# cat res1 | sort -k2 -n >> headder 
# mv headder 49_56_63_70_77_genomes_busco.txt
# rm completeness complete_buscos complete_single_copy_buscos complete_and_duplicated_buscos fragmented_buscos missing_buscos total_buscos_searched file_names res1 

####Run checkM on query Ps genomes with busco filter 
# data=$(cat /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/49_56_63_70_77_genomes_busco.txt | sed 's/%//g' | awk '($2>=99) {print $1}' | grep -v "Genome")
# files=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/*
# temp_files=/mnt/shared/scratch/jconnell/temp_checkm_450_query
# outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/49_56_63_70_77_checkM
# mkdir -p $temp_files $outdir
# for x in $(echo $data); do
#   cp --symbolic ${files}/${x}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs/*.fa $temp_files
# done 
# scriptdir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# sbatch $scriptdir/checkM.sh $temp_files $outdir

####Blast query genomes - set up master file for all genomes 
# for x in 241278_241370 241185_241277 241464_241556 241371_241463 241557_241649; do
#   for y in $(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}); do
#     genomes=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/${x}/${y}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs/*.fa
#     name=$(basename $y _contigs_unmasked.fa)
#     outfile=/home/jconnell/pseudomonas/blastmaster
#     echo ">"${name} >> ${outfile}
#     samtools faidx ${genomes} contig_1 >> ${outfile}
#   done 
# done 
# cat ${outfile} | grep -v ">contig_1" > /home/jconnell/pseudomonas/blastmaster.txt
# rm /home/jconnell/pseudomonas/blastmaster

####Run Blastn
# data=/home/jconnell/pseudomonas/blastmaster.txt
# outdir=/home/jconnell/blast_query 
# mkdir -p $outdir
# blastdb=/home/jconnell/niab/pseudomonas/blast/prokaryote/ref_prok_rep_genomes
# progdir=/home/jconnell/git_repos/niab_repos/pseudomonas_analysis
# sbatch $progdir/blast.sh $data $outdir $blastdb

### Run pyani on reference genomes 
