#!/usr/bin/env bash 
#SBATCH -p long 
#SBATCH -J PSGDB
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

####Script to document the commands used to create Ps gene database 

####Download DC3000 P. syringae pv. tomato DC3000 genome (NCBI accession no. NC_004578)
# source miniconda3/bin/activate ncbi_datasets 
# mkdir -p pseudomonas/gene_database/DC3000
# cd pseudomonas/gene_database/DC3000
# datasets download genome accession "GCF_000007805.1" --include genome,gff3,cds
# unzip ncbi_dataset.zip
# rm ncbi_dataset.zip

# ####Create NCBI database from DC3000 CDS
# source activate blast 
# cds=/home/jconnell/pseudomonas/gene_database/DC3000/ncbi_dataset/data/GCF_000007805.1
# makeblastdb -in $cds/cds_from_genomic.fna -parse_seqids -blastdb_version 5  -out $cds/DC3000_DB/DC3000 -dbtype nucl 

# ####Filter Pss genomes for blast search 
# for x in $(cat 49_56_63_70_77_passing_filter_IDres.txt | awk '($3 == "syringae") {print $1}' | cut -d "_" -f1); do 
# 	genomes=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/*/${x}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs
# 	blastn -query ${genomes}/*.fa -db $cds/DC3000_DB/DC3000 -out pseudomonas/gene_database/blast_results/"$x"_hits_table.txt -evalue 1e-5 -outfmt 6
# done 
# ####Filter Pss with additional species added (spacies names)
# for x in $(cat pseudomonas/gene_database/species_names); do
#     data=$(cat 49_56_63_70_77_passing_filter_IDres.txt | cut -f1-2 | grep -w $x | cut -f1 | cut -d "_" -f1)
#     for y in $data; do
#     	genomes=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/*/${y}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs
#     	blastn -query ${genomes}/*.fa -db $cds/DC3000_DB/DC3000 -out pseudomonas/gene_database/blast_results/"$y"_hits_table.txt -evalue 1e-10 -outfmt 6
# 	done
# done  

# ####Control strains 



# ####Create tables of gene matches per strain and concatenate 
# for x in $(ls /home/jconnell/pseudomonas/gene_database/blast_results); do
# 	cat /home/jconnell/pseudomonas/gene_database/blast_results/${x} \
# 	| awk '{print $2}' | grep -Eo '[0-9]+$' | sort -nu | awk '{print "g"$0}' \
# 	> /home/jconnell/pseudomonas/gene_database/temp_files/"$x"_dc3000_blast_genes.txt
# done 

# paste $(for x in $(ls /home/jconnell/pseudomonas/gene_database/temp_files | cut -c 1-6); do \
#  echo /home/jconnell/pseudomonas/gene_database/temp_files/"$x"_hits_table.txt_dc3000_blast_genes.txt; done) \
#  >> /home/jconnell/pseudomonas/gene_database/Ps_genes_DC3000.txt

# for x in $(ls /home/jconnell/pseudomonas/gene_database/temp_files | cut -c 1-6); do 
#     echo ${x} >> /home/jconnell/pseudomonas/gene_database/Ps_names
# done 
# cat Ps_names | paste -s > Ps_names.txt
# cat Ps_genes_DC3000.txt >> Ps_names.txt	
# paste DC3000_genes.txt Ps_names.txt > unordered_DB.txt


# for x in $(cat gene_table_matched.txt | cut -f1 | sed 's/X//g'); do 
# 	cat 49_56_63_70_77_passing_filter_IDres.txt | cut -f1 | grep $x >> realnames
#done 

####Type 3 effectors 
#source activate blast 
# effectors=/home/jconnell/pseudomonas/michelle_effectors/Pss/attempt1/t3es.fasta
# for x in $(cat pseudomonas/gene_database/species_names); do
#     data=$(cat /home/jconnell/49_56_63_70_77_passing_filter_IDres.txt | cut -f1-2 | grep -w $x | cut -f1 | cut -d "_" -f1)
#     for y in $data; do
#     	strain=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/*/${y}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs/
#     	gdb=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/temp_dbs/${y}
#     	outdir=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/effector_blast/${y}
#     	mkdir -p $outdir
# 		makeblastdb -in $strain/*.fa -parse_seqids -blastdb_version 5  -out $gdb/"$y"_db/"$y" -dbtype nucl 
# 		tblastn -query $effectors -db $gdb/"$y"_db/"$y" -out $outdir/"$y"_hits_table -evalue 1e-5 -outfmt 6
# 	done 
# done 

####Filter results 
#FIlter file 1 : unique effectors for each strain 
#Filter file 2 : ≥ 70% identity and ≥ 40% query length
#Filter file 3 : E value,1e^-5, 20% minimum identity 
#header : qseqid sseqid %ident length mismatch gapopen qstart qend sstart send evalue bitscore
#>40% query length 


####Get list of files as headers 

a=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/effector_blast)
echo $a > /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/effector_blast_headers

####Filter 1 unique effectors 
temp_loc=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/temp_effector_lists
mkdir -p $temp_loc
infiles=$(ls /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/effector_blast)
for x in $infiles; do 
	cat /home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/effector_blast/$x/*\
	 | awk '{print $1}' | sort | uniq -c | awk '{print $2}' > $temp_loc/"$x".txt
done 
####Combine lists 
paste $(for x in $(ls ./); do echo "$x"; done) >> ../effector_list1.txt 