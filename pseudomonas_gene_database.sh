#!/usr/bin/env bash 
#SBATCH -p long 
#SBATCH -J PSGDB
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

####Script to document the commands used to create Ps gene database 

####Download DC3000 P. syringae pv. tomato DC3000 genome (NCBI accession no. NC_004578)
source miniconda3/bin/activate ncbi_datasets 
mkdir -p pseudomonas/gene_database/DC3000
cd pseudomonas/gene_database/DC3000
datasets download genome accession "GCF_000007805.1" --include genome,gff3,cds
unzip ncbi_dataset.zip
rm ncbi_dataset.zip

###Create NCBI database from DC3000 CDS
source miniconda3/bin/activate blast 
cds=/home/jconnell/pseudomonas/gene_database/DC3000/ncbi_dataset/data/GCF_000007805.1
makeblastdb -in $cds/cds_from_genomic.fna -parse_seqids -blastdb_version 5  -out $cds/DC3000_DB/DC3000 -dbtype nucl 

###Filter Pss genomes for blast search 
for x in $(cat 49_56_63_70_77_passing_filter_IDres.txt | awk '($3 == "syringae") {print $1}' | cut -d "_" -f1); do 
	genomes=/home/jconnell/pseudomonas/pseudomonas_syringae_additional_sequencing/pseudomonas_data/*/${x}/spades/filtered_contigs/renamed_files_contigs/ncbi_contigs
	blastn -query ${genomes}/*.fa -db $cds/DC3000_DB/DC3000 -out pseudomonas/gene_database/blast_results/"$x"_hits_table.txt -evalue 1e-5 -outfmt 6
done 

###Create tables of gene matches per strain and concatenate 
for x in $(ls /home/jconnell/pseudomonas/gene_database/blast_results); do
	cat /home/jconnell/pseudomonas/gene_database/blast_results/${x} \
	| awk '{print $2}' | grep -Eo '[0-9]+$' | sort -nu | awk '{print "g"$0}' \
	> /home/jconnell/pseudomonas/gene_database/temp_files/"$x"_dc3000_blast_genes.txt
done 

paste $(for x in $(ls /home/jconnell/pseudomonas/gene_database/temp_files | cut -c 1-6); do \
 echo /home/jconnell/pseudomonas/gene_database/temp_files/"$x"_hits_table.txt_dc3000_blast_genes.txt; done) \
 >> /home/jconnell/pseudomonas/gene_database/Ps_genes_DC3000.txt

for x in $(ls /home/jconnell/pseudomonas/gene_database/temp_files | cut -c 1-6); do 
    echo ${x} >> /home/jconnell/pseudomonas/gene_database/Ps_names
done 
cat Ps_names | paste -s Ps_names.txt
cat Ps_genes_DC3000.txt >> Ps_names.txt	