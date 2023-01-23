#!/usr/bin/env bash
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -J edit_salmon

#Edit gene names in quant file to be g1 onwards 

directory=/home/jconnell/andrea_rna_seq/salmon/Pss
files=$(ls -l /home/jconnell/andrea_rna_seq/salmon/Pss | awk '{print $9}')


for x in ${files}; do 
	mkdir -p ${directory}/${x}/edit_data
	cat ${directory}/${x}/quant.sf | cut -f1 | grep -Eo '[0-9]+$' | awk '{print "g" $0}' > ${directory}/${x}/edit_data/gene_names
	cat ${directory}/${x}/quant.sf | cut -f2- > ${directory}/${x}/edit_data/quant.sf
	cd  ${directory}/${x}/edit_data 
	paste gene_names quant.sf > q
	rm gene_names quant.sf
	mv q quant.sf 
done 
