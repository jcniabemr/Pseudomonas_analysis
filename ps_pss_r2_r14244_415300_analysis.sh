#!/usr/bin/env bash
#SBATCH -p medium
#SABTCH -J RNA_SEQ_Pseudomonas
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

####Script to analyse 4 different pseudomonas strains that underwent RNA Seq

# Pss9644
# RefSeq GCF_023277945.1
# https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_023277945.1/
# RNA_1 RNA_9_2607 RNA_14 RNA_16 RNA_18 RNA_19
####
# R2leaf
# RefSeq GCF_002905795.2
# https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002905795.2/
# RNA_4 RNA_6 RNA_10_2607  RNA_17 RNA_23 RNA_24
####
# R15244
# RefSeq GCF_002905685.2
# https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002905685.2/
# RNA_2 RNA_3 RNA_5_2607 RNA_11 RNA_12 RNA_15
####
# R15300
# RefSeq GCF_002905875.2
# https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002905875.2/
# RNA_7 RNA_8_2607 RNA_13 RNA_20 RNA_21 RNA_22

####Note
#RNA samples in /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria
#ref genomes in /home/jconnell/pseudomonas/genome_data_andrea

####Pre trimming fastQC
# v=$(ls -l /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData | awk '{print $9}')
# for f in $(echo $v); do
#   for RawData in $(ls /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData/${f}/*.fq.gz); do
#       OutDir=/home/jconnell/niab/andrea_rna_seq/${f}
#       mkdir -p $OutDir
#       ProgDir=/home/jconnell/niab/andrea_rna_seq/scripts
#       sbatch $ProgDir/fastqc_ud.sh $RawData $OutDir
#   done
# done  

####Collect QC data
# ls -l /home/jconnell/niab/andrea_rna_seq | awk '{print $9}' | grep -v "scripts\|result_table_untrim.txt" > filenames
# ls -l /home/jconnell/niab/andrea_rna_seq | grep -v "scripts\|result_table_untrim.txt" | awk '{print $9 "-F"}' > newnames_f
# ls -l /home/jconnell/niab/andrea_rna_seq | grep -v "scripts\|result_table_untrim.txt" | awk '{print $9 "-R"}' > newnames_r
# sed -i 1d newnames_f
# sed -i 1d newnames_r
# #Get read counts 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/$name/*_1_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  f_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  f_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_f_pre
# 	done 
# done 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/$name/*_2_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  r_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  r_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_r_pre
# 	done 
# done  
# paste newnames_f f_read_counts_pre f_GC_counts_pre result_f_pre >> final_result_f_pre
# paste newnames_r r_read_counts_pre r_GC_counts_pre result_r_pre >> final_result_r_pre
# cat final_result_r_pre >> final_result_f_pre
# echo -e "Read name"'\t'"Read count"'\t'"GC content" > headder1
# cat /home/jconnell/niab/andrea_rna_seq/$name/*_2_*/summary.txt | cut -f2 | paste -s > headder2
# paste headder1	headder2 > headders
# sort -V final_result_f_pre > final_result
# cat final_result >> headders	
# mv headders	/home/jconnell/niab/andrea_rna_seq/result_table_untrim.txt
# rm newnames_f filenames newnames_r f_read_counts_pre r_read_counts_pre f_GC_counts_pre r_GC_counts_pre result_f_pre	result_r_pre final_result_f_pre	final_result_r_pre headder1 headder2 headders final_result	

####Trim_reads
 # IlluminaAdapters=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
 # v=$(ls -l /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData | awk '{print $9}')
 # for strain in $(echo $v); do
 #  for StrainPath in /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData; do
 #   ReadsF=$(ls $StrainPath/$strain/*.fq.gz | grep "_1.fq.gz")
 #   ReadsR=$(ls $StrainPath/$strain/*.fq.gz | grep "_2.fq.gz")
 #   F_outdir=/home/jconnell/niab/andrea_rna_seq/trimmed_reads/$strain/F
 #   R_outdir=/home/jconnell/niab/andrea_rna_seq/trimmed_reads/$strain/R
 #   mkdir -p $F_outdir $R_outdir
 #   ProgDir=/home/jconnell/niab/andrea_rna_seq/scripts
 #   sbatch $ProgDir/fastq_mcf.sh $ReadsF $ReadsR $IlluminaAdapters $F_outdir $R_outdir
 #  done   
 # done

####Post trimming fastQC
# v=$(ls -l /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData | awk '{print $9}')
# for f in $(echo $v); do
#   for RawData in $(ls /home/avadillo/projects/niab/andrea_rna_seq/trimmed_reads/${f}/*/*.fq.gz); do
#       OutDir=/home/avadillo/projects/niab/andrea_rna_seq/trimmed_reads_QC/${f}
#       mkdir -p $OutDir
#       ProgDir=/home/avadillo/projects/niab/andrea_rna_seq/scripts
#       sbatch $ProgDir/fastqc_ud.sh $RawData $OutDir
#   done
# done  

####Collect QC data past trim 
# ls -l /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC | awk '{print $9}' > filenames
# ls -l /home/jconnell/projects/niab/andrea_rna_seq/trimmed_reads | awk '{print $9 "-F"}' > newnames_f
# ls -l /home/jconnell/projects/niab/andrea_rna_seq/trimmed_reads | awk '{print $9 "-R"}' > newnames_r
# sed -i 1d newnames_f
# sed -i 1d newnames_r
# #Get read counts 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_1_trim_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  f_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  f_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_f_pre
# 	done 
# done 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_2_trim_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  r_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  r_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_r_pre
# 	done 
# done  
# paste newnames_f f_read_counts_pre f_GC_counts_pre result_f_pre >> final_result_f_pre
# paste newnames_r r_read_counts_pre r_GC_counts_pre result_r_pre >> final_result_r_pre
# cat final_result_r_pre >> final_result_f_pre
# echo -e "Read name"'\t'"Read count"'\t'"GC content" > headder1
# cat /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_2_*/summary.txt | cut -f2 | paste -s > headder2
# paste headder1	headder2 > headders
# sort -V final_result_f_pre > final_result
# cat final_result >> headders	
# mv headders	/home/jconnell/niab/andrea_rna_seq/result_table_trim_QC.txt
# rm newnames_f filenames newnames_r f_read_counts_pre r_read_counts_pre f_GC_counts_pre r_GC_counts_pre result_f_pre	result_r_pre final_result_f_pre	final_result_r_pre headder1 headder2 headders final_result	

####rRNA decontamination  
# a=$(ls -l /home/jconnell/niab/andrea_rna_seq/trimmed_reads | awk '{print $9}')
# for x in $a; do 
#     for RawData in /home/jconnell/niab/andrea_rna_seq/trimmed_reads/${x}; do
#         FileF=$RawData/F/*.fq.gz
#         FileR=$RawData/R/*.fq.gz
#         echo $FileF
#         echo $FileR
#         Ref=/home/jconnell/projects/niab/johnc/RNA_seq_data/ribokmers.fa.gz
#         outdir=/home/jconnell/niab/andrea_rna_seq/decontaminated_reads/${x}
#         mkdir -p $outdir
#         ProgDir=/home/jconnell/niab/andrea_rna_seq/scripts
#         sbatch $ProgDir/bbduk.sh $Ref $FileF $FileR $x $outdir
#     done
# done


##############RUN Salmon  *fixed*

#Pss9644 Done 
# for strain in RNA_1 RNA_9_2607 RNA_14 RNA_16 RNA_18 RNA_19; do 
# Transcriptome=/home/jconnell/niab/andrea_rna_seq/ref_genomes/GCF_023277945.1/ncbi_dataset/data/GCF_023277945.1/cds_from_genomic.fna
#     for RawData in /home/jconnell/niab/andrea_rna_seq/decontaminated_reads/$strain; do
#         FileF=$RawData/F/*cleaned.fq.gz
#         FileR=$RawData/R/*cleaned.fq.gz
#         echo $FileF
#         echo $FileR
#         OutDir=/home/jconnell/niab/andrea_rna_seq/salmon/Pss/$strain
#         mkdir -p $OutDir
#         ProgDir=/home/jconnell/andrea_rna_seq/scripts
#         sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
#     done
# done



#R2leaf 
for strain in RNA_10_2607 RNA_24 RNA_6 RNA_17 RNA_23 RNA_4; do 
#Transcriptome=/home/jconnell/andrea_rna_seq/ref_genomes/GCF_002905795.2/ncbi_dataset/data/GCF_002905795.2/cds_from_genomic.fna
Transcriptome=/home/jconnell/andrea_rna_seq/ref_genomes/R2Leaf_genebank_seq/ncbi_dataset/data/GCA_002905875.2/cds_from_genomic.fna
    for RawData in /home/jconnell/niab/andrea_rna_seq/decontaminated_reads/$strain; do
        index=/home/jconnell/niab/johnc/transcripts_index
        FileF=$RawData/F/*cleaned.fq.gz
        FileR=$RawData/R/*cleaned.fq.gz
        echo $FileF
        echo $FileR
        OutDir=/home/jconnell/niab/andrea_rna_seq/salmon/R2Leaf/$strain
        mkdir -p $OutDir
        ProgDir=/home/jconnell/andrea_rna_seq/scripts
        sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
done

#R15244 Done 
# for strain in RNA_2 RNA_3 RNA_5_2607 RNA_11 RNA_12 RNA_15; do 
# Transcriptome=/home/jconnell/niab/andrea_rna_seq/ref_genomes/GCF_002905685.2/ncbi_dataset/data/GCF_002905685.2/cds_from_genomic.fna
#     for RawData in /home/jconnell/niab/andrea_rna_seq/decontaminated_reads/$strain; do
#         index=/home/jconnell/niab/johnc/transcripts_index
#         FileF=$RawData/F/*cleaned.fq.gz
#         FileR=$RawData/R/*cleaned.fq.gz
#         echo $FileF
#         echo $FileR
#         OutDir=/home/jconnell/niab/andrea_rna_seq/salmon/R15244/$strain
#         mkdir -p $OutDir
#         ProgDir=/home/jconnell/andrea_rna_seq/scripts
#         sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir $index   
#     done
# done

#R15300 Done
# for strain in RNA_7 RNA_8_2607 RNA_13 RNA_20 RNA_21 RNA_22 RNA_21 RNA_20 ; do 
# Transcriptome=/home/jconnell/niab/andrea_rna_seq/ref_genomes/GCF_002905875.2/ncbi_dataset/data/GCF_002905875.2/cds_from_genomic.fna
#     for RawData in /home/jconnell/niab/andrea_rna_seq/decontaminated_reads/$strain; do
#         FileF=$RawData/F/*cleaned.fq.gz
#         FileR=$RawData/R/*cleaned.fq.gz
#         echo $FileF
#         echo $FileR
#         OutDir=/home/jconnell/projects/niab/andrea_rna_seq/salmon/R15300/$strain
#         mkdir -p $OutDir
#         ProgDir=/home/jconnell/andrea_rna_seq/scripts
#         sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
#     done
# done


####Create TX2GENE file 
cds=/home/jconnell/andrea_rna_seq/ref_genomes/R2Leaf_genebank_seq/ncbi_dataset/data/GCA_002905875.2/cds_from_genomic.fna
cat $cds | grep ">" | awk '{print $1}' | sed 's/>//g' > trans
cat $cds | grep ">" |  awk '{print $1}' | grep -Eo '[0-9]+$' | awk '{print "g"$0}' > gene_names
paste trans gene_names > tx2gene.txt
rm trans gene_names

####Blast effectors agaisnt genome and extract from DEG data 
#Create blast nucelotide database of CDS file 
makeblastdb -in cds_from_genomic.fna -parse_seqids -blastdb_version 5  -out R2Leaf_DB/R2Leaf -dbtype nucl 
#Run blast of type 3 effectors against database 
tblastn -query t3es.fasta -db R2Leaf_DB/R2Leaf -out t3es_R2leaf_hits_table.txt -evalue 1e-30 -outfmt 6

####Make a unique table of effectors for each CDS hit with blast %
#Create unique CDS list 
es t3es_R2leaf_hits_table.txt  | cut -f2 | sort | uniq -c | sort -nr | awk '{print $2}' | grep -v "Sub" > unique_cds_all
#Subset effector hits by cds taking unique + %match + match lenth
for x in $(cat unique_cds_all); do 
    cat t3es_R2leaf_hits_table.txt | grep $x | sort | uniq -c | awk '{print $2" "$4"% "$5"bp"}' >> "$x"_unique_hits.txt
done
#Begin to creat table 
es unique_cds_all | paste -s > table.txt
#Paste all data togetger in order of unique_cds_all
paste $(for x in $(cat unique_cds_all); do echo "$x"_unique_hits.txt; done) >> data 
#Combine all data
cat data >> table.txt

####Extract unique hits from l2FC table for plotting 
for x in $(cat unique_cds_all | grep -Eo '[0-9]+$' | awk '{print "g" $0}'); do
    cat R2Leaf_all_data.txt | grep -w $x >> effector_l2FC_results.txt
done 