#!/usr/bin/env python 

####Script to extract gene names and CDS number from CDS file downloaded from NCBI  

####Import libs
import argparse

####Parse files 
ap = argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help="CDS file")
parse = ap.parse_args()

####Open files
infile = open(parse.i)

cds_list = []
cds_dict = {}

for x in infile:
    if ">" in x:
    	cds_list.append(x)
    else:
    	pass 

for x in cds_list:
	x = x.split(" ")
	list1 = x[1] + '\t' + x[2]
	gene = list1.split('\t')[0]
	gene = gene.replace("[","").replace("]","").replace("locus_tag=","")
	cds_num = x[0].split("_")[-1]
	if cds_num not in cds_dict:
		cds_dict[cds_num] = gene
	else:
		cds_dict[cds_num].append(gene)

for x in cds_dict:
    print (x, ":", cds_dict[x])
