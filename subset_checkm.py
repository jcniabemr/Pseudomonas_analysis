#!/usr/bin/env python

#Script to filter checkM data after subsetting
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help='Subset checkM file')
ap.add_argument('-o',required=True,type=str,help='outfile')
args = ap.parse_args()
infile = open(args.i)
outfile = open(args.o, 'w')

for x in infile:
	x = x.replace("\n","")
	x = x.split("\t")
	if float(x[1]) >= 95:
		if float(x[2]) <= 5:
			if int(x[6]) >= 40000:
				data = x[0:]
				out = "\t".join(data)
				outfile.write(out + "\n")
infile.close()
outfile.close()
