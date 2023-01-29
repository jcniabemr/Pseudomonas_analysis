#!/usr/bin/env python

#Script to filter checkM data after subsetting
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help='Subset checkM file')
ap.add_argument('-po',required=True,type=str,help='outfile passing filter')
ap.add_argument('-fo',required=True,type=str,help='data failing filter')
args = ap.parse_args()

infile = open(args.i, 'r')
outfile = open(args.po, 'w')
outfile2 = open(args.fo, 'w')


fail = []

for x in infile:
	x = x.replace("\n","")
	x = x.split(",")
	if float(x[22]) >= 95:
		if float(x[24]) <= 5:
			if int(x[44]) >= 40000:
				data = x[0:58]
				out = "\t".join(data)
				outfile.write(out + "\n")
	if float(x[22]) < 95:
		fail.append(x[0] + x[22])
	elif float(x[24]) > 5:
		fail.append(x[0] +  "Contamination > 5"+ x[24])
	elif int(x[44]) < 40000:
		fail.append(x[0] + "N50 <" + x[44]) 

out = "\n".join(fail)
outfile2.write(out)

infile.close()
outfile.close()
outfile2.close()


