#!/usr/bin/env python 

#Script function: filter pyani results to give higghest matches for each strain 

#Import packages 
import argparse 

#Parse files 
ap=argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help='input pyani file') 
args=ap.parse_args()

#Open files 
infile = open(args.i)

strain_dict = {}

for x in infile:
    x = x.replace("\n","")
    x = x.split("\t")
    strain = x[0].split("_")[0]
    if strain not in strain_dict:
    	strain_dict	[strain] = x[1:]
    else:
    	strain_dict	[strain].append[x[1:]]

for x in strain_dict:
    print (x, ":", strain_dict[x])
