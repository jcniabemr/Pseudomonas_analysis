#!/usr/bin/env Rscript 

####Script to align data for Pseudomonas gene table
####Only requied argument is <infile -i> a combined gene table made during pseudomonas_gene_database.sh

####Import required libs 
library(optparse)
####Set option list 
option_list = list(
    make_option(c("-i","--i"), help="infile"),
    make_option(c("-o","--o"), help="output file", default="gene_table_matched")
    )
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$i)){
	print_help(opt_parser)
	stop("Please provide gene table",, call.=FALSE)
}
####Read in data
data <- read.table(opt$i, header = T, sep = '\t', fill = TRUE)
####Creat custom match function and apply to cols of data matching all data to col1
reformed <- apply(data, 2, FUN=function(x) {
	x[match(data[,1] ,x)]
	})
####Switch NA > 0 and transform 
reformed[is.na(reformed)] <-0
final <- t(reformed)
####Set name and write table 
filename <- paste(opt$o,".txt", sep = "")
write.table(final, file = filename, col.names = FALSE, sep = '\t', quote = FALSE)



