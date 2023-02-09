#!/usr/bin/env R

#Script to quickly plot dendogram from pyani tsv output 

library(optparse)

option_list = list(
  make_option(c("-i","--i"), help="infile name"), 
  make_option(c("-o","--o"), help="output file name", default="dendogram.png"),
  make_option(c("-t","--t"), help="title", default="Dendogram"),
  make_option(c("-x","--x"), help="x lab", default="Genome Assembly"),
  make_option(c("-y","--y"), help="y lab", default="Height"),
  make_option(c("-m","--m"), help="method", default="single")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$i)){
  print_help(opt_parser)
  stop("please provide tsv file after -i flag", call.=FALSE)
}
data <- read.table(opt$i, row.names = 1, header = T, sep = '\t')
dist_data <- as.dist(data)
cluster <- hclust(dist_data, method=opt$m) 
method <- paste("(","method","-",opt$m,")", sep ="") 
filename <- paste(opt$o,".png", sep="")
jpeg(file=filename)
plot(cluster, main = opt$t, ylab=opt$y, xlab = opt$x, sub = method)
dev.off


