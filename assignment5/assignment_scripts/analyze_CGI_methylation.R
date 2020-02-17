
#  This script takes the average CGI methylation bed file, plots the distribution of average CGI methylation levels,
#  and save the plot as <averageCGI methylation bed basename>._distribution.png

library(tidyverse)
library(calecopal)

avg_cgi_metyl = read.table('/home/chase/code/cmatkhan/genomics/assignment5/WGBS_CGI_methylation.bed', sep='\t')

by_cgi = aggregate(avg_cgi_metyl$V5, by=list(Category=avg_cgi_metyl$V4), FUN=sum)

hist(by_cgi$x, breaks = 50, 
     probability = TRUE,
     col = cal_palette(name = "bigsur", n=1),
     border = "black",
     xlab = 'CGI Methylation Level',
     ylab = 'Density',
     main = 'Average CGI methylation')
lines(density(by_cgi$x), col = cal_palette("kelp1"), lwd=3)
