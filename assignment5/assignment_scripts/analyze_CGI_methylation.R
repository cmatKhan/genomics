
#  This script takes the average CGI methylation bed file, plots the distribution of average CGI methylation levels,
#  and save the plot as <averageCGI methylation bed basename>._distribution.png

library(tidyverse)
library(calecopal)

avg_cgi_metyl = read.table('/home/chase/code/cmatkhan/genomics/assignment5/WGBS_CGI_methylation.bed', sep='\t')

hist(avg_cgi_metyl$V5, breaks = 50, 
     probability = TRUE,
     col = cal_palette(name = "bigsur", n=1),
     border = "black",
     xlab = 'CGI Methylation Level',
     ylab = 'Density',
     main = 'Average CGI methylation')
lines(density(avg_cgi_metyl$V5), col = cal_palette("kelp1"), lwd=3)

# > df = read_csv('methylation_coverage.csv')
# Parsed with column specification:
#         cols(
#                 methyl_level = col_double(),
#                 coverage = col_double()
#         )
# > head(df)
# 
# > hist(df$methyl_level)
# > hist(df$coverage)
# > hist(df$coverage, breaks = 100)
# > hist(df$coverage, breaks = 100, probability = TRUE)
# > hist(df$coverage, breaks = 1000, probability = TRUE)
# > hist(df$coverage, breaks = 1000, probability = TRUE, xlim = c(0,200))
# > hist(df$coverage, breaks = 3, probability = TRUE, xlim = c(0,200))
# > hist(df$coverage, breaks = 30, probability = TRUE, xlim = c(0,200))
# > hist(df$coverage, breaks = 10000, probability = TRUE, xlim = c(0,200))