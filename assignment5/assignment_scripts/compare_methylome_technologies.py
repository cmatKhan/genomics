"""
Write a script called ​compare_methylome_technologies.py​ that

1.Creates the following scatter plots

    a.MeDIP-seq RPKM vs. MRE-seq RPKM.
    Save the plot as ​<MeDIP-seq RPKM bedbasename>​_vs_​<MRE-seq RPKM bed basename>​.png​

    b.MeDIP-seq RPKM vs. WGBS average DNA methylation level.
    Save the plot as<MeDIP-seq RPKM bed basename>​_vs_​<WGBS methylation level bedbasename>​.png​

    c.MRE-seq RPKM vs. WGBS average DNA methylation level.
    Save the plot as<MRE-seq RPKM bed basename>​_vs_​<WGBS methylation level bedbasename>​.png​

2.Calculate the correlation for each comparison. Print the correlation to stdout. Hint: use scipy.stats​.

Run ​compare_methylome_technologies.py​ on ​MeDIP_CGI_RPKM.bed​,​​MRE_CGI_RPKM.bed​,and ​WGBS_CGI_methylation.bed​.
Save the output figures as

MeDIP_CGI_RPKM_vs_MRE_CGI_RPKM.png​

MeDIP_CGI_RPKM_vs_WGBS_CGI_methylation.png​

MRE_CGI_RPKM_vs_WGBS_CGI_methylation.png​

Paste the command to generate the figures,the correlations, and justify which correlation statistic you used in your README.

"""