library(tidyverse)
library(scales)

sample_list <- read.delim('/home/chase/code/cmatkhan/genomics/assignment4/sample_list.txt', header = FALSE, sep='\n')
sample_list <- sample_list$V1

raw_fltr_lib_size <- read.delim('/home/chase/code/cmatkhan/genomics/assignment4/filtered_lib_size.txt' , header = FALSE, sep = '\n')
raw_fltr_lib_size <- raw_fltr_lib_size$V1
names(raw_fltr_lib_size) <- sample_list
raw_fltr_df = data.frame(raw_fltr_lib_size) / 1000000

norm_lib_size <- read.delim('/home/chase/code/cmatkhan/genomics/assignment4/norm_lib_size.txt', header = FALSE, sep='\n')
norm_lib_size <- norm_lib_size$V1
names(norm_lib_size) <- sample_list
norm_df = data.frame(norm_lib_size) / 1000000
norm_df$sample <- factor(sample_list)

rab30 <- read_csv('/home/chase/code/cmatkhan/genomics/assignment4/RAB30_exp.csv')
before_after = c(rab30[rab30$X1=='Before_mean', ]$RAB30, rab30[rab30$X1=='After_mean', ]$RAB30)

raw <- ggplot(raw_fltr_df, aes(x=rownames(raw_fltr_df), y=raw_fltr_df$raw_fltr_lib_size)) +
  geom_bar(stat='identity') +
  ggtitle('Filtered Raw Counts') +
  xlab("sample") +
  ylab("Library Size (mil)") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("library_size.png" width = 10, height = 10)

norm <- ggplot(norm_df, aes(x=sample, y=norm_lib_size)) +
  geom_bar(stat='identity') +
  ggtitle('Normalized Raw Counts') +
  xlab("sample") +
  ylab("Library Size (mil)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("library_size_normalzied.png" width = 10, height = 10)

rab30_boxplot <- barplot(x, main = "Normalized Count of RAB30", xlab = 'Sample', ylab = 'Count', ylim = c(0,200))


print(raw)
ggsave("library_size.png", width = 10, height = 10)
print(norm)
ggsave("library_size_normalzied.png", width = 10, height = 10)
print(rab30_boxplot)
ggsave("mean_expression.png", width = 10, height = 10)
