# load data
  library(tidyverse)
  library(ggpubr)
  setwd("C:/Users/sfhar/Dropbox (Personal)/Personal Files/CATFISH/tcga")
  tcga_skcm <- read_tsv(file = "TCGA_SKCM_SNVs.txt")
  tcga_skcm

# process and visualize
  tcga_skcm %>%
    count(IID) %>%
    ggplot(aes(n))+
    geom_histogram()
  tcga_skcm2 <- tcga_skcm %>%
    mutate(snp = paste0(chr,"-",pos,".",ref,">",alt)) %>%
    #filter(alt_count/depth > 0.15) %>%
    select(IID,snp) %>%
    unique()
  tcga_skcm2_n <- count(tcga_skcm2, snp) %>%
    arrange(desc(n))

plotA <- tcga_skcm2_n  %>%
  mutate(count = case_when(n < 10 ~ as.character(n),
                           n >= 10 & n < 16 ~ "10-16",
                           n >= 16 ~ "16-493")) %>%
  count(count) %>%
  mutate(count = factor(count, levels = c(as.character(1:9), "10-16", "16-493"))) %>%
  ggplot(aes(x = count, y = n)) +
  geom_col()+
  labs(x = "Number of melanomas sharing somatic mutation (n=463)",
       y = "Number of SNVs") +
  geom_text(aes(label = n), vjust = -0.5) + 
  theme_classic()
plotA

perm <- 1
output <- data.frame()
set.seed(12345)

for(perm in 1:100){
  tcga_skcm <- filter(tcga_snvs, IID %in% sample(SKCMs, 16))
  tcga_skcm2 <- tcga_skcm %>%
    mutate(snp = paste0(chr,"-",pos,".",ref,">",alt)) %>%
    #filter(alt_count/depth > 0.15) %>%
    select(IID,snp) %>%
    unique()
  tcga_skcm2_n <- count(tcga_skcm2, snp) %>%
    arrange(desc(n))
  result <- tcga_skcm2_n  %>%
    mutate(count = factor(n, levels = 1:16)) %>%
    count(count, .drop = FALSE) %>%
    mutate(rep = perm)
  output <- rbind(output, result)
}
output_mean <- output %>%
  group_by(count) %>%
  summarize(mean = mean(n), sd = sd(n))

output 
output_mean 
plotB <- output %>%
  ggplot(aes(x = count, y = n))+
  geom_boxplot(outlier.shape = NA) +
  ylim(0,16000)+
  # geom_text(data = output_mean, aes(x = count, y = 15000, label = round(mean, 2)))+
  geom_text(data = output_mean, aes(x = count, y = mean, label = round(mean, 2)), vjust = -0.5) +
  geom_text(x = 1, y = 15000, label = "mean =")+
  labs(x = "Number of melanomas sharing somatic mutation (n=16 randomly sampled)",
       y = "Number of SNVs"#,
       #title = "100 permutations, randomly sampling 16 melanomas from TCGA"
       ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.05))) +
  theme_classic()
plotB

plotAB <- ggarrange(plotA, plotB, nrow = 1, ncol = 2, labels = c("A","B"))

ggsave("TCGA_SKCM_fig.png",
       plot = plotAB,
       width = 14,
       height = 5,
       units = "in")
