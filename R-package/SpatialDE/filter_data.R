library(dplyr)
library(readr)
mouseob_small <- read_csv("../../Analysis/MouseOB/data/Rep11_MOB_0.csv")
mouseob_small <- mouseob_small %>%
  semi_join(locations)
rn <- mouseob_small$X1
mouseob_small <- as.matrix(mouseob_small[,-1])
rownames(mouseob_small) <- rn
mouseob_small <- mouseob_small[,sample(seq_len(ncol(mouseob_small)), 100)]
