library(readr)
library(tidyr)
library(dplyr)

# scores is a csv file with REAP scores. The first column contains the protein
# name. The rest of the columns are individual patients and cells within a column
# are the REAP scores for the given patient for the respective protein. The first row
# contains patient sample names.
scores <- read_csv("change_me.csv")
scores <- as.data.frame(scores)
rownames(scores) <- scores[,1]
scores <- scores[,-1]

# clinical is a csv file. The first column contains the various clinical conditions
# being assessed. The rest of the columns are individual patients and cells within
# a column are a binary indication of whether a given patient has a given clinical
# condition. The first row contains patient sample names.
clinical <- read_csv("change_me.csv")
conditions <- clinical[,1]
clinical <- clinical[,-1]
samples <- colnames(clinical)

# filter clinical so that only conditions that are present in at least one patient
# are included
conditions.present <- rowSums(clinical) > 0
clinical <- clinical[conditions.present,]
conditions.present <- unlist(conditions)[conditions.present]

# loop through and test for significant differences in REAP score for every
# protein x condition pair
full.pvalue <- array()
for(i in 1:length(conditions.present)){
  # generates list of TRUE/NA values for whether a condition is present in the row
  conditions.tf <- clinical[i,] > 0
  
  # generates vector containing positive sample names
  positive <- samples[conditions.tf]
  
  # generates vector containing negative sample names
  negative <- samples[!conditions.tf]
  
  # filters scores into separate variables based on whether they are in patients
  # with or without a given clinical condition
  scores.pos <- as.data.frame(scores[,positive])
  scores.neg <- as.data.frame(scores[,negative])
  
  # perform wilcox test for each protein
  pvalue.select <- tibble()
  for(i in 1:nrow(scores)){
    mw.test <- wilcox.test(as.numeric(scores.pos[i,]), as.numeric(scores.neg[i,]))
    pvalue.select[i,1] <- mw.test$p.value
  }
  full.pvalue <- cbind(full.pvalue, pvalue.select)
}
full.pvalue <- full.pvalue[,-1]
colnames(full.pvalue) <- conditions.present
rownames(full.pvalue) <- rownames(scores)

# removes any pathways that have less than 3 patients positive (score >= 3) for the pathway
scores.cut <- sapply(scores, function(x) ifelse(x >= 3, x <- 1, x <- 0))
rownames(scores.cut) <- rownames(scores)
scores.sum <- as.data.frame(rowSums(scores.cut))
scores.sum <- filter(scores.sum, scores.sum[,1] >= 3)
sub.pvalue <- subset(full.pvalue, rownames(full.pvalue) %in% rownames(scores.sum))

# adjust the p values using FDR
sub.pvalue.adj <- sapply(sub.pvalue, function(x) p.adjust(x, method = "fdr"))
rownames(sub.pvalue.adj) <- rownames(sub.pvalue)

write.csv(sub.pvalue.adj, file = "output.csv")
