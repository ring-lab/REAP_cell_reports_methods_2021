library(readr)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(forcats)
library(FactoMineR)

##### Input Variables ####

# Set the directory
setwd("change_me")

# Set the input file. Input file is a csv of read counts for yeast clones.
# Each row is a unique yeast clone (besides the first row with column names). 
# The first column contains protein names, the second column contains the 
# barcode sequence associated with the yeast clone, the third column contains 
# the counts from the pre-selection library, and the rest of the columns 
# contain counts from post-selection libraries.
input_bc_file <- "change_me.csv"

# Set a name for the output folder
folder_name <- "change_me"

##### Import ####

# Read in the CSV file
bc_counts <- read_csv(input_bc_file, col_names = TRUE)
bc_counts_head <- read_csv(input_bc_file, n_max = 1, col_names = FALSE)
colnames(bc_counts)[1] <- "Protein"

# Remove the column containing barcode sequences
bc_counts <- bc_counts[,-2]
bc_counts_head <- bc_counts_head[,-2]

# Create vector of samples
samples <- unlist(bc_counts_head[-1])

# Create and move to output folder
dir.create(folder_name)
setwd(folder_name)

##### Collapse Barcodes ####

# Get a tibble of all the unique identified genes
genes <- unique(bc_counts[,1])

collapsed <- tibble()
barcode_num <- tibble()
# Looping over all the unique genes
for (i in 1:nrow(genes)){
  # Take all the rows for the same protein
  # Remove the column with the protein name
  t <- bc_counts %>% subset(bc_counts$Protein == genes$Protein[i]) %>% 
    select(-Protein)
  # Add the sum of all the columns (all the different 
  # barcodes for a given protein) to the "collapsed" array
  collapsed <- collapsed %>% bind_rows(colSums(t))
}

##### Calculating Total Enrichment ####

# function that iterates over the groups performing an exact test for difference between
# each group and group 1 (pre-library) and adding the logFC between each to a matrix
calcexp <- function(n){
  exp <- n$genes
  groups <- as.vector(unique(n$samples[,"group"]))
  #creates a non-repetitive list of numeric groups present in analysis
  for (i in 2:length(unique(samples))) {
    et <- exactTest(n, pair = c(samples[1],unique(samples)[i])) 
    #conducts a pairwise test comparing iterative samples to the first sample
    exp <- cbind(exp, et$table[,'logFC']) 
    #adds logFC of each comparison to a matrix
    groupname <- groups[i]
    #creates group name for the matrix columns
    colnames(exp)[i] <- groupname
  }
  return(exp)
}

# calculate fold change
expmatrix <- collapsed %>% DGEList(genes = genes, group = samples) %>%
  calcNormFactors() %>%
  estimateDisp() %>%
  calcexp()

# creates fold enrichment matrix
FE.noneg<- expmatrix %>%
  tibble() %>% 
  select(-Protein) %>%
  mutate_all(function(x) if_else(x < 0,0,x)) %>%
  bind_cols(genes) %>%
  relocate(Protein)

##### Calculating Clonal Enrichment #### 

# Minimum log fold change to be counted as an enriched barcode
enriched_barcode_fc <- 2

# Separate uncollapsed gene names and counts
bc_counts2_genes <- bc_counts %>% select(Protein)
bc_counts2 <- bc_counts %>% select(-Protein)

# calculate fold change for each barcode independently
expmatrix2 <- bc_counts2 %>% DGEList(genes = bc_counts2_genes, group = samples) %>%
  calcNormFactors() %>%
  estimateDisp() %>%
  calcexp()

# create fold enrichment matrix
FE.noneg.b <- expmatrix2 %>%
  tibble() %>% 
  select(-Protein) %>%
  mutate_all(function(x) if_else(x < 0,0,x))

# converts logical matrix of barcode enrichemnt where enrichment = 1
enrich.logical <- FE.noneg.b %>%
  mutate_all(function(x) if_else(x > enriched_barcode_fc,1,0)) %>%
  bind_cols(bc_counts2_genes) %>%
  relocate(Protein)

# create clonal enrichment matrix
frequency.b <- tibble()
for (i in 1:nrow(genes)){
  # Looping over all the unique genes
  t <- enrich.logical %>% filter(enrich.logical$Protein == genes$Protein[i]) %>%
    select(-Protein)
  # Take all the rows for the same protein. Remove the column with the protein name
  sums <- colSums(t)
  freq <- sums/barcode_num$Number[i]
  frequency.b <- frequency.b %>% bind_rows(freq)
}
frequency.b <- bind_cols(genes, frequency.b)

##### Calculating REAP Score ####

# create tibble containing barcode number correction factor
bc_num <- barcode_num %>%
  select(-Protein) %>%
  mutate_all(function(x) if_else(x <= 5, log(x + 0.5)/1.705, 1))

# create tibble containing protein frequency correction factor
num_frequency_l <- tibble(collapsed[,1]/sum(collapsed[,1])) %>%
  log10() %>%
  mutate_all(function(x) if_else(x < -6,-6,x)) %>%
  mutate_all(function(x) if_else(x == 0,-6,x)) %>%
  mutate_all(function(x) if_else(x <= -4, log(x + 7.1)/1.16, 1))

# create tibble of scores
score <- select(frequency.b,-Protein)^2*
  select(FE.noneg,-Protein)*
  unlist(num_frequency_l)*
  unlist(bc_num)
score <- bind_cols(genes, score)

write.csv(score, file = "score.csv", row.names = FALSE)

