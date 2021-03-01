library(parallel)
library(data.table)
library(tidyverse)
#v1.2.13 available on github
#https://github.com/samuelbfernandes/simplePHENOTYPES/tree/0eabefa0bc79c7a0abd2884b11498166f29e4daf
library(simplePHENOTYPES) 
source("./run_simulation_P.R")
source("./run_gemma_P.R")
MAF <-  c(0.05, 0.4)
h2 <- matrix(c(0.3, 0.3, 0.8, 0.3, 0.8, 0.8), 3)
# zipF <- list.files(path = "./data", pattern = "*.zip", full.names = TRUE)
# plyr::ldply(.data = zipF, .fun = unzip, exdir = "./data")
# plyr::ldply(.data = "./data/maize_kinship.zip", .fun = unzip, exdir = "./data")
# plyr::ldply(.data = "./data/soy_kinship.zip", .fun = unzip, exdir = "./data")
data <- c("ames_500_mac5_ld09_numeric.txt",
          "ames_1000_mac5_ld09_numeric.txt",
          "ames_2815_mac5_ld09_numeric.txt",
          "soy_500_mac5_ld09_numeric.txt",
          "soy_1000_mac5_ld09_numeric.txt",
          "soy_2815_mac5_ld09_numeric.txt")
#seed = round(runif(12, 0, 2000000))
seed <- c(1014956, 613537, 853815, 1386204, 170272, 317557, 549061, 544610, 1231659, 859343, 1303311, 1135476)
grid <- expand.grid(maf = MAF, data = data, stringsAsFactors = F)
grid <- grid[order(grid$data, grid$maf),]
settings <-  data.frame(i = 1:12, seed = seed, grid)
settings <- split(settings, settings$i)
out <- "./simulation/P/"
for (x in settings) {
  cat(paste("Setting", x$i, "\n"))
  run_simulation_P(geno_file = x$data, list = x, out = out, h2 = h2, rep = 10)
  path  <-
    paste0(out, gsub("_mac5.*", "", x$data), "_MAF", x$maf, "/")
  cat(paste("GEMMA: Setting", x$i, "\n"))
  run_gemma_P(path = path, pop = x$data, n_cores = detectCores())
}

# correlation

setwd("./simulation/P")
x <- dir()
correl <- tibble(
  specie = character(3600),
  h2 = character(3600),
  sample_size = character(3600),
  ld = character(3600),
  maf = character(3600),
  type = character(3600),
  cor = numeric(3600))
a <- 1
for (i in x){ 
  y <-  dir(paste0("./",i))[grepl("Simulated_Data", dir(paste0("./",i)))]
  for (ii in y){
    temp <- fread(paste0("./", i, "/", ii ), data.table = F)
    correl[a, "h2"] <- gsub(".*Herit_|.fam","", ii)
    correl[a, "cor"] <- cor(temp[,6], temp[,7])
    correl[a, "type"] <- "P"
    correl[a, "specie"] <- gsub("_.*","", i)
    correl[a, "sample_size"] <- gsub(".*_", "", gsub("_MAF.*","", i))
    correl[a, "maf"] <- gsub(".*MAF|_.*","", i)
    a <- a + 1
    print(a)
  }
}

fwrite(correl, "correlation_P.txt", sep = "\t", quote = F, row.names = F)
