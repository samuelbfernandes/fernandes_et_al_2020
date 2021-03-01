library(parallel)
library(data.table)
library(tidyverse)
#v1.2.13 available on github
#https://github.com/samuelbfernandes/simplePHENOTYPES/tree/0eabefa0bc79c7a0abd2884b11498166f29e4daf
library(simplePHENOTYPES) 
source("./run_simulation_ST.R")
source("./run_gemma_ST.R")
MAF <-  c(0.05, 0.4)
h2 <- matrix(c(0.3, 0.3, 0.8, 0.3, 0.8, 0.8), 3)
custom_a <- list(
  list(trait_1 = c(0.1),
       trait_2 = c(0.1))
)
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
#seed <- round(runif(12, 0, 2000000))
seed <- c(531017,  744248, 1145707, 1816416,  403364, 1796779, 1889351, 1321596, 1258228,  123573,  411949,  353114)
grid <- expand.grid(maf = MAF, effect = custom_a, data = data, stringsAsFactors = F)
grid <- grid[order(grid$data,  grid$maf),]
settings <-  data.frame(i = 1:12, seed = seed, grid)
settings <- split(settings, settings$i)
out <- "./simulation/ST/"
for (x in settings) {
  pop <- x$data
  run_simulation_ST(geno_file = pop, list = x, out = out, h2 = h2, rep = 100)
  path  <-
    paste0(out, gsub("_mac5.*", "", pop),"_MAF", x$maf,"/")
  cat(paste("GEMMA: Setting", x$i, "\n"))
  run_gemma_ST(path = path, pop = pop, n_cores = detectCores())
}

#-------------
x <- dir("./simulation/ST/")
a <- 1
correl <- vector("list", 12)
for (i in x){ 
  temp <- fread(paste0("./simulation/ST/", i, "/correlation.txt"), data.table = F)
  colnames(temp) <- c("0.3_0.3", "0.3_0.8", "0.8_0.8")
  correl[[a]] <-  tibble(
    specie = gsub("_.*","", i),
    sample_size = gsub(".*_", "", gsub("_MAF.*","", i)),
    ld = "",
    maf = gsub(".*MAF|_.*","", i),
    type = "ST", 
    pivot_longer(temp, names_to = "h2", values_to = "cor", cols = `0.3_0.3`:`0.8_0.8`) )
  print(a)
  a <- a + 1
}
correl <- rbindlist(correl) 
correl <- correl %>% select(specie,	h2,	sample_size,	ld,	maf,	type,	cor)
fwrite(correl, "correlation_ST.txt", sep = "\t", quote = F, row.names = F)
