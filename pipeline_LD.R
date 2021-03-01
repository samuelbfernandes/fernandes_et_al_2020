library(parallel)
library(data.table)
library(tidyverse)
#v1.2.13 available on github
#https://github.com/samuelbfernandes/simplePHENOTYPES/tree/0eabefa0bc79c7a0abd2884b11498166f29e4daf
library(simplePHENOTYPES) 
source("./run_simulation_LD.R")
source("./run_gemma_LD.R")
MAF <-  c(0.05, 0.4)
h2 <- matrix(c(0.3, 0.3, 0.8, 0.3, 0.8, 0.8), 3)
ld <- c(0.99)
ld_type <- c("indirect","direct")
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
#seed <- round(runif(24, 0, 2000000))
seed <- c(603993,1584586,663709,1264425,227020,688405,1391618,1089958,1988011,1175267,967792,51371,1724582,1035043, 617519,1166864,1406062,455126,644204,816707,770381,802332,1837742,1145493)

grid <- expand.grid(ld = ld, maf = MAF, ld_type = ld_type, data = data, stringsAsFactors = F)
grid <- grid[order(grid$data, rev(grid$ld_type), grid$ld, grid$maf),]
settings <-  data.frame(i = 1:24, seed = seed, grid)
settings <- split(settings[order(settings$i, decreasing = T),], settings$i)
out <- "./simulation/LD/"
for (x in settings) {
  cat(paste("Setting", x$i, "\n"))
  run_simulation(geno_file = x$data, list = x, out = out, h2 = h2, ld_type = x$ld_type, rep = 3000)
  if(x$ld_type == "direct"){
    temp <-  fread(paste0(out, gsub("../../|_mac.*","", x$data),
                          "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,
                          "/LD_summary_Additive.txt"), data.table = F)
    temp2 <- temp[order(abs(temp$Actual_LD), decreasing = T),][1:100,]
    files <- dir(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type))
    files_remove <- unlist(lapply(paste0("Simulated_Data__Rep", setdiff(temp$rep, temp2$rep), "_Herit"), function(x)files[grepl(x, files)]))
    unlink(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/",files_remove), recursive = T, force = T)
  } else {
    temp <-  fread(paste0(out, gsub("../../|_mac.*","", x$data),
                          "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,
                          "/LD_summary_additive.txt"), data.table = F)
    temp2 <- temp[order(abs(temp$LD_between_QTNs), decreasing = T),][1:100,]
    files <- dir(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type))
    files_remove <- unlist(lapply(paste0("Simulated_Data__Rep", setdiff(temp$rep, temp2$rep), "_Herit"), function(x)files[grepl(x, files)]))
    unlink(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/",files_remove), recursive = T, force = T)
  }
  
  path  <-
    paste0(out, gsub("_mac5.*", "", x$data), "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/")
  cat(paste("GEMMA: Setting", x$i, "\n"))
  run_gemma_LD(path = path, pop = x$data, n_cores = detectCores())
}

new <- list(settings[[15]], settings[[13]])
new[[1]]$i <- 1
new[[2]]$i <- 2
#running additional settings 
out <- "./simulation/additional/"
x <- new[[1]]
run_simulation(geno_file = x$data, rep = 100,  list = x, out = out, h2 = h2, ld_type = x$ld_type)
temp <-  fread(paste0(out, gsub("../../|_mac.*","", x$data),
                      "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,
                      "/LD_summary_additive.txt"), data.table = F)
temp2 <- temp[order(abs(temp$LD_between_QTNs), decreasing = T),][1,]
files <- dir(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type))
files_remove <- unlist(lapply(paste0("Simulated_Data__Rep", setdiff(temp$rep, temp2$rep), "_Herit"), function(x)files[grepl(x, files)]))
unlink(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/",files_remove), recursive = T, force = T)
path  <-
  paste0(out, gsub("_mac5.*", "", x$data), "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/")
run_gemma(path = path, pop = x$data, n_cores = detectCores())

x <- new[[2]]
run_simulation(geno_file = x$data, rep = 100,  list = x, out = out, h2 = h2, ld_type = x$ld_type)
temp <-  fread(paste0(out, gsub("../../|_mac.*","", x$data),
                      "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,
                      "/LD_summary_Additive.txt"), data.table = F)
temp2 <- temp[order(abs(temp$Actual_LD), decreasing = T),][1:5,]
temp2 <- temp2[-3,]
files <- dir(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type))
files_remove <- unlist(lapply(paste0("Simulated_Data__Rep", setdiff(temp$rep, temp2$rep), "_Herit"), function(x)files[grepl(x, files)]))
unlink(paste0(out, gsub("../../|_mac.*","", x$data),"_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/",files_remove), recursive = T, force = T)

path  <-
  paste0(out, gsub("_mac5.*", "", x$data), "_ld_", x$ld, "_MAF", x$maf, "_", x$ld_type,"/")
run_gemma(path = path, pop = x$data, n_cores = detectCores())

# correlation
library(tidyverse)
setwd("./simulation/LD")
x <- dir()
correl <- tibble(
  specie = character(14400),
  h2 = character(14400),
  sample_size = character(14400),
  ld = character(14400),
  maf = character(14400),
  type = character(14400),
  cor = numeric(14400))
a <- 1
for (i in x){ 
  y <-  dir(paste0("./",i))[grepl("Simulated_Data", dir(paste0("./",i)))]
  for (ii in y){
    temp <- fread(paste0("./", i, "/", ii ), data.table = F)
    correl[a, "h2"] <- gsub(".*Herit_|.fam","", ii)
    correl[a, "cor"] <- cor(temp[,6], temp[,7])
    correl[a, "type"] <- gsub(".*_","", i)
    correl[a, "specie"] <- gsub("_.*","", i)
    correl[a, "sample_size"] <- gsub(".*_", "", gsub("_ld.*","", i))
    correl[a, "ld"] <- gsub(".*ld_|_MAF.*","", i)
    correl[a, "maf"] <- gsub(".*MAF|_.*","", i)
    a <- a + 1
  }
}

fwrite(correl, "correlation_LD.txt", sep = "\t", quote = F, row.names = F)
