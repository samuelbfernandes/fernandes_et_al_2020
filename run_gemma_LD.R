source("./gemma.R")
run_gemma_LD <- function(path = getwd(), pop, n_cores){
  files <- dir(path)
  files_fam <- files[grepl("Simulated_Data_", files)]
  bfiles <- files[grepl(".bed|.bim", files)]
  dir.create(paste0(path, "gemma"))
  n_files <- length(files_fam)
  files_fam <- split(files_fam, 1:n_files)
  invisible(mclapply(files_fam, function(x) {
    genoname <- gsub(".fam", "", x)
    prefix <- round(runif(1, min = 1, max = 10000))
    file.rename(paste0(path, x), paste0(path, "gemma/", x))
    file.copy(paste0(path, bfiles[1]), paste0(path, "B", prefix, bfiles[1]))
    file.copy(paste0(path, bfiles[2]), paste0(path, "B", prefix, bfiles[2]))
    file.rename(paste0(path, "B", prefix, bfiles[grepl(".bed", bfiles)]), paste0(path,"gemma/", genoname, ".bed"))
    file.rename(paste0(path, "B", prefix, bfiles[grepl(".bim", bfiles)]), paste0(path,"gemma/", genoname, ".bim"))
    gemma(
      geno = genoname,
      d = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenD.txt"),
      u = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenU.txt"),
      trait_col = c(1, 2),
      miss = 0.999,
      maf = 0.0001,
      r2 = 0.9999999,
      out_name = gsub("Simulated_Data__", "MT", genoname),
      path_out = paste0(path,"gemma" )
    )
    gemma(
      geno = genoname,
      d = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenD.txt"),
      u = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenU.txt"),
      trait_col = 1,
      miss = 0.001,
      maf = 0.001,
      r2 = 0.999999,
      out_name = gsub("Simulated_Data__", "T1_", genoname),
      path_out = paste0(path,"gemma" )
    )
    gemma(
      geno = genoname,
      d = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenD.txt"),
      u = paste0("../../../../data/kin_", gsub("_|_mac5.*", "", pop), ".eigenU.txt"),
      trait_col = 2,
      miss = 0.001,
      maf = 0.001,
      r2 = 0.999999,
      out_name = gsub("Simulated_Data__", "T2_", genoname),
      path_out = paste0(path,"gemma" )
    )
    unlink(c(
      paste0(path, "gemma/", genoname, ".bed"),
      paste0(path, "gemma/", genoname, ".bim"),
      paste0(path, "gemma/", genoname, ".fam")
    ),
    recursive = T,
    force = T)
  },
  mc.preschedule = FALSE,
  mc.cores = n_cores))
  unlink(paste0(path, files[grepl(".bed|.bim|.fam", files)]),
         recursive = T,
         force = T)
}
