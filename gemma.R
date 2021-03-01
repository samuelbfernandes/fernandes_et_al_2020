gemma <- function(gemma_path = "PATH/TO/GEMMA", path_out = getwd(),
                  geno = NULL,
                  d = NULL,
                  u = NULL,
                  lmm = 2,
                  miss = 0.05,
                  maf = 0.01,
                  r2 = 0.9999,
                  trait_col = 1,
                  out_name = NULL,
                  verbose = FALSE) {
  home_path <- getwd()
  setwd(path_out) 
  system(
    command = paste( gemma_path, "--bfile", geno,
                     "-lmm",
                     lmm,
                     "-miss",
                     miss,
                     "-maf",
                     maf,
                     "-r2",
                     r2,
                     "-n",
                     paste(trait_col, collapse = " "), "-d",
                     d,
                     "-u",
                     u,
                     "-o",
                     out_name
    ),
    ignore.stdout = !verbose, ignore.stderr = !verbose
  )
  setwd(home_path) }