run_simulation_P <-  function(geno_file = NULL, list = NULL, out, h2, rep) {
  pref <- gsub("../../|_mac.*","", geno_file)
  create_phenotypes(
    geno_file = paste0("data/", geno_file),
    add_QTN_num = 1,
    h2 = h2,
    add_effect = c(0.1, 0.1),
    rep = rep,
    seed = list$seed + as.numeric(gsub(".*_","",pref)),
    output_format = "gemma",
    architecture = "pleiotropic",
    ntraits = 2,
    output_dir = paste0(out, pref, "_MAF", list$maf, collapse = "_"),
    model = "A",
    home_dir = getwd(),
    constraints = list(
      maf_above = list$maf - (list$maf * 0.1),
      maf_below = list$maf + (list$maf * 0.1)
    ),
    quiet = T,
    remove_QTN = TRUE,
    vary_QTN = TRUE,
    warning_file_saver = FALSE,
    verbose = F
  )
}
