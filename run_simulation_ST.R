run_simulation_ST <-  function(geno_file = NULL, list = NULL, out, h2, rep) {
  pref <- gsub("../../|_mac.*","", geno_file)
  new_dir <- paste0(out, pref, "_MAF", list$maf)
  dir.create(new_dir)
  qtn <- vector("list", nrow(h2))
  for(i in 1:nrow(h2)){
    create_phenotypes(
      geno_file = paste0("data/", geno_file),
      add_QTN_num = 1,
      h2 = h2[i,1],
      add_effect = 0.1,
      rep = rep * 2,
      seed = list$seed + as.numeric(gsub(".*_","",pref))+i,
      output_format = "multi-file",
      ntraits = 1,
      output_dir = paste0(new_dir,"/t1",i),
      model = "A",
      home_dir = getwd(),
      constraints = list(
        maf_above = list$maf - (list$maf * 0.1),
        maf_below = list$maf + (list$maf * 0.1)
      ),
      quiet = T,
      vary_QTN = T,
      verbose = F,
      out_geno = "plink"
    )
    create_phenotypes(
      geno_file = paste0("data/", geno_file),
      add_QTN_num = 1,
      h2 = h2[i,2],
      add_effect = 0.1,
      rep = rep * 2,
      seed = round((list$seed + as.numeric(gsub(".*_","",pref)))/2) + i,
      output_format = "multi-file",
      ntraits = 1,
      output_dir = paste0(new_dir,"/t2",i),
      model = "A",
      home_dir = getwd(),
      constraints = list(
        maf_above = list$maf - (list$maf * 0.1),
        maf_below = list$maf + (list$maf * 0.1)
      ),
      quiet = T,
      vary_QTN = T,
      verbose = F
    )
    qtn1 <- fread(paste0(new_dir,"/t1", i,"/Additive_selected_QTNs.txt"), data.table = F)
    qtn1$snp_type <- "QTN_for_trait_1"
    qtn1$h2 <- paste0(i, "_",h2[i,1])
    qtn2 <- fread(paste0(new_dir,"/t2", i,"/Additive_selected_QTNs.txt"), data.table = F)
    qtn2$snp_type <- "QTN_for_trait_2"
    qtn2$h2 <-  paste0(i, "_",h2[i,2])
    qtn1$chr <- as.numeric(gsub("Gm", "", qtn1$chr))
    qtn2$chr <- as.numeric(gsub("Gm", "", qtn2$chr))
    same_chr <- qtn1$chr == qtn2$chr
    qtn1 <- qtn1[!same_chr,]
    qtn2 <- qtn2[!same_chr,]
    qtn[[i]] <- rbind(qtn1, qtn2)
  }
  qtn <- data.table::rbindlist(qtn)
  qtn <- qtn[order(qtn$rep), c("rep",	"snp_type",	"snp",	"allele",	"chr",	"pos",	"cm",	"maf", "h2")]
  
  file.rename(paste0(new_dir,"/t11/",gsub(".txt",".bed", geno_file)), paste0(new_dir,"/",gsub(".txt",".bed", geno_file)))
  file.rename(paste0(new_dir, "/t11/",gsub(".txt",".bim", geno_file)), paste0(new_dir,"/",gsub(".txt",".bim", geno_file)))
  
  qtn <- qtn[qtn$rep %in%  names(table(qtn$rep)==6)[table(qtn$rep)==6], ]
  reps <- unique(qtn$rep)
  reps <- reps[1:100]
  qtn <- qtn[qtn$rep %in% reps,]
  qtn$rep <- rep(1:100, each = 6)
  fwrite(qtn, file = paste0(new_dir,"/Additive_selected_QTNs.txt"),
         quote = F, sep = "\t", row.names = F)
  a <- 1
  cor <- data.frame(matrix(NA, 100, 3))
  colnames(cor) <- c("h1", "h2", "h3")
  for (i in seq_len(rep*2)){
    if (!i %in% reps ) {next}
    t1 <- fread(paste0(new_dir,"/t11/Simulated_Data_Rep", i, "_Herit_",h2[1,1],".txt"), data.table = F)
    t2 <- fread(paste0(new_dir,"/t21/Simulated_Data_Rep", i, "_Herit_",h2[1,2],".txt"), data.table = F)
    cor$h1[a] <- cor(t1[,2], t2[,2])
    fam1 <- data.frame(t1[,1], t1[,1], 0, 0, 0, t1[,2])
    colnames(fam1) <- paste0("V", 1:6)
    fam1 <- merge(fam1, t2, by.x = "V1", by.y = "<Trait>", sort = F)
    colnames(fam1)[7] <- "V7"
    fwrite(fam1, file = paste0(new_dir,"/","Simulated_Data__Rep", a,"_Herit_",h2[1,1],"_",h2[1,2],".fam"),
           quote = F, sep = "\t", row.names = F, col.names = F)
    
    t3 <- fread(paste0(new_dir, "/t12/Simulated_Data_Rep", i, "_Herit_",h2[2,1],".txt"), data.table = F)
    t4 <- fread(paste0(new_dir, "/t22/Simulated_Data_Rep", i, "_Herit_",h2[2,2],".txt"), data.table = F)
    cor$h2[a] <- cor(t3[,2], t4[,2])   
    fam2 <- data.frame(t3[,1], t3[,1], 0, 0, 0, t3[,2])
    colnames(fam2) <- paste0("V", 1:6)
    fam2 <- merge(fam2, t4, by.x = "V1", by.y = "<Trait>", sort = F)
    colnames(fam2)[7] <- "V7"
    fwrite(fam2, file = paste0(new_dir,"/","Simulated_Data__Rep", a,"_Herit_",h2[1,2],"_",h2[2,2],".fam"),
           quote = F, sep = "\t", row.names = F, col.names = F)
    
    t5 <- fread(paste0(new_dir,"/t13/Simulated_Data_Rep", i, "_Herit_",h2[3,1],".txt"), data.table = F)
    t6 <- fread(paste0(new_dir,"/t23/Simulated_Data_Rep", i, "_Herit_",h2[3,2],".txt"), data.table = F)
    cor$h3[a] <- cor(t5[,2], t6[,2])
    fam3 <- data.frame(t5[,1], t5[,1], 0, 0, 0, t5[,2])
    colnames(fam3) <- paste0("V", 1:6)
    fam3 <- merge(fam3, t6, by.x = "V1", by.y = "<Trait>", sort = F)
    colnames(fam3)[7] <- "V7"
    fwrite(fam3, file = paste0(new_dir,"/","Simulated_Data__Rep", a,"_Herit_",h2[3,1],"_",h2[3,2],".fam"),
           quote = F, sep = "\t", row.names = F, col.names = F)
    a <- a+1
  }
  fwrite(cor, file = paste0(new_dir,"/correlation.txt"),
         quote = F, sep = "\t", row.names = F)
  unlink(c(paste0(new_dir, "/t11"),
           paste0(new_dir, "/t12"),
           paste0(new_dir, "/t21"),
           paste0(new_dir, "/t22"),
           paste0(new_dir, "/t13"),
           paste0(new_dir, "/t23")),
         recursive = T,
         force = T)
}
