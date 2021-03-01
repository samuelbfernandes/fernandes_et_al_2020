library(data.table)
library(tidyverse)
library(SNPRelate)
#---------- data summary -----------
path = "./simulation/P/"
folder <- dir(path)[grep("ames_|soy_", dir(path))]
lfolder <- length(folder)
dist_threshold_ames <- 10000
dist_threshold_soy <- 1000000
count_total <- vector("list", lfolder)
total_33 <- tibble(
  specie = character(lfolder),
  sample_size = character(lfolder),
  maf = character(lfolder),
  h2 = "0.3_0.3",
  mt = numeric(lfolder),
  t1 = numeric(lfolder),
  t2 = numeric(lfolder),
  t1_and_t2 = numeric(lfolder),
  af1 = numeric(lfolder))
total_38 <- total_33
total_38$h2 <- "0.3_0.8"
total_88 <- total_33
total_88$h2 <- "0.8_0.8"
for (z in seq_along(folder)){
  gemma_results <- dir(paste0(path,folder[z], "/gemma/output"))
  gemma_gwas_MT <- gemma_results[!grepl("log.txt",gemma_results) & grepl("MTRep",gemma_results)]
  gemma_gwas_t1 <- gemma_results[!grepl("log.txt",gemma_results) & grepl("T1_",gemma_results)]
  gemma_gwas_t2 <- gemma_results[!grepl("log.txt",gemma_results) & grepl("T2_",gemma_results)]
  n_gwas <- length(gemma_gwas_MT)
  qtn <- fread(
    paste0(path, folder[z], "/Additive_selected_QTNs.txt"),
    header = TRUE,
    sep = "\t",
    data.table = F
  )
  d <- dir()[grepl(".gds", dir())]
  gds <- d[grepl(gsub("_MAF.*", "", folder[z]), d)]
  count <- tibble(
    specie = gsub("_.*", "", folder[z]),
    sample_size = gsub(".*_", "", gsub("_MAF.*", "", folder[z])),
    maf = gsub(".*MAF", "", folder[z]),
    h2 = character(n_gwas),
    rep = numeric(n_gwas),
    mt = numeric(n_gwas),
    t1 = numeric(n_gwas),
    t2 = numeric(n_gwas),
    t1_and_t2 = numeric(n_gwas),
    af1 = numeric(n_gwas))
  count_all <- count
  dist <- ifelse(gsub("_.*","",folder[z]) == "ames", dist_threshold_ames, dist_threshold_soy)
  total_33$specie[z] <- gsub("_.*", "", folder[z])
  total_33$sample_size[z] <- gsub(".*_", "", gsub("_MAF.*", "", folder[z]))
  total_33$maf[z] <- gsub(".*MAF", "", folder[z])
  total_38$specie[z] <- total_33$specie[z]
  total_38$sample_size[z] <- total_33$sample_size[z]
  total_38$maf[z] <- total_33$maf[z]
  total_88$specie[z] <- total_33$specie[z]
  total_88$sample_size[z] <- total_33$sample_size[z]
  total_88$maf[z] <- total_33$maf[z]
  
  total_all_33 <- total_33
  total_all_38 <- total_38
  total_all_88 <- total_88
  for (n in 1:n_gwas){
    rep <- gsub("MTRep|_Herit.*", "", gemma_gwas_MT[n])
    h2 <- gsub(".*Herit_|.assoc.*", "", gemma_gwas_MT[n])
    gemma_mt <- fread(
      paste0(path, folder[z], "/gemma/output/", gemma_gwas_MT[n]),
      header = TRUE,
      sep = "\t",
      select = c("chr", "rs", "ps", "p_lrt"),
      data.table = F
    )
    gemma_mt$fdr <- p.adjust(gemma_mt$p_lrt , method = "BH")
    gemma_t1 <- fread(
      paste0(path, folder[z], "/gemma/output/", paste0("T1_", gsub("MT", "", gemma_gwas_MT[n]))),
      header = TRUE,
      sep = "\t",
      select = c("chr", "rs", "ps", "p_lrt"),
      data.table = F
    )
    gemma_t1$fdr <- p.adjust(gemma_t1$p_lrt , method = "BH")
    gemma_t2 <- fread(
      paste0(path, folder[z], "/gemma/output/",  paste0("T2_", gsub("MT", "", gemma_gwas_MT[n]))),
      header = TRUE,
      sep = "\t",
      select = c("chr", "rs", "ps", "p_lrt"),
      data.table = F
    )
    gemma_t2$fdr <- p.adjust(gemma_t2$p_lrt , method = "BH")
    genofile <- snpgdsOpen(gds)
    if (total_33$specie[z] == "ames") {
      snps <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
    } else {
      snps <- read.gdsn(index.gdsn(genofile, "snp.id"))
    }
    qtn1 <- qtn[qtn$rep == rep, ]
    qtn1_genotype <-
      read.gdsn(
        index.gdsn(genofile, "genotype"),
        start = c(1, which(snps%in%qtn1$snp)),
        count = c(-1, 1)
      )
    snpgdsClose(genofile)
    gemma_mt <-
      gemma_mt %>%
      arrange(fdr) %>%
      filter(fdr < 0.05 & chr == qtn1$chr) %>%
      mutate(
        dist_qtn1 = abs(ps - qtn1$pos),
        rep = rep,
        h2 = h2,
        model = "MT"
      )
    gemma_t1 <-
      gemma_t1 %>%
      arrange(fdr) %>%
      filter(fdr < 0.05 & chr == qtn1$chr) %>%
      mutate(
        dist_qtn1 = abs(ps - qtn1$pos),
        rep = rep,
        h2 = h2,
        model = "T1"
      )
    gemma_t2 <-
      gemma_t2 %>%
      arrange(fdr) %>%
      filter(fdr < 0.05 & chr == qtn1$chr) %>%
      mutate(
        dist_qtn1 = abs(ps - qtn1$pos),
        rep = rep,
        h2 = h2,
        model = "T2"
      )
    mt <- gemma_mt %>% filter(dist_qtn1 <= dist)
    t1 <- gemma_t1 %>% filter(dist_qtn1 <= dist)
    t2 <- gemma_t2 %>% filter(dist_qtn1 <= dist)

    count$mt[n] <- as.numeric(nrow(mt) > 0)
    count$t1[n] <- as.numeric(nrow(t1) > 0)
    count$t2[n] <- as.numeric(nrow(t2) > 0)
    count$t1_and_t2[n] <- as.numeric(nrow(t1) > 0 & nrow(t2) > 0)
    count$h2[n] <- h2
    count$rep[n] <- rep
    count$af1[n] <- qtn1$maf

    mt2 <- gemma_mt %>% filter(dist_qtn1 >= dist)
    t12 <- gemma_t1 %>% filter(dist_qtn1 >= dist)
    t22 <- gemma_t2 %>% filter(dist_qtn1 >= dist)
    
    count_all$mt[n] <- as.numeric(nrow(mt2) > 0)
    count_all$t1[n] <- as.numeric(nrow(t12) > 0)
    count_all$t2[n] <- as.numeric(nrow(t22) > 0)
    count_all$t1_and_t2[n] <- as.numeric(nrow(t12) > 0 & nrow(t22) > 0)
    count_all$h2[n] <- h2
    count_all$rep[n] <- rep
    count_all$af1[n] <- qtn1$maf
    
    print(n)
    #gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
  }
  print(folder[z])
  count_h33 <-  count %>% filter(h2 == "0.3_0.3")
  count_h38 <-  count %>% filter(h2 == "0.3_0.8")
  count_h88 <-  count %>% filter(h2 == "0.8_0.8")
  
  count_all_h33 <-  count_all %>% filter(h2 == "0.3_0.3")
  count_all_h38 <-  count_all %>% filter(h2 == "0.3_0.8") 
  count_all_h88 <-  count_all %>% filter(h2 == "0.8_0.8") 
  
  total_33$mt[z] <- count_h33 %>% select(mt) %>% sum
  total_33$t1[z] <- count_h33 %>% select(t1) %>% sum
  total_33$t2[z] <- count_h33 %>% select(t2) %>% sum
  total_33$t1_and_t2[z] <- count_h33 %>% select(t1_and_t2) %>% sum
  total_33$af1[z] <-  count_h33 %>% summarise(af1 = median(af1))
  #----
  total_38$mt[z] <- count_h38 %>% select(mt) %>% sum
  total_38$t1[z] <- count_h38 %>% select(t1) %>% sum
  total_38$t2[z] <- count_h38 %>% select(t2) %>% sum
  total_38$t1_and_t2[z] <- count_h38 %>% select(t1_and_t2) %>% sum
  total_38$af1[z] <-  count_h38 %>% summarise(af1 = median(af1))
  #----
  total_88$mt[z] <- count_h88 %>% select(mt) %>% sum
  total_88$t1[z] <- count_h88 %>% select(t1) %>% sum
  total_88$t2[z] <- count_h88 %>% select(t2) %>% sum
  total_88$t1_and_t2[z] <- count_h88 %>% select(t1_and_t2) %>% sum
  total_88$af1[z] <-  count_h88 %>% summarise(af1 = median(af1))
  count_total[[z]] <- count 
  
  #----------
  total_all_33$mt[z] <- count_all_h33 %>% select(mt) %>% sum
  total_all_33$t1[z] <- count_all_h33 %>% select(t1) %>% sum
  total_all_33$t2[z] <- count_all_h33 %>% select(t2) %>% sum
  total_all_33$t1_and_t2[z] <- count_all_h33 %>% select(t1_and_t2) %>% sum
  total_all_33$af1[z] <-  count_all_h33 %>% summarise(af1 = median(af1))
  #----
  total_all_38$mt[z] <- count_all_h38 %>% select(mt) %>% sum
  total_all_38$t1[z] <- count_all_h38 %>% select(t1) %>% sum
  total_all_38$t2[z] <- count_all_h38 %>% select(t2) %>% sum
  total_all_38$t1_and_t2[z] <- count_all_h38 %>% select(t1_and_t2) %>% sum
  total_all_38$af1[z] <-  count_all_h38 %>% summarise(af1 = median(af1))

  total_all_88$mt[z] <- count_all_h88 %>% select(mt) %>% sum
  total_all_88$t1[z] <- count_all_h88 %>% select(t1) %>% sum
  total_all_88$t2[z] <- count_all_h88 %>% select(t2) %>% sum
  total_all_88$t1_and_t2[z] <- count_all_h88 %>% select(t1_and_t2) %>% sum
  total_all_88$af1[z] <-  count_all_h88 %>% summarise(af1 = median(af1))
}

count_total <- rbindlist(count_total)
total <- bind_rows(total_33, total_38, total_88)
total_all <- bind_rows(total_all_33, total_all_38, total_all_88)
# fwrite(count_total, "total_count_all_P.txt", sep = "\t", quote = F, row.names = F)
# fwrite(total, "total_count_P.txt", sep = "\t", quote = F, row.names = F)

path = "./simulation/P/"
folder <- dir(path)[grep("ames_|soy_", dir(path))]
dist_threshold_ames <- 1000
dist_threshold_soy <- 10000
ld1 <- list()
ld2 <- list()
for (z in seq_along(folder)){
  ld1[[z]] <- tibble(
    specie = gsub("_.*", "", folder[z]),
    sample_size = gsub(".*_", "", gsub("_MAF.*", "", folder[z])),
    maf = gsub(".*MAF", "", folder[z]))
  ld2[[z]] <- tibble(
    specie = gsub("_.*", "", folder[z]),
    sample_size = gsub(".*_", "", gsub("_MAF.*", "", folder[z])),
    maf = gsub(".*MAF", "", folder[z]))
  
  gemma_results <- dir(paste0(path,folder[z], "/gemma/output"))
  gemma_gwas_MT <- gemma_results[!grepl("log.txt",gemma_results) & grepl("MTRep",gemma_results)]
  qtn <- fread(
    paste0(path, folder[z], "/Additive_selected_QTNs.txt"),
    header = TRUE,
    sep = "\t",
    data.table = F
  )
  d <- dir()[grepl(".gds", dir())]
  gds <- d[grepl(gsub("_MAF.*", "", folder[z]), d)]
  dist <- ifelse(gsub("_.*","",folder[z]) == "ames", dist_threshold_ames, dist_threshold_soy)
  ldtemp1 <- as.data.frame(matrix(NA, nrow = 100, ncol = 41 ))
  ldtemp2 <- as.data.frame(matrix(NA, nrow = 100, ncol = 41 ))
  for (n in 1:100){
    rep <- gsub("MTRep|_Herit.*", "", gemma_gwas_MT[n])
    genofile <- snpgdsOpen(gds)
    if (gsub("_.*","",folder[z]) == "ames") {
      snps <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
    } else {
      snps <- read.gdsn(index.gdsn(genofile, "snp.id"))
    }
    qtn$h2 <- gsub("_.*", "", qtn$h2)
    for(i in unique(qtn$h2)){
      qtn1 <- qtn[qtn$rep == rep & qtn$snp_type=="QTN_for_trait_1" & qtn$h2 == i, ]
      qtn_position <- which(snps%in%qtn1$snp)      
      snps_in_ld <-  snps[-qtn_position][(qtn_position-20):(qtn_position+20)]
      qtn1_genotype <-
        read.gdsn(
          index.gdsn(genofile, "genotype"),
          start = c(1, qtn_position),
          count = c(-1, 1)
        )
      a <- 1
      for(j in snps_in_ld) {
        qtn2_genotype <-
          read.gdsn(
            index.gdsn(genofile, "genotype"),
            start = c(1, which(snps%in%j)),
            count = c(-1, 1)
          )
        ldtemp1[n, a] <-
          abs(SNPRelate::snpgdsLDpair(qtn1_genotype, qtn2_genotype, method = "composite"))
        a <- a+1
      }
      
      qtn2 <- qtn[qtn$rep == rep & qtn$snp_type=="QTN_for_trait_2" & qtn$h2 == i, ]
      qtn_position <- which(snps%in%qtn2$snp)
      snps_in_ld <-  snps[-qtn_position][(qtn_position-20):(qtn_position+20)]
      qtn1_genotype <-
        read.gdsn(
          index.gdsn(genofile, "genotype"),
          start = c(1, qtn_position),
          count = c(-1, 1)
        )
      a <- 1
      for(j in snps_in_ld) {
        qtn2_genotype <-
          read.gdsn(
            index.gdsn(genofile, "genotype"),
            start = c(1, which(snps%in%j)),
            count = c(-1, 1)
          )
        ldtemp2[n, a] <-
          abs(SNPRelate::snpgdsLDpair(qtn1_genotype, qtn2_genotype, method = "composite"))
        a <- a+1
      }
    }
    snpgdsClose(genofile)
    print(n)
  }
  ld1[[z]] <- tibble(ld1[[z]], ldtemp1)
  ld2[[z]] <- tibble(ld2[[z]], ldtemp2)
  print(folder[z])
}

ld1 <- rbindlist(ld1)
ld2 <- rbindlist(ld2)
ld <- rbind(ld1, ld2)
ld$qtn <- rep(c(1, 2), each = 1200)
ld <- ld[, c(1:5, 45, 6:44)]
#fwrite(ld, "ld_QTN_P.txt", sep = "\t", quote = F, row.names = F)