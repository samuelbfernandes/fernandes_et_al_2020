library(ggplot2)
library(cowplot)
library(viridis)
library(tidyverse)
library(data.table)
st  <- fread("total_count_all_ST.txt", data.table = F)
ld  <- fread("total_count_all_LD.txt", data.table = F) %>% 
  mutate(maf = gsub("_.*", "", maf)) %>% group_by(specie, sample_size, ld, maf, type, h2) %>% mutate(rep = order(rep)) %>%
  ungroup()
p  <- fread("total_count_all_P.txt", data.table = F)
cor_st  <- fread("correlation_ST.txt", data.table = F)
cor_ld  <- fread("correlation_LD.txt", data.table = F) %>% 
  mutate(maf = as.character(maf)) %>% select(specie, sample_size, maf, h2, ld, type, cor)
cor_p  <- fread("correlation_P.txt", data.table = F)

cor_st$rep <- st$rep
cor_ld$rep <- ld$rep
cor_p$rep <- p$rep

st <- st %>% left_join(cor_st %>% select(specie, sample_size, maf, rep, h2, cor), by = c("specie", "sample_size", "maf","rep", "h2"))
p  <- p %>% left_join(cor_p %>% select(specie, sample_size, maf, rep, h2, cor), by = c("specie", "sample_size", "maf", "rep", "h2"))
ld <- ld %>% right_join(cor_ld, by = c("specie", "sample_size", "ld", "maf", "type", "h2", "rep"))

#---- Fig MAF_LD (FIG2)----
st1 <- st %>%  select(-chr1, -chr2, -dist_QTN) %>%  mutate(ld_QTN = (ld_QTN)**2, cor = abs(cor)) %>%
  pivot_longer(names_to = "model", values_to = "count", cols = mt:t1_and_t2) %>%
  pivot_longer(names_to = "af", values_to = "af_value", cols = af1:cor) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multivariate", "Univariate"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT","mtt1" = "MTinT1","mtt2" = "MTinT2","mtt1t2" = "MTinT1&T2", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2" ),
         x = recode(model, "MTinT1" = "T1", "MTinT2" = "T2",  "MTinT1&T2" = "T1&T2"),
         model = fct_relevel(model, "MT", "MTinT1", "MTinT2", "MTinT1&T2", "T1", "T2", "T1&T2"),
         x = fct_relevel(x, "T1", "T2", "T1&T2"),
         sample_size = fct_relevel(as.character(sample_size), "500", "1000", "2815"),
         maf = as.character(maf),
         alpha = "MAF") %>%
  rename("Specie" = "specie", "MAF" = "maf")


ld1 <- ld %>%  select(-dist_QTN) %>% pivot_longer(names_to = "model", values_to = "count", cols = mt:t1_and_t2) %>%
  mutate(ld_QTN = (ld_QTN)**2, cor = abs(cor)) %>%
  pivot_longer(names_to = "af", values_to = "af_value", cols = af1:cor) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multivariate", "Univariate"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         type = recode(type, "direct" = "Direct", "indirect" = "Indirect"),
         maf = recode(maf, "0.05_direct" = "0.05", "0.05_indirect" = "0.05","0.4_direct" = "0.4", "0.4_indirect" = "0.4"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT","mtt1" = "MTinT1","mtt2" = "MTinT2","mtt1t2" = "MTinT1&T2", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2" ),
         x = recode(model, "MTinT1" = "T1", "MTinT2" = "T2",  "MTinT1&T2" = "T1&T2"),
         model = fct_relevel(model, "MT", "MTinT1", "MTinT2", "MTinT1&T2", "T1", "T2", "T1&T2"),
         x = fct_relevel(x, "T1", "T2", "T1&T2"),
         sample_size = as.factor(sample_size),
         sample_size = fct_relevel(sample_size, "500", "1000", "2815"),
         ld = as.character(ld), maf = as.character(maf),
         alpha = ifelse((type == "Direct" & af == "af1") | af == "ld_QTN" |  af == "cor", "MAF", "No MAF")) %>%
  rename("Specie" = "specie", "LD" = "ld", "MAF" = "maf")


p1 <- p %>%  pivot_longer(names_to = "model", values_to = "count", cols = mt:t1_and_t2) %>%
  pivot_longer(names_to = "af", values_to = "af_value", cols = af1:cor) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multivariate", "Univariate"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2"),
         model = fct_relevel(model, "MT","T1", "T2", "T1&T2"),
         sample_size = fct_relevel(as.character(sample_size), "500", "1000", "2815"),
         maf = as.character(maf),
         alpha = "MAF") %>%
  rename("Specie" = "specie", "MAF" = "maf")

all_af <- bind_rows(LD = ld1 %>% select(Specie, sample_size, LD, MAF, type, h2, af, af_value,  x, alpha, count),
                    ST = st1 %>% mutate(type = "Independent", LD = "") %>% select(Specie, sample_size, LD, MAF, type, h2, af, af_value,  x, alpha, count),
                    P = p1 %>% mutate(x = "MT", type = "Pleiotropy", LD = "") %>%
                      select(Specie, sample_size, LD, MAF, type, h2, af, af_value,  x, alpha, count), .id = "architecture") %>%
  mutate(vars = paste0(type, LD),
         vars = recode(vars, "Direct0.1" = "Direct LD: 0.01", "Direct0.99" = "Direct LD: 0.98", "Indirect0.1" = "Indirect LD: 0.01", "Indirect0.99" = "Indirect LD: 0.98"),
         vars = fct_relevel(vars, "Independent", "Indirect LD: 0.01", "Indirect LD: 0.98", "Direct LD: 0.01", "Direct LD: 0.98", "Pleiotropy"),
         af = fct_relevel(af, "af1", "af2", "ld_QTN", "cor"))

labs <- c(expression(QTN['T1']),
          expression(QTN['T2']), 
          expression('LD Between QTNs'[('r'^'2')]),
          expression('|r|'))
labs2 <- c("Yes", "No")

all_af$background <- ifelse(all_af$architecture == "LD", "A", "B")

all_af %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "1000") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c("#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
    ) + guides( fill = guide_legend(title ="Observed MAF, LD, and |r| distributions:", nrow = 1, byrow = TRUE, order = 2))

ggsave(file="./figures/MAF_LD.pdf", width = 30, height = 25, units = "cm", dpi=300,family="Times")

#---- supplementary figure 1-----
a <- 
  all_af %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "1000") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c( "#0000EE", "#006400","#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + theme(legend.position = "none")

b <- all_af %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "1000") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c( "#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed MAF, LD, and |r| distributions:", nrow = 1, byrow = TRUE, order = 2)
  )
plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.35)
)

ggsave(file="./figures/sup1_maf_ld_1000.pdf", width = 30, height = 40, units = "cm", dpi=300,family="Times")

#---- Sup1 2 -----
a2 <- all_af %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "500") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual(values=c("#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) +   theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + theme(legend.position = "none")

b2 <- all_af %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "2815") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c( "#0000EE", "#006400","#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed MAF, LD, and |r| distributions:", nrow = 1, byrow = TRUE, order = 2)
  )
plot_grid(
  a2, b2,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.36)
)

ggsave(file="./figures/sup1_maf_ld_h2_03.pdf", width = 30, height = 40, units = "cm", dpi=300,family="Times")

#---- Sup3 -----
a3 <- all_af %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "500") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual(values=c("#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) +  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + theme(legend.position = "none")

b3 <- all_af %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "2815") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c( "#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed MAF, LD, and |r| distributions:", nrow = 1, byrow = TRUE, order = 2)
  )
plot_grid(
  a3, b3,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.36)
)

ggsave(file="./figures/sup1_maf_ld_h2_0308.pdf", width = 30, height = 40, units = "cm", dpi=300,family="Times")
#---- Sup4 -----
a4 <- all_af %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "500") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual( values=c( "#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) +  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + theme(legend.position = "none")

b4 <- all_af %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "2815") %>%
  ggplot(aes(x = sample_size, fill = af, y = af_value, alpha = alpha)) +
  scale_fill_manual(values=c( "#0000EE", "#006400", "#FFE4C4", "#FF6EB4"), labels = labs) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) + ggtitle(~""*underline("                                         QTNs in Linkage                                              ")
  ) + theme(plot.title = element_text(size = 16, hjust = .48, vjust = -4)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed MAF, LD, and |r| distributions:", nrow = 1, byrow = TRUE, order = 2)
  )
plot_grid(
  a4, b4,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.36)
)

ggsave(file="./figures/sup1_maf_ld_h2_08.pdf", width = 30, height = 40, units = "cm", dpi=300,family="Times")

#Fig 3 QTN detection----------
st  <- fread("total_count_ST.txt", data.table = F)
ld  <- fread("total_count_LD.txt", data.table = F) 
p  <- fread("total_count_P.txt", data.table = F)

st1 <- st %>%  select(-chr1, -chr2) %>%  mutate(ld_QTN = (ld_QTN)**2) %>%
  pivot_longer(names_to = "model", values_to = "count", cols = mt:t1_and_t2) %>%
  pivot_longer(names_to = "af", values_to = "af_value", cols = af1:ld_QTN) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multi", "Uni"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT","mtt1" = "MTinT1","mtt2" = "MTinT2","mtt1t2" = "MTinT1&T2", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2" ),
         x = recode(model, "MTinT1" = "T1", "MTinT2" = "T2",  "MTinT1&T2" = "T1&T2"),
         model = fct_relevel(model, "MT", "MTinT1", "MTinT2", "MTinT1&T2", "T1", "T2", "T1&T2"),
         x = fct_relevel(x, "T1", "T2", "T1&T2"),
         sample_size = fct_relevel(as.character(sample_size), "500", "1000", "2815"),
         maf = as.character(maf),
         alpha = ifelse(af != "ld_QTN", "MAF", "No MAF")) %>%
  rename("Specie" = "specie", "MAF" = "maf")


ld1 <- ld %>%  pivot_longer(names_to = "model", values_to = "count", cols = mt:t1_and_t2) %>%
  mutate(ld_QTN = (ld_QTN)**2) %>%
  pivot_longer(names_to = "af", values_to = "af_value", cols = af1:ld_QTN) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multi", "Uni"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         type = recode(type, "direct" = "Direct", "indirect" = "Indirect"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT","mtt1" = "MTinT1","mtt2" = "MTinT2","mtt1t2" = "MTinT1&T2", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2" ),
         x = recode(model, "MTinT1" = "T1", "MTinT2" = "T2",  "MTinT1&T2" = "T1&T2"),
         model = fct_relevel(model, "MT", "MTinT1", "MTinT2", "MTinT1&T2", "T1", "T2", "T1&T2"),
         x = fct_relevel(x, "T1", "T2", "T1&T2"),
         sample_size = fct_relevel(as.character(sample_size), "500", "1000", "2815"),
         ld = as.character(ld), maf = as.character(maf),
         alpha = ifelse((type == "Direct" & af == "af1") | af == "ld_QTN", "MAF", "No MAF")) %>%
  rename("Specie" = "specie", "LD" = "ld", "MAF" = "maf")


p1 <- p %>%  mutate(mtt1=NA, mtt2 =NA, mtt1t2=NA ) %>% select(-af1) %>% pivot_longer(names_to = "model", values_to = "count", cols = mt:mtt1t2) %>%
  mutate(ntraits = ifelse(grepl("mt", model), "Multi", "Uni"),
         h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         model = recode(model, "mt" = "MT","mtt1" = "MTinT1","mtt2" = "MTinT2","mtt1t2" = "MTinT1&T2", "t1_and_t2" = "T1&T2", "t1" = "T1", "t2" = "T2" ),
         x = recode(model, "MTinT1" = "T1", "MTinT2" = "T2",  "MTinT1&T2" = "T1&T2"),
         model = fct_relevel(model, "MT", "MTinT1", "MTinT2", "MTinT1&T2", "T1", "T2", "T1&T2"),
         x = fct_relevel(x, "T1", "T2", "T1&T2"),
         sample_size = fct_relevel(as.character(sample_size), "500", "1000", "2815"),
         maf = as.character(maf),
         alpha = "MAF") %>%
  rename("Specie" = "specie", "MAF" = "maf")

all <- bind_rows(LD = ld1 %>% select(Specie, sample_size, LD, MAF, type, h2, model, count, ntraits, x),
                 ST = st1 %>% mutate(type = "Independent", LD = "") %>% select(Specie, sample_size, LD, MAF, type, h2, model, count, ntraits, x),
                 P = p1 %>% mutate(type = "Pleiotropy", LD = "") %>%
                   select(Specie, sample_size, LD, MAF, type, h2, model, count, ntraits, x), .id = "architecture") %>%
  mutate(vars = paste0(type, LD),
         vars = recode(vars, "Direct0.1" = "Direct LD: 0.01", "Direct0.99" = "Direct LD: 0.98", "Indirect0.1" = "Indirect LD: 0.01", "Indirect0.99" = "Indirect LD: 0.98"),
         vars = fct_relevel(vars, "Independent", "Indirect LD: 0.01", "Indirect LD: 0.98", "Direct LD: 0.01", "Direct LD: 0.98", "Pleiotropy"))

all$background <- ifelse(all$architecture == "LD", "A", "B")
all <- all %>% mutate(MAF = ifelse(MAF == 0.4, "0.40", "0.05"))  
#------------------
xlabs <-  c(expression(0.3-0.3),
            expression(0.3-0.8),
            expression(0.8-0.8))

maf1 <-
  all %>% filter(Specie == "Maize",  ifelse(x == "MT" & architecture != "P", F, T)) %>% 
  filter( sample_size == "1000", MAF == "0.05") %>%
  ggplot(aes(x = h2, y = count, group = x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels = xlabs) +
  facet_grid(ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "none"
  )

maf2 <-
  all %>% filter(Specie == "Maize",  ifelse(x == "MT" & architecture != "P", F, T)) %>% 
  filter( sample_size == "1000", MAF == "0.40") %>%
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                               ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "bottom"
  ) + guides(fill=FALSE,  colour = guide_legend(nrow = 1, override.aes = list(size=2)))

plot_grid(
  maf1, maf2,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.15)
)


ggsave(file="./figures/h2_maf.pdf", width = 32, height = 30, units = "cm", dpi=300,family="Times")
#--------------------

maf1 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>% 
  filter(sample_size == "500", MAF == "0.05") %>%
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "none"
  )

maf2 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>%
  filter(sample_size == "500", MAF == "0.05") %>%
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "bottom"
  ) + guides(fill=FALSE, colour = guide_legend(nrow = 1, override.aes = list(size=2)))

plot_grid(
  maf1, maf2,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file="./figures/sup2_h2_maf_500.pdf", width = 35, height = 35, units = "cm", dpi=300,family="Times")


#--------------------
maf1 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>% 
  filter(sample_size == "1000", MAF == "0.05") %>% 
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "none"
  )

maf2 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>%
  filter(sample_size == "1000", MAF == "0.40") %>% 
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "bottom"
  ) + guides(fill=FALSE, colour = guide_legend(nrow = 1, override.aes = list(size=2)))

plot_grid(
  maf1, maf2,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file="./figures/sup2_h2_maf_1000.pdf", width = 35, height = 35, units = "cm", dpi=300,family="Times")


#--------------------
maf1 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>%
  filter(sample_size == "2815", MAF == "0.05") %>%
  ggplot( aes(x=h2, y=count, group=x,  color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "none"
  )

maf2 <-
  all %>% filter(ifelse(x == "MT" &
                          architecture != "P", F, T)) %>% 
  filter(sample_size == "2815", MAF == "0.40") %>%
  ggplot( aes(x=h2, y=count, group=x, color = x)) +
  geom_line() +
  geom_point() +
  scale_fill_manual(values=c("#000000","#FFFFFF")) +
  scale_color_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c("#008B00", "#0000FF", "#EEEE00", "#68228B"))+
  scale_x_discrete(labels=xlabs)+
  facet_grid(Specie + ntraits ~ vars, labeller = labeller(h2 = label_parsed)) +
  xlab("Heritability") + ylim(0, 100) +
  ylab("Values (%)") + theme_bw(base_size = 20) + 
  ggtitle( ~ "" * underline("                                             QTNs in Linkage                                                 ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "horizontal",
    legend.position = "bottom"
  ) + guides(fill=FALSE, colour = guide_legend(nrow = 1, override.aes = list(size=2)))

plot_grid(
  maf1, maf2,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file="./figures/sup2_h2_maf_2815.pdf", width = 35, height = 35, units = "cm", dpi=300,family="Times")

#-----------Sup3 ---------
all %>% filter(
  h2 == "h^2:~0.3~-~0.8",
  Specie == "Maize",
  sample_size == "1000",
  ifelse(x == "MT" & architecture != "P", F, T),
  (vars  == "Independent" | vars == "Pleiotropy" | vars == "Direct LD: 0.01" |vars  == "Direct LD: 0.98" ),
  !(vars  == "Direct LD: 0.98" & x == "T2"),
  !(vars  == "Direct LD: 0.01" & x == "T2"),
  x != "T1&T2"
) %>% 
  ggplot(aes(x = x, fill = MAF, y = count)) +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(ntraits ~ vars, labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100)+
  ggtitle( ~ "" * underline("             QTNs in Linkage          ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="MAF:", nrow = 1, byrow = TRUE, order = 2) ,vjust = 0.5 )

ggsave(file="./figures/Maize_MAF.pdf", width = 30, height = 20, units = "cm", dpi=300,family="Times")
#--------------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.3",Specie == "Maize",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"), x != "T1&T2")%>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.position = "none"
  )
b <- all %>% filter(h2 == "h^2:~0.3~-~0.3",Specie == "Soybean",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"), x != "T1&T2")%>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="MAF:", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1,  1.01)
)

ggsave(file="./figures/sup3_MAF_h2_0303.pdf", width = 35, height = 40, units = "cm", dpi=300,family="Times")

#--------------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.8",Specie == "Maize",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"),
                    x != "T1&T2")%>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )

b <- all %>% filter(h2 == "h^2:~0.3~-~0.8",Specie == "Soybean",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"), x != "T1&T2")%>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="MAF:", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1,  1.01)
)

ggsave(file="./figures/sup3_MAF_h20308.pdf", width = 35, height = 40, units = "cm", dpi=300,family="Times")

#--------------------
a <- all %>% filter(h2 == "h^2:~0.8~-~0.8",Specie == "Maize",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"),  x != "T1&T2") %>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )

b <- all %>% filter(h2 == "h^2:~0.8~-~0.8",Specie == "Soybean",
                    ifelse(x == "MT" & architecture != "P", F, T), (vars  == "Independent"| vars == "Pleiotropy" | vars == "Direct LD: 0.01" | vars  == "Direct LD: 0.98"), !(vars  == "Direct LD: 0.98" & x == "T2"), !(vars  == "Direct LD: 0.01" & x == "T2"),  x != "T1&T2")%>%
  ggplot(aes(x = x, fill = MAF, y = count)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=c("0.05", "0.40"), values=c("#FFFF00", "#551A8B", "#000000","#FFFFFF")) +
  facet_grid(ntraits+ sample_size ~ vars , labeller = labeller(MAF = label_both), scales = "free", space="free") +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  ggtitle( ~ "" * underline("               QTNs in Linkage            ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .405, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="MAF:", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.01)
)

ggsave(file="./figures/sup3_MAF_h2_0808.pdf", width = 35, height = 40, units = "cm", dpi=300,family="Times")

#-----------------
all %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "1000", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./figures/LD_effect.pdf", width = 33, height = 20, units = "cm", dpi=300,family="Times")

#-----------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "500", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )
b <- all %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == "2815", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.15)
)

ggsave(file="./figures/sup4_LD_effect_h2_0308.pdf", width = 33, height = 40, units = "cm", dpi=300,family="Times")

#-----------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "500", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )

b <- all %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "500", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.15)
)

ggsave(file="./figures/sup4_LD_effect_h2_500.pdf", width = 33, height = 40, units = "cm", dpi=300,family="Times")

#-----------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "1000", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )

b <- all %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "1000", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.15)
)

ggsave(file="./figures/sup4_LD_effect_h2_1000.pdf", width = 33, height = 40, units = "cm", dpi=300,family="Times")

#-----------------
a <- all %>% filter(h2 == "h^2:~0.3~-~0.3", sample_size == "2815", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "none"
  )

b <- all %>% filter(h2 == "h^2:~0.8~-~0.8", sample_size == "2815", ifelse(x == "MT" & architecture != "P", F, T))%>%
  ggplot(aes(x = ntraits, fill = x, y = count)) +
  scale_fill_manual(breaks=c("T1", "T2", "T1&T2", "MT"), values=c( "#551A8B", "#458B74", "#00CD00", "#FFFF00")) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Values (%)") +  theme_bw(base_size = 20) +
  ggtitle( ~ "" * underline("                                         QTNs in Linkage                                              ")) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

plot_grid(
  a, b,
  labels = "AUTO", ncol = 1, rel_heights = c(1, 1.15)
)

ggsave(file="./figures/sup4_LD_effect_h2_2815.pdf", width = 33, height = 40, units = "cm", dpi=300,family="Times")

#--- summary of LD -------

st  <- fread("ld_QTN_ST.txt", data.table = F)
ld  <- fread("ld_QTN_LD.txt", data.table = F)
p  <- fread("ld_QTN_P.txt", data.table = F)

st1 <- st %>%   pivot_longer(names_to = "snps", values_to = "ld_QTN", cols = V1:V41) %>%  mutate(ld_QTN = (ld_QTN)**2) %>%
  mutate(specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         sample_size = as.factor(sample_size),
         sample_size = fct_relevel(sample_size, "500", "1000", "2815"),
         maf = as.factor(maf),
         alpha = "MAF") %>%
  rename("Specie" = "specie", "MAF" = "maf")


ld1 <- ld %>%   pivot_longer(names_to = "snps", values_to = "ld_QTN", cols = V1:V41) %>%  mutate(ld_QTN = (ld_QTN)**2) %>%
  mutate(type = recode(type, "direct" = "Direct", "indirect" = "Indirect"),
         maf = as.factor(maf),
         specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         sample_size = as.factor(sample_size),
         sample_size = fct_relevel(sample_size, "500", "1000", "2815"),
         ld = as.factor(ld),
         alpha = ifelse((type == "Direct" & qtn == "1"), "MAF", "No MAF")) %>%
  rename("Specie" = "specie", "LD" = "ld", "MAF" = "maf")


p1 <- p  %>%   pivot_longer(names_to = "snps", values_to = "ld_QTN", cols = V1:V41) %>%  mutate(ld_QTN = (ld_QTN)**2) %>%
  mutate(specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
         sample_size = as.factor(sample_size),
         sample_size = fct_relevel(sample_size, "500", "1000", "2815"),
         maf = as.factor(maf),
         alpha = "MAF") %>%
  rename("Specie" = "specie", "MAF" = "maf")

all_af <- bind_rows(LD = ld1,
                    ST = st1 %>% mutate(type = "Independent", LD = ""),
                    P = p1 %>% mutate(qtn = 1, type = "Pleiotropy", LD = ""), .id = "architecture") %>%
  mutate(vars = paste0(type, LD),
         vars = recode(vars, "Direct0.1" = "Direct LD: 0.01", "Direct0.99" = "Direct LD: 0.81", "Indirect0.1" = "Indirect LD: 0.01", "Indirect0.99" = "Indirect LD: 0.81"),
         vars = fct_relevel(vars, "Independent", "Indirect LD: 0.01", "Indirect LD: 0.81", "Direct LD: 0.01", "Direct LD: 0.81", "Pleiotropy"),
         qtn = as.factor(qtn))

labs <- c(expression(QTN['T1']),
          expression(QTN['T2']))
labs2 <- c("Yes", "No")

all_af %>% filter(sample_size == "500") %>%
  ggplot(aes(x = sample_size, fill = qtn, y = ld_QTN, alpha = alpha)) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Linkage Disequilibrium") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_fill_manual(values = c("#0000EE", "#006400"), labels = labs) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed LD in the QTN region", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./figure/LD_QTN_500.pdf", width = 30, height = 25, units = "cm", dpi=300,family="Times")


all_af %>% filter(sample_size == "1000") %>%
  ggplot(aes(x = sample_size, fill = qtn, y = ld_QTN, alpha = alpha)) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Linkage Disequilibrium") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_fill_manual(values = c("#0000EE", "#006400"), labels = labs) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed LD in the QTN region", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./figure/LD_QTN_1000.pdf", width = 30, height = 25, units = "cm", dpi=300,family="Times")



all_af %>% filter(sample_size == "2815") %>%
  ggplot(aes(x = sample_size, fill = qtn, y = ld_QTN, alpha = alpha)) +
  geom_boxplot(lwd = 0.05, position = position_dodge2(preserve = "single"), outlier.colour = c("#ADADAD")) +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("") +  ylab("Linkage Disequilibrium") +  theme_bw(base_size = 20) + ylim(0, 1) +
  scale_fill_manual(values = c("#0000EE", "#006400"), labels = labs) +
  scale_alpha_manual(values = c(0.7, 0.1), labels = labs2) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, linetype = "solid", colour ="black")) +
  theme(strip.background =  element_rect(fill = "white", colour = "white")) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    axis.ticks = element_blank()
  ) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(.5, "lines")) +
  theme(strip.text.y = element_text( size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(
    legend.title.align=0.5,
    legend.key.height=unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical"
  ) + guides(
    fill = guide_legend(title ="Observed LD in the QTN region", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./figure/LD_QTN_2815.pdf", width = 30, height = 25, units = "cm", dpi=300,family="Times")

#---- summary of distance ------
ld  <- fread("summary_all_LD.txt", data.table = F) %>% 
  mutate(maf = gsub("_.*", "", maf)) %>% group_by(specie, sample_size, ld, maf, type, h2) %>% mutate(rep = order(rep)) %>%
  ungroup()

all_af <- ld %>% select(cols = -c(mt:ld_QTN )) %>%
  mutate(
    h2 = recode(h2, "0.3_0.3" = "h^2:~0.3~-~0.3", "0.3_0.8" = "h^2:~0.3~-~0.8", "0.8_0.8" = "h^2:~0.8~-~0.8"),
    type = recode(type, "direct" = "Direct", "indirect" = "Indirect"),
    specie = recode(specie, "ames" = "Maize", "soy" = "Soybean"),
    sample_size = as.factor(sample_size),
    sample_size = fct_relevel(sample_size, "500", "1000", "2815"),
    ld = as.character(ld),
    maf = ifelse(maf == 0.4, "0.40", "0.05")) %>%
  rename("Specie" = "specie", "LD" = "ld", "MAF" = "maf") %>% 
  mutate(vars = paste0(type, LD),
         vars = recode(vars, "Direct0.1" = "Direct LD: 0.01", "Direct0.99" = "Direct LD: 0.98", "Indirect0.1" = "Indirect LD: 0.01", "Indirect0.99" = "Indirect LD: 0.98"),
         vars = fct_relevel(vars,"Indirect LD: 0.01", "Indirect LD: 0.98", "Direct LD: 0.01", "Direct LD: 0.98"))


all_af %>% filter(h2 == "h^2:~0.3~-~0.8", sample_size == 1000) %>%
  ggplot(aes(x = log(dist_QTN))) +
  geom_histogram(colour = "black", fill = "green") +
  facet_grid(Specie + MAF ~  vars, labeller = labeller(MAF = label_both)) +
  xlab("Distance between QTNs log(bp)") +  ylab("") +  theme_bw(base_size = 20) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = .48, vjust = -4, family = "Times"),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(title ="", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./figure/Distance.pdf", width = 30, height = 25, units = "cm", dpi=300,family="Times")
