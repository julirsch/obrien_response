#' ---
#' title: "Compare studies across gene sets"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' ### Compare GATA1 target gene set across erythroid maturation and for DBA/healthy for each study.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(BuenColors)

#' Read in processed data and a subset of gene sets.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.SCAN.GSE22552 <- readRDS("../processed/combined_SCAN_GSE22552.rds")
combined.SCAN.GSE24759 <- readRDS("../processed/combined_SCAN_GSE24759.rds")
combined.SCAN.GSE41817 <- readRDS("../processed/combined_SCAN_GSE41817.rds")
combined.SCAN.GSE89540 <- readRDS("../processed/combined_SCAN_GSE89540.rds")
combined.SCAN.BN.GSE22552 <- readRDS("../processed/combined_SCAN_BN_GSE22552.rds")
combined.SCAN.BN.GSE24759 <- readRDS("../processed/combined_SCAN_BN_GSE24759.rds")
combined.SCAN.BN.GSE89540 <- readRDS("../processed/combined_SCAN_BN_GSE89540.rds")
combined.RMA.GSE22552 <- readRDS("../processed/combined_RMA_GSE22552.rds")
combined.RMA.GSE24759 <- readRDS("../processed/combined_RMA_GSE24759.rds")
combined.RMA.GSE41817 <- readRDS("../processed/combined_RMA_GSE41817.rds")
combined.RMA.GSE89540 <- readRDS("../processed/combined_RMA_GSE89540.rds")
combined.RMA.BN.GSE22552 <- readRDS("../processed/combined_RMA_BN_GSE22552.rds")
combined.RMA.BN.GSE24759 <- readRDS("../processed/combined_RMA_BN_GSE24759.rds")
combined.RMA.BN.GSE89540 <- readRDS("../processed/combined_RMA_BN_GSE89540.rds")
geneset <- read.table("../data/GeneSets.txt",header=T,sep="\t",na.strings="",stringsAsFactors=F)
geneset.l <-lapply(geneset, function (x) x[!is.na(x)])
geneset1 <- geneset$Transfac_GATA1
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage","group1","group2")

#' GSE22552 RMA
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE22552)"
combined.RMA.GSE22552.temp <- combined.RMA.GSE22552
combined.RMA.GSE22552.melt <- melt(as.matrix(combined.RMA.GSE22552.temp))
combined.RMA.GSE22552.melt <- list(combined.RMA.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.RMA.GSE22552.melt$group1 <- factor(combined.RMA.GSE22552.melt$group1, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
combined.RMA.GSE22552.melt.mean <- combined.RMA.GSE22552.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.GSE22552.melt.mean <- combined.RMA.GSE22552.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.GSE22552.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) + 
  scale_color_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw()
kruskal.test(val ~ group1,data=combined.RMA.GSE22552.melt.mean)

#' GSE24759 RMA
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=10, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE24759)"
combined.RMA.GSE24759.temp <- combined.RMA.GSE24759
combined.RMA.GSE24759.melt <- melt(as.matrix(combined.RMA.GSE24759.temp))
combined.RMA.GSE24759.melt <- list(combined.RMA.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.RMA.GSE24759.melt$group1 <- factor(combined.RMA.GSE24759.melt$group1, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos")), ordered=T)
combined.RMA.GSE24759.melt.mean <- combined.RMA.GSE24759.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.GSE24759.melt.mean <- combined.RMA.GSE24759.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.GSE24759.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) + 
  scale_color_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw()
kruskal.test(val ~ group1,data=combined.RMA.GSE24759.melt.mean)

#' GSE41817 RMA
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.RMA.GSE41817.temp <- combined.RMA.GSE41817
combined.RMA.GSE41817.melt <- melt(as.matrix(combined.RMA.GSE41817.temp))
combined.RMA.GSE41817.melt <- list(combined.RMA.GSE41817.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.RMA.GSE41817.melt.mean <- combined.RMA.GSE41817.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.GSE41817.melt.mean <- combined.RMA.GSE41817.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.GSE41817.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4)]) +
  theme_bw()
kruskal.test(val ~ factor(group1),data=combined.RMA.GSE41817.melt.mean)

#' GSE89540 RMA
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.RMA.GSE89540.temp <- combined.RMA.GSE89540 %>%
  dplyr::select(contains("44"))
combined.RMA.GSE89540.melt <- melt(as.matrix(combined.RMA.GSE89540.temp))
combined.RMA.GSE89540.melt <- list(combined.RMA.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.RMA.GSE89540.melt.mean <- combined.RMA.GSE89540.melt %>%
  filter(group2 == 44) %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.GSE89540.melt.mean <- combined.RMA.GSE89540.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.GSE89540.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) +
  theme_bw()
kruskal.test(val ~ factor(group1),data=combined.RMA.GSE89540.melt.mean)

#' GSE22552 SCAN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE22552)"
combined.SCAN.GSE22552.temp <- combined.SCAN.GSE22552
combined.SCAN.GSE22552.melt <- melt(as.matrix(combined.SCAN.GSE22552.temp))
combined.SCAN.GSE22552.melt <- list(combined.SCAN.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.SCAN.GSE22552.melt$group1 <- factor(combined.SCAN.GSE22552.melt$group1, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
combined.SCAN.GSE22552.melt.mean <- combined.SCAN.GSE22552.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.GSE22552.melt.mean <- combined.SCAN.GSE22552.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.GSE22552.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) + 
  scale_color_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw()
kruskal.test(val ~ group1,data=combined.SCAN.GSE22552.melt.mean)

#' GSE24759 SCAN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=8, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE24759)"
combined.SCAN.GSE24759.temp <- combined.SCAN.GSE24759
combined.SCAN.GSE24759.melt <- melt(as.matrix(combined.SCAN.GSE24759.temp))
combined.SCAN.GSE24759.melt <- list(combined.SCAN.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.SCAN.GSE24759.melt$group1 <- factor(combined.SCAN.GSE24759.melt$group1, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos")), ordered=T)
combined.SCAN.GSE24759.melt.mean <- combined.SCAN.GSE24759.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.GSE24759.melt.mean <- combined.SCAN.GSE24759.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.GSE24759.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) + 
  scale_color_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw()
kruskal.test(val ~ group1,data=combined.SCAN.GSE24759.melt.mean)

#' GSE41817 SCAN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.SCAN.GSE41817.temp <- combined.SCAN.GSE41817
combined.SCAN.GSE41817.melt <- melt(as.matrix(combined.SCAN.GSE41817.temp))
combined.SCAN.GSE41817.melt <- list(combined.SCAN.GSE41817.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.SCAN.GSE41817.melt.mean <- combined.SCAN.GSE41817.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.GSE41817.melt.mean <- combined.SCAN.GSE41817.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.GSE41817.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4)]) +
  theme_bw()
kruskal.test(val ~ factor(group1),data=combined.SCAN.GSE41817.melt.mean)

#' GSE89540 SCAN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.SCAN.GSE89540.temp <- combined.SCAN.GSE89540 %>%
  dplyr::select(contains("44"))
combined.SCAN.GSE89540.melt <- melt(as.matrix(combined.SCAN.GSE89540.temp))
combined.SCAN.GSE89540.melt <- list(combined.SCAN.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.SCAN.GSE89540.melt.mean <- combined.SCAN.GSE89540.melt %>%
  filter(group2 == 44) %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.GSE89540.melt.mean <- combined.SCAN.GSE89540.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.GSE89540.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) +
  theme_bw()
kruskal.test(val ~ factor(group1),data=combined.SCAN.GSE89540.melt.mean)

#' GSE22552 RMA.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE22552)"
combined.RMA.BN.GSE22552.temp <- combined.RMA.BN.GSE22552
combined.RMA.BN.GSE22552.melt <- melt(as.matrix(combined.RMA.BN.GSE22552.temp))
combined.RMA.BN.GSE22552.melt <- list(combined.RMA.BN.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.RMA.BN.GSE22552.melt$group1 <- factor(combined.RMA.BN.GSE22552.melt$group1, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
combined.RMA.BN.GSE22552.melt.mean <- combined.RMA.BN.GSE22552.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.BN.GSE22552.melt.mean <- combined.RMA.BN.GSE22552.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.BN.GSE22552.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) + 
  scale_color_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw()

#' GSE24759 RMA.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=8, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE24759)"
combined.RMA.BN.GSE24759.temp <- combined.RMA.BN.GSE24759
combined.RMA.BN.GSE24759.melt <- melt(as.matrix(combined.RMA.BN.GSE24759.temp))
combined.RMA.BN.GSE24759.melt <- list(combined.RMA.BN.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.RMA.BN.GSE24759.melt$group1 <- factor(combined.RMA.BN.GSE24759.melt$group1, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos")), ordered=T)
combined.RMA.BN.GSE24759.melt.mean <- combined.RMA.BN.GSE24759.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.BN.GSE24759.melt.mean <- combined.RMA.BN.GSE24759.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.BN.GSE24759.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) + 
  scale_color_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw()

#' GSE89540 RMA.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.RMA.BN.GSE89540.temp <- combined.RMA.BN.GSE89540 %>%
  dplyr::select(contains("44"))
combined.RMA.BN.GSE89540.melt <- melt(as.matrix(combined.RMA.BN.GSE89540.temp))
combined.RMA.BN.GSE89540.melt <- list(combined.RMA.BN.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.RMA.BN.GSE89540.melt.mean <- combined.RMA.BN.GSE89540.melt %>%
  filter(group2 == 44) %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.RMA.BN.GSE89540.melt.mean <- combined.RMA.BN.GSE89540.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.RMA.BN.GSE89540.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) +
  theme_bw()

#' GSE22552 SCAN.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE22552)"
combined.SCAN.BN.GSE22552.temp <- combined.SCAN.BN.GSE22552
combined.SCAN.BN.GSE22552.melt <- melt(as.matrix(combined.SCAN.BN.GSE22552.temp))
combined.SCAN.BN.GSE22552.melt <- list(combined.SCAN.BN.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.SCAN.BN.GSE22552.melt$group1 <- factor(combined.SCAN.BN.GSE22552.melt$group1, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
combined.SCAN.BN.GSE22552.melt.mean <- combined.SCAN.BN.GSE22552.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.BN.GSE22552.melt.mean <- combined.SCAN.BN.GSE22552.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.BN.GSE22552.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) + 
  scale_color_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw()

#' GSE24759 SCAN.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=8, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by erythroid maturation (GSE24759)"
combined.SCAN.BN.GSE24759.temp <- combined.SCAN.BN.GSE24759
combined.SCAN.BN.GSE24759.melt <- melt(as.matrix(combined.SCAN.BN.GSE24759.temp))
combined.SCAN.BN.GSE24759.melt <- list(combined.SCAN.BN.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.SCAN.BN.GSE24759.melt$group1 <- factor(combined.SCAN.BN.GSE24759.melt$group1, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos")), ordered=T)
combined.SCAN.BN.GSE24759.melt.mean <- combined.SCAN.BN.GSE24759.melt %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.BN.GSE24759.melt.mean <- combined.SCAN.BN.GSE24759.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.BN.GSE24759.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) + 
  scale_color_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw()

#' GSE89540 SCAN.BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=2.5, fig.align='center', fig.cap = "Gene set expression by DBA genotype (O'Brien et al.)"
combined.SCAN.BN.GSE89540.temp <- combined.SCAN.BN.GSE89540 %>%
  dplyr::select(contains("44"))
combined.SCAN.BN.GSE89540.melt <- melt(as.matrix(combined.SCAN.BN.GSE89540.temp))
combined.SCAN.BN.GSE89540.melt <- list(combined.SCAN.BN.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .) 
combined.SCAN.BN.GSE89540.melt.mean <- combined.SCAN.BN.GSE89540.melt %>%
  filter(group2 == 44) %>%
  group_by(Var1,group1) %>%
  summarize(val=mean(value))
combined.SCAN.BN.GSE89540.melt.mean <- combined.SCAN.BN.GSE89540.melt.mean %>%
  filter(Var1 %in% geneset1)
ggplot(combined.SCAN.BN.GSE89540.melt.mean, aes(x=val, fill = group1, color = group1)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) + 
  scale_color_manual(values = jdb_palette("Rushmore")[c(5,4,3)]) +
  theme_bw()


