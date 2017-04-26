#' ---
#' title: "Analyze GSE89540 with reference data using CIBERSORT"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Check if DBA status and genotype are confounded by maturation stage using CIBERSORT.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(BuenColors)
sem<-function(x) {sd(x)/sqrt(length(x))}

#' Read in results from CIBERSORT, which has an online submission portal (https://cibersort.stanford.edu).
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
#Read in 
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage","group1","group2")
CIBERSORT.GSE22552 <- read.table("../processed/CIBERSORT.Output1.txt",sep="\t",header=T)
CIBERSORT.GSE22552 <- merge(CIBERSORT.GSE22552,samples,by.x="Input.Sample",by.y="name")
CIBERSORT.GSE22552$group2 <- as.factor(CIBERSORT.GSE22552$group2)
CIBERSORT.GSE24759<- read.table("../processed/CIBERSORT.Output2.txt",sep="\t",header=T)
CIBERSORT.GSE24759 <- merge(CIBERSORT.GSE24759,samples,by.x="Input.Sample",by.y="name")
CIBERSORT.GSE24759$group2 <- as.factor(CIBERSORT.GSE24759$group2)

#' Deconvolve CD235a+/CD235a- mixtures from GSE22552 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=4, fig.height=2.5, fig.align='center', fig.cap = "O'Brien et al.CD235a+/CD235a- mixtures from GSE22552 cell types"
CIBERSORT.GSE22552.temp <- CIBERSORT.GSE22552 %>%
  filter(group1 == "Control") %>%
  dplyr::select(CFU_E:PRO_E,group2)
CIBERSORT.GSE22552.melt <- melt(CIBERSORT.GSE22552.temp)
CIBERSORT.GSE22552.melt$variable <- factor(CIBERSORT.GSE22552.melt$variable, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
CIBERSORT.GSE22552.sum <- CIBERSORT.GSE22552.melt %>% 
  group_by(group2,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE22552.plot <- CIBERSORT.GSE22552.sum %>% 
  group_by(group2) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE22552.plot, aes(x=group2, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

#' Deconvolve DBA genotypes for CD235a- mixtures from GSE22552 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=4, fig.height=4, fig.align='center', fig.cap = "O'Brien et al. CD235a- DBA genotype mixtures from GSE22552 cell types"
CIBERSORT.GSE22552.temp <- CIBERSORT.GSE22552 %>%
  filter(group2 == "44") %>%
  dplyr::select(CFU_E:PRO_E,group1)
CIBERSORT.GSE22552.melt <- melt(CIBERSORT.GSE22552.temp)
CIBERSORT.GSE22552.melt$variable <- factor(CIBERSORT.GSE22552.melt$variable, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
CIBERSORT.GSE22552.sum <- CIBERSORT.GSE22552.melt %>% 
  group_by(group1,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE22552.plot <- CIBERSORT.GSE22552.sum %>% 
  group_by(group1) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE22552.plot, aes(x=group1, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

#' Deconvolve DBA genotypes for CD235a+ mixtures from GSE22552 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=4, fig.height=4, fig.align='center', fig.cap = "O'Brien et al. CD235a+ DBA genotype mixtures from GSE22552 cell types"
CIBERSORT.GSE22552.temp <- CIBERSORT.GSE22552 %>%
  filter(group2 == "235") %>%
  dplyr::select(CFU_E:PRO_E,group1)
CIBERSORT.GSE22552.melt <- melt(CIBERSORT.GSE22552.temp)
CIBERSORT.GSE22552.melt$variable <- factor(CIBERSORT.GSE22552.melt$variable, levels= rev(c("CFU_E","PRO_E","INT_E","LATE_E")), ordered=T)
CIBERSORT.GSE22552.sum <- CIBERSORT.GSE22552.melt %>% 
  group_by(group1,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE22552.plot <- CIBERSORT.GSE22552.sum %>% 
  group_by(group1) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE22552.plot, aes(x=group1, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("flame_light")[c(2,4,6,8)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

#' Deconvolve CD235a+/CD235a- mixtures from CIBERSORT.GSE24759 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=4, fig.height=4, fig.align='center', fig.cap = "O'Brien et al.CD235a+/CD235a- mixtures from CIBERSORT.GSE24759 cell types"
CIBERSORT.GSE24759.temp <- CIBERSORT.GSE24759 %>%
  filter(group1 == "Control") %>%
  dplyr::select(CD34neg_CD71lo_GlyApos:CD34pos_CD71pos_GlyAneg,group2)
CIBERSORT.GSE24759.melt <- melt(CIBERSORT.GSE24759.temp)
CIBERSORT.GSE24759.melt$variable <- factor(CIBERSORT.GSE24759.melt$variable, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71neg_GlyApos","CD34neg_CD71lo_GlyApos")), ordered=T)
CIBERSORT.GSE24759.sum <- CIBERSORT.GSE24759.melt %>% 
  group_by(group2,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE24759.plot <- CIBERSORT.GSE24759.sum %>% 
  group_by(group2) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE24759.plot, aes(x=group2, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

#' Deconvolve DBA genotypes for CD235a- mixtures from GSE24759 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=5, fig.height=4, fig.align='center', fig.cap = "O'Brien et al. CD235a- DBA genotype mixtures from GSE24759 cell types"
CIBERSORT.GSE24759.temp <- CIBERSORT.GSE24759 %>%
  filter(group2 == "44") %>%
  dplyr::select(CD34neg_CD71lo_GlyApos:CD34pos_CD71pos_GlyAneg,group1)
CIBERSORT.GSE24759.melt <- melt(CIBERSORT.GSE24759.temp)
CIBERSORT.GSE24759.melt$variable <- factor(CIBERSORT.GSE24759.melt$variable, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71neg_GlyApos","CD34neg_CD71lo_GlyApos")), ordered=T)
CIBERSORT.GSE24759.sum <- CIBERSORT.GSE24759.melt %>% 
  group_by(group1,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE24759.plot <- CIBERSORT.GSE24759.sum %>% 
  group_by(group1) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE24759.plot, aes(x=group1, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

#' Deconvolve DBA genotypes for CD235a+ mixtures from GSE24759 cell types.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=5, fig.height=4, fig.align='center', fig.cap = "O'Brien et al. CD235a+ DBA genotype mixtures from GSE24759 cell types"
CIBERSORT.GSE24759.temp <- CIBERSORT.GSE24759 %>%
  filter(group2 == "235") %>%
  dplyr::select(CD34neg_CD71lo_GlyApos:CD34pos_CD71pos_GlyAneg,group1)
CIBERSORT.GSE24759.melt <- melt(CIBERSORT.GSE24759.temp)
CIBERSORT.GSE24759.melt$variable <- factor(CIBERSORT.GSE24759.melt$variable, levels= rev(c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71neg_GlyApos","CD34neg_CD71lo_GlyApos")), ordered=T)
CIBERSORT.GSE24759.sum <- CIBERSORT.GSE24759.melt %>% 
  group_by(group1,variable) %>%
  summarize(percent = mean(value)) %>%
  arrange(desc(variable))
CIBERSORT.GSE24759.plot <- CIBERSORT.GSE24759.sum %>% 
  group_by(group1) %>%
  mutate(pos = cumsum(percent)) %>%
  ungroup() %>%
  mutate(lower = pos - sem(percent), upper = pos + sem(percent))
ggplot(CIBERSORT.GSE24759.plot, aes(x=group1, y=percent, fill=variable)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.2, width = .1, col = "black") +
  scale_fill_manual(values = jdb_palette("brewer_celsius")[c(9,7,5,3,1)]) +
  theme_bw() +
  xlab("") + 
  ylab ("")

