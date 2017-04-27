#' ---
#' title: "Analyze GSE89540 with reference data using SCAN (pre-processing)"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' SCAN normalizes each sample separately and has been shown to be *at least* as robust as RMA. RMA is normally our choice for analyzing singular datasets, but since these are across multiple platforms and there is batch normalization, it probably makes more sense to SCAN + BN rather than RMA (in groups) + BN as BN may do some strange things between RMA-normalized groups. (Note, we considered fRMA but fRMA reference datasets were not available for HuGene 2.0 ST microarrays).
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries
library(SCAN.UPC)
library(dplyr)
library(annotate)
library(sva)
library(reshape2)
library(tidyr)

#' Read in data for each experiment on GEO and save as RDS object so we only have to do this once
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE, eval = FALSE
SCAN.GSE22552 <- SCAN("../data/GSE22552_RAW/*")
saveRDS(SCAN.GSE22552,"../processed/SCAN.GSE22552.rds")
SCAN.GSE24759 <- SCAN("../data/GSE24759_RAW/*", annotationPackageName="pd.ht.hg.u133a")
saveRDS(SCAN.GSE24759,"../processed/SCAN.GSE24759.rds")
SCAN.GSE34268 <- SCAN("../data/GSE34268_RAW/*")
saveRDS(SCAN.GSE41599,"../processed/SCAN.GSE41599.rds")
SCAN.GSE41817 <- SCAN("../data/GSE41817_RAW/*")
saveRDS(SCAN.GSE41817,"../processed/SCAN.GSE41817.rds")
SCAN.GSE89540 <- SCAN("../data/GSE89540_RAW/*")
saveRDS(SCAN.GSE89540,"../processed/SCAN.GSE89540.rds")

#' Load in SCAN-normalized expression
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = FALSE, eval = TRUE
SCAN.GSE22552 <- readRDS("../processed/SCAN.GSE22552.rds")
SCAN.GSE24759 <- readRDS("../processed/SCAN.GSE24759.rds")
SCAN.GSE41599 <- readRDS("../processed/SCAN.GSE41599.rds")
SCAN.GSE41817 <- readRDS("../processed/SCAN.GSE41817.rds")
SCAN.GSE89540 <- readRDS("../processed/SCAN.GSE89540.rds")

#' Load and add gene level annotation for each microarray study, remove probesets without gene annotations, take max across probesets 
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
# GSE22552
library(hgu133plus2.db)
annotation(SCAN.GSE22552) <- "hgu133plus2.db"
ID <- featureNames(SCAN.GSE22552)
Symbol <- getSYMBOL(ID,"hgu133plus2.db")
fData(SCAN.GSE22552) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE22552 <- as.data.frame(exprs(SCAN.GSE22552))
exprs.GSE22552$gene <- as.vector(fData(SCAN.GSE22552)$Symbol)
exprs.GSE22552 <- exprs.GSE22552 %>%
  filter(gene != "<NA>") %>%
  group_by(gene) %>%
  summarise_each(funs(max))
# GSE24759
library(hthgu133a.db)
annotation(SCAN.GSE24759) <- "hthgu133a.db"
ID <- featureNames(SCAN.GSE24759)
Symbol <- getSYMBOL(ID,"hthgu133a.db")
fData(SCAN.GSE24759) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE24759 <- as.data.frame(exprs(SCAN.GSE24759))
exprs.GSE24759$gene <- as.vector(fData(SCAN.GSE24759)$Symbol)
exprs.GSE24759 <- exprs.GSE24759 %>%
  filter(gene != "<NA>") %>%
  group_by(gene) %>%
  summarise_each(funs(max))
# GSE41599
library(hgu133a.db)
annotation(SCAN.GSE41599) <- "hgu133a.db"
ID <- featureNames(SCAN.GSE41599)
Symbol <- getSYMBOL(ID,"hgu133a.db")
fData(SCAN.GSE41599) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE41599 <- as.data.frame(exprs(SCAN.GSE41599))
exprs.GSE41599$gene <- as.vector(fData(SCAN.GSE41599)$Symbol)
exprs.GSE41599 <- exprs.GSE41599 %>%
  filter(gene != "<NA>") %>%
  group_by(gene) %>%
  summarise_each(funs(max))
# GSE41817
library(hgu133a.db)
annotation(SCAN.GSE41817) <- "hgu133a.db"
ID <- featureNames(SCAN.GSE41817)
Symbol <- getSYMBOL(ID,"hgu133a.db")
fData(SCAN.GSE41817) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE41817 <- as.data.frame(exprs(SCAN.GSE41817))
exprs.GSE41817$gene <- as.vector(fData(SCAN.GSE41817)$Symbol)
exprs.GSE41817 <- exprs.GSE41817 %>%
  filter(gene != "<NA>") %>%
  group_by(gene) %>%
  summarise_each(funs(max))
# GSE89540
library(hugene20sttranscriptcluster.db)
annotation(SCAN.GSE89540) <- "hugene20sttranscriptcluster.db"
ID <- featureNames(SCAN.GSE89540)
Symbol <- getSYMBOL(ID,"hugene20sttranscriptcluster.db")
fData(SCAN.GSE89540) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE89540 <- as.data.frame(exprs(SCAN.GSE89540))
exprs.GSE89540$gene <- as.vector(fData(SCAN.GSE89540)$Symbol)
exprs.GSE89540 <- exprs.GSE89540 %>%
  filter(gene != "<NA>") %>%
  group_by(gene) %>%
  summarise_each(funs(max))

#' Merge datasets together by gene name and convert from GSM to human interpretable sample names
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
combined <- list(exprs.GSE22552, exprs.GSE24759, exprs.GSE41599, exprs.GSE41817, exprs.GSE89540) %>%
  Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by="gene"), .)
names(combined) <- gsub(".CEL.gz","",gsub("_.*","",names(combined)))
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage","group1","group2")
samples$GSM <- factor(samples$GSM,levels=names(combined),ordered=T)
samples <- samples[order(samples$GSM),]
samples <- samples %>% filter(samples$GSM %in% names(combined))
names(combined) <- c("gene",samples$name)
combined <- data.frame(combined)
row.names(combined) <- combined$gene
combined <- combined[,-1]
combined.t <- data.frame(t(combined))
names(combined.t) <- row.names(combined)

#' Batch normalize! Can choose model with only intercept or include covariate for a rough sorted "stage". Really should only include studies with multiple stages (e.g. two references and O'Brien et al) otherwise completely confounded.
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
keep <- samples %>% 
  filter(GSE %in% c("GSE22552", "GSE24759", "GSE89540")) %>% 
  filter(group1 != "")
combined.temp <- data.frame(combined) %>% 
  dplyr::select(one_of(keep$name))
combined.temp.t <- data.frame(t(combined.temp))
names(combined.temp.t) <- row.names(combined.temp)
modcombat <- model.matrix(~1, data=combined.temp.t)
#modcombat <- model.matrix(~1 + samples$stage, data=combined.t)
combined.bn <- ComBat(dat=combined.temp, batch=as.factor(keep$GSE), mod=modcombat, par.prior=TRUE, prior.plots=F)
row.names(combined.bn) <- row.names(combined)
combined.bn.t <- data.frame(t(combined.bn))

#' Split both normalized and normalized + batch normalized into data frames for each study.
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
# GSE22552 BN
keep <- samples %>% 
  filter(GSE == "GSE22552") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.bn.GSE22552 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
# GSE24759 BN
keep <- samples %>% 
  filter(GSE == "GSE24759") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.bn.GSE24759 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
# GSE89540 BN
keep <- samples %>% 
  filter(GSE == "GSE89540") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.bn.GSE89540 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
# GSE22552
keep <- samples %>% 
  filter(GSE == "GSE22552") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.GSE22552 <- data.frame(combined) %>% 
  dplyr::select(one_of(keep))
# GSE24759
keep <- samples %>% 
  filter(GSE == "GSE24759") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.GSE24759 <- data.frame(combined) %>% 
  dplyr::select(one_of(keep))
# GSE41599
keep <- samples %>% 
  filter(GSE == "GSE41599") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.GSE41599 <- data.frame(combined) %>% 
  dplyr::select(one_of(keep))
# GSE41817
keep <- samples %>% 
  filter(GSE == "GSE41817") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.GSE41817 <- data.frame(combined) %>% 
  dplyr::select(one_of(keep))
# GSE89540
keep <- samples %>% 
  filter(GSE == "GSE89540") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.GSE89540 <- data.frame(combined) %>% 
  dplyr::select(one_of(keep))

#' Done with pre-processing. Save both normalized and normalized + batch normalized data. Un-log for CIBERSORT.
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
saveRDS(combined,"../processed/combined_SCAN.rds")
saveRDS(combined.GSE22552,"../processed/combined_SCAN_GSE22552.rds")
saveRDS(combined.GSE24759,"../processed/combined_SCAN_GSE24759.rds")
saveRDS(combined.GSE41599,"../processed/combined_SCAN_GSE41599.rds")
saveRDS(combined.GSE41817,"../processed/combined_SCAN_GSE41817.rds")
saveRDS(combined.GSE89540,"../processed/combined_SCAN_GSE89540.rds")
saveRDS(combined.bn,"../processed/combined_SCAN_BN.rds")
saveRDS(combined.bn.GSE22552,"../processed/combined_SCAN_BN_GSE22552.rds")
saveRDS(combined.bn.GSE24759,"../processed/combined_SCAN_BN_GSE24759.rds")
saveRDS(combined.bn.GSE89540,"../processed/combined_SCAN_BN_GSE89540.rds")

#' Identify most variable genes between groups in the two reference datasets. Write signature gene sets for overlaps.
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
# Functions for doing one vs all t-tests
t.test.onevsall <- function(df,group) {
  df.ret <- df %>% 
    mutate(group3 = ifelse(group1 == group,1,0)) %>%
    group_by(Var1) %>%
    do({model = t.test(value ~ group3, data = .);
    data.frame(pval = model$p.value)}) %>%
    arrange(pval) %>%
    head(100)
  return(df.ret)
}
t.test.multi <- function(df) {
  ret <- list()
  for (i in levels(factor(df$group1))) {
    ret[[i]]<- t.test.onevsall(df,i)
  }
  return(ret)
}
# GSE22552
combined.bn.GSE22552.melt <- melt(as.matrix(combined.bn.GSE22552))
combined.bn.GSE22552.melt <- list(combined.bn.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.bn.GSE22552.aov <- combined.bn.GSE22552.melt %>% 
  group_by(Var1) %>% 
  do({model = summary(lm(value ~ group1, data = .));
    data.frame(fval = model$fstatistic[1])}) %>%
  arrange(desc(fval))
combined.bn.GSE22552.t <- t.test.multi(combined.bn.GSE22552.melt)
combined.bn.GSE22552.t.vect <- as.vector(unique(bind_rows(combined.bn.GSE22552.t)$Var1))
# GSE24759
combined.bn.GSE24759.melt <- melt(as.matrix(combined.bn.GSE24759))
combined.bn.GSE24759.melt <- list(combined.bn.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.bn.GSE24759.aov <- combined.bn.GSE24759.melt %>% 
  group_by(Var1) %>% 
  do({model = summary(lm(value ~ group1, data = .));
  data.frame(fval = model$fstatistic[1])}) %>%
  arrange(desc(fval))
combined.bn.GSE24759.t <- t.test.multi(combined.bn.GSE24759.melt)
combined.bn.GSE24759.t.vect <- as.vector(unique(bind_rows(combined.bn.GSE24759.t)$Var1))

# Find overlaps
#include <- table(c(as.vector(combined.bn.GSE22552.aov$Var1),as.vector(combined.bn.GSE24759.aov$Var1)))
#include <- names(include[include > 1])
#include <- unique(c(include,as.vector(head(combined.bn.GSE22552.aov$Var1,100)),as.vector(head(combined.bn.GSE24759.aov$Var1,100))))
# Create and write GSE22552 signature matrix
combined.bn.GSE22552.out <- combined.bn.GSE22552.melt %>%
  filter(Var1 %in% combined.bn.GSE22552.t.vect) %>%
  group_by(group1,Var1) %>%
  summarize(mean = mean(value)) %>%
  spread(group1,mean)
write.table(combined.bn.GSE22552.out,"../processed/combined_SCAN_BN_SIG_GSE22552.txt",quote=F,sep="\t",row.names=F)
# Create and write GSE24759 signature matrix
combined.bn.GSE24759.out <- combined.bn.GSE24759.melt %>%
  filter(Var1 %in% combined.bn.GSE24759.t.vect) %>%
  group_by(group1,Var1) %>%
  summarize(mean = mean(value)) %>%
  spread(group1,mean)
write.table(combined.bn.GSE24759.out,"../processed/combined_SCAN_BN_SIG_GSE24759.txt",quote=F,sep="\t",row.names=F)






