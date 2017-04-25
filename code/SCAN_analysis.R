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

#' Read in data for each experiment on GEO and save as RDS object so we only have to do this once
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE, eval = FALSE
SCAN.GSE22552 <- SCAN("../data/GSE22552_RAW/*")
saveRDS(SCAN.GSE22552,"../processed/SCAN.GSE22552.rds")
SCAN.GSE24759 <- SCAN("../data/GSE24759_RAW/*", annotationPackageName="pd.ht.hg.u133a")
saveRDS(SCAN.GSE24759,"../processed/SCAN.GSE24759.rds")
SCAN.GSE34268 <- SCAN("../data/GSE34268_RAW/*")
saveRDS(SCAN.GSE34268,"../processed/SCAN.GSE34268.rds")
SCAN.GSE41599 <- SCAN("../data/GSE41599_RAW/*")
saveRDS(SCAN.GSE41599,"../processed/SCAN.GSE41599.rds")
SCAN.GSE41817 <- SCAN("../data/GSE41817_RAW/*")
saveRDS(SCAN.GSE41817,"../processed/SCAN.GSE41817.rds")
SCAN.GSE89540 <- SCAN("../data/GSE89540_RAW/*")
saveRDS(SCAN.GSE89540,"../processed/SCAN.GSE89540.rds")

#' Load in SCAN-normalized expression
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = FALSE, eval = TRUE
SCAN.GSE22552 <- readRDS("../processed/SCAN.GSE22552.rds")
SCAN.GSE24759 <- readRDS("../processed/SCAN.GSE24759.rds")
SCAN.GSE34268 <- readRDS("../processed/SCAN.GSE34268.rds")
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
# GSE34268
library(hgu133a.db)
annotation(SCAN.GSE34268) <- "hgu133a.db"
ID <- featureNames(SCAN.GSE34268)
Symbol <- getSYMBOL(ID,"hgu133a.db")
fData(SCAN.GSE34268) <- data.frame(ID=ID,Symbol=Symbol)
exprs.GSE34268 <- as.data.frame(exprs(SCAN.GSE34268))
exprs.GSE34268$gene <- as.vector(fData(SCAN.GSE34268)$Symbol)
exprs.GSE34268 <- exprs.GSE34268 %>%
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
#combined <- list(exprs.GSE22552, exprs.GSE24759, exprs.GSE34268, exprs.GSE41599, exprs.GSE41817, exprs.GSE89540) %>%
combined <- list(exprs.GSE22552, exprs.GSE24759, exprs.GSE89540) %>%
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

#' Batch normalize! Can choose model with only intercept or include covariate for a rough sorted "stage"
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
modcombat <- model.matrix(~1, data=combined.t)
#modcombat <- model.matrix(~1 + samples$stage, data=combined.t)
combined.bn <- ComBat(dat=combined, batch=as.factor(samples$GSE), mod=modcombat, par.prior=TRUE, prior.plots=F)
row.names(combined.bn) <- row.names(combined)
combined.bn.t <- data.frame(t(combined.bn))

#' Done with pre-processing. Save both normalized and normalized + batch normalized data.
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
combined <- 2^combined
combined.bn <- 2^combined.bn
saveRDS(combined,"../processed/combined_SCAN.rds")
write.table(combined[,samples[samples$GSE %in% "GSE22552",]$name][,-c(13:16)],"../processed/combined_SCAN_GSE22552.txt",quote=F,sep="\t",col.names=NA)
write.table(combined[,samples[samples$GSE %in% "GSE24759",]$name],"../processed/combined_SCAN_GSE24759.txt",quote=F,sep="\t",col.names=NA)
write.table(combined[,samples[samples$GSE %in% "GSE34268",]$name],"../processed/combined_SCAN_GSE34268.txt",quote=F,sep="\t",col.names=NA)
write.table(combined[,samples[samples$GSE %in% "GSE41599",]$name],"../processed/combined_SCAN_GSE41599.txt",quote=F,sep="\t",col.names=NA)
write.table(combined[,samples[samples$GSE %in% "GSE41817",]$name],"../processed/combined_SCAN_GSE41817.txt",quote=F,sep="\t",col.names=NA)
write.table(combined[,samples[samples$GSE %in% "GSE89540",]$name],"../processed/combined_SCAN_GSE89540.txt",quote=F,sep="\t",col.names=NA)
saveRDS(combined.bn,"../processed/combined_SCAN_BN.rds")
write.table(combined.bn[,samples[samples$GSE %in% "GSE22552",]$name][,-c(13:16)],"../processed/combined_SCAN_BN_GSE22552.txt",quote=F,sep="\t",col.names=NA)
write.table(combined.bn[,samples[samples$GSE %in% "GSE24759",]$name],"../processed/combined_SCAN_BN_GSE24759.txt",quote=F,sep="\t",col.names=NA)
write.table(combined.bn[,samples[samples$GSE %in% "GSE34268",]$name],"../processed/combined_SCAN_BN_GSE34268.txt",quote=F,sep="\t",col.names=NA)
write.table(combined.bn[,samples[samples$GSE %in% "GSE41599",]$name],"../processed/combined_SCAN_BN_GSE41599.txt",quote=F,sep="\t",col.names=NA)
write.table(combined.bn[,samples[samples$GSE %in% "GSE41817",]$name],"../processed/combined_SCAN_BN_GSE41817.txt",quote=F,sep="\t",col.names=NA)
write.table(combined.bn[,samples[samples$GSE %in% "GSE89540",]$name],"../processed/combined_SCAN_BN_GSE89540.txt",quote=F,sep="\t",col.names=NA)



keep <- samples %>% 
  filter(GSE == "GSE22552") %>% 
  filter(group1 != "") %>% 
  .$name 
combined.bn.GSE22552 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
combined.bn.GSE22552.melt <- melt(as.matrix(combined.bn.GSE22552))
combined.bn.GSE22552.melt <- list(combined.bn.GSE22552.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.bn.GSE22552.aov <- combined.bn.GSE22552.melt %>% 
  group_by(Var1) %>% 
  do({model = summary(lm(value ~ group1, data = .));
    data.frame(fval = model$fstatistic[1])}) %>%
  arrange(desc(fval)) %>%
  head(1000)


keep <- samples %>% 
  filter(GSE == "GSE24759") %>% 
  filter(group1 != "") %>% 
  .$name
combined.bn.GSE24759 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
combined.bn.GSE24759.melt <- melt(as.matrix(combined.bn.GSE24759))
combined.bn.GSE24759.melt <- list(combined.bn.GSE24759.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.bn.GSE24759.aov <- combined.bn.GSE24759.melt %>% 
  group_by(Var1) %>% 
  do({model = summary(lm(value ~ group1, data = .));
  data.frame(fval = model$fstatistic[1])}) %>%
  arrange(desc(fval)) %>%
  head(1000)


keep <- samples %>% 
  filter(GSE == "GSE34268") %>% 
  filter(group1 != "") %>% 
  .$name
combined.bn.GSE34268 <- data.frame(combined.bn) %>% 
  dplyr::select(one_of(keep))
combined.bn.GSE34268.melt <- melt(as.matrix(combined.bn.GSE34268))
combined.bn.GSE34268.melt <- list(combined.bn.GSE34268.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
combined.bn.GSE34268.aov <- combined.bn.GSE34268.melt %>% 
  group_by(Var1) %>% 
  do({model = summary(lm(value ~ group1, data = .));
  data.frame(fval = model$fstatistic[1])}) %>%
  arrange(desc(fval)) %>%
  head(1000)


combined.bn.GSE34268.out <- combined.bn.GSE34268.melt %>%
  filter(Var1 %in% combined.bn.GSE34268.aov$Var1) %>%
  group_by(group1,Var1) %>%
  summarize(mean = mean(value)) %>%
  spread(group1,mean)

combined.bn.GSE22552.out <- combined.bn.GSE22552.melt %>%
  filter(Var1 %in% names(include)) %>%
  group_by(group1,Var1) %>%
  summarize(mean = mean(value)) %>%
  spread(group1,mean)

combined.bn.GSE24759.out <- combined.bn.GSE24759.melt %>%
  filter(Var1 %in% names(include)) %>%
  group_by(group1,Var1) %>%
  summarize(mean = mean(value)) %>%
  spread(group1,mean)



write.table(combined.bn.GSE24759.out,"../processed/combined_SCAN_BN_SIG_GSE24759.txt",quote=F,sep="\t",row.names=F)


