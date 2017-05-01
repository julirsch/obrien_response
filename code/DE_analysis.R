#' ---
#' title: "Analyze GSE89540 with reference data using SCAN (gene sets)"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Investigate the expression of key genes potentially involed in DBA pathogenesis
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries
library(oligo)
library(limma)
library(dplyr)
library(reshape2)

#' Read in processed data and a subset of gene sets.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.SCAN.GSE22552 <- readRDS("../processed/combined_SCAN_GSE22552.rds")
combined.SCAN.GSE24759 <- readRDS("../processed/combined_SCAN_GSE24759.rds")
combined.SCAN.GSE41817 <- readRDS("../processed/combined_SCAN_GSE41817.rds")
combined.SCAN.GSE89540 <- readRDS("../processed/combined_SCAN_GSE89540.rds")
combined.RMA.GSE22552 <- readRDS("../processed/combined_RMA_GSE22552.rds")
combined.RMA.GSE24759 <- readRDS("../processed/combined_RMA_GSE24759.rds")
combined.RMA.GSE41817 <- readRDS("../processed/combined_RMA_GSE41817.rds")
combined.RMA.GSE89540 <- readRDS("../processed/combined_RMA_GSE89540.rds")
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage","group1","group2")

#' Perform differential expression analysis for SCAN normalized gene sets
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.SCAN.GSE41817.melt <- melt(as.matrix(combined.SCAN.GSE41817))
combined.SCAN.GSE41817.melt <- list(combined.SCAN.GSE41817.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.SCAN.GSE41817)) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA")
cont_matrix <- makeContrasts(DBA - Control,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.SCAN.GSE41817))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.SCAN.GSE41817.table <- topTable(fit, n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.SCAN.GSE41817.table,"../processed/combined.SCAN.GSE41817.table.rds")

# GSE89540
combined.SCAN.GSE89540.temp <- combined.SCAN.GSE89540 %>%
  dplyr::select(contains("44"))
combined.SCAN.GSE89540.melt <- melt(as.matrix(combined.SCAN.GSE89540.temp))
combined.SCAN.GSE89540.melt <- list(combined.SCAN.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.SCAN.GSE89540.temp)) %>%
  filter(group2 == 44) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA","DBA_GATA1"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA","DBA_GATA1")
cont_matrix <- makeContrasts(DBA - Control,
                             DBA_GATA1 - Control,
                             DBA_GATA1 - DBA,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.SCAN.GSE89540.temp))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.SCAN.GSE89540.table.DBAvControl <- topTable(fit, coef="DBA - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.SCAN.GSE89540.table.DBA_GATA1vControl <- topTable(fit, coef="DBA_GATA1 - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.SCAN.GSE89540.table.DBA_GATA1vDBA <- topTable(fit, coef="DBA_GATA1 - DBA",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.SCAN.GSE89540.table.DBAvControl,"../processed/combined.SCAN.GSE89540.table.DBAvControl.rds")
saveRDS(combined.SCAN.GSE89540.table.DBA_GATA1vControl,"../processed/combined.SCAN.GSE89540.table.DBA_GATA1vControl.rds")
saveRDS(combined.SCAN.GSE89540.table.DBA_GATA1vDBA,"../processed/combined.SCAN.GSE89540.table.DBA_GATA1vDBA.rds")

#' Perform differential expression analysis for RMA normalized gene sets
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.RMA.GSE41817.melt <- melt(as.matrix(combined.RMA.GSE41817))
combined.RMA.GSE41817.melt <- list(combined.RMA.GSE41817.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.RMA.GSE41817)) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA")
cont_matrix <- makeContrasts(DBA - Control,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.RMA.GSE41817))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.RMA.GSE41817.table <- topTable(fit, n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.RMA.GSE41817.table,"../processed/combined.RMA.GSE41817.table.rds")

# GSE89540
combined.RMA.GSE89540.temp <- combined.RMA.GSE89540 %>%
  dplyr::select(contains("44"))
combined.RMA.GSE89540.melt <- melt(as.matrix(combined.RMA.GSE89540.temp))
combined.RMA.GSE89540.melt <- list(combined.RMA.GSE89540.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.RMA.GSE89540.temp)) %>%
  filter(group2 == 44) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA","DBA_GATA1"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA","DBA_GATA1")
cont_matrix <- makeContrasts(DBA - Control,
                             DBA_GATA1 - Control,
                             DBA_GATA1 - DBA,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.RMA.GSE89540.temp))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.RMA.GSE89540.table.DBAvControl <- topTable(fit, coef="DBA - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.table.DBA_GATA1vControl <- topTable(fit, coef="DBA_GATA1 - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.table.DBA_GATA1vDBA <- topTable(fit, coef="DBA_GATA1 - DBA",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.RMA.GSE89540.table.DBAvControl,"../processed/combined.RMA.GSE89540.table.DBAvControl.rds")
saveRDS(combined.RMA.GSE89540.table.DBA_GATA1vControl,"../processed/combined.RMA.GSE89540.table.DBA_GATA1vControl.rds")
saveRDS(combined.RMA.GSE89540.table.DBA_GATA1vDBA,"../processed/combined.RMA.GSE89540.table.DBA_GATA1vDBA.rds")


#' Check for differential expression of key genes.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
combined.RMA.GSE41817.table["ALAS2",]
combined.RMA.GSE41817.table["GATA1",]
combined.RMA.GSE41817.table["HBB",]
combined.RMA.GSE89540.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.table.DBAvControl["HBB",]
combined.RMA.GSE89540.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.table.DBAvControl["HBB",]
combined.RMA.GSE89540.table.DBA_GATA1vControl["ALAS2",]
combined.RMA.GSE89540.table.DBA_GATA1vControl["GATA1",]
combined.RMA.GSE89540.table.DBA_GATA1vControl["HBB",]
combined.RMA.GSE89540.table.DBA_GATA1vDBA["ALAS2",]
combined.RMA.GSE89540.table.DBA_GATA1vDBA["GATA1",]
combined.RMA.GSE89540.table.DBA_GATA1vDBA["HBB",]
combined.SCAN.GSE41817.table["ALAS2",]
combined.SCAN.GSE41817.table["GATA1",]
combined.SCAN.GSE41817.table["HBB",]
combined.SCAN.GSE89540.table.DBAvControl["ALAS2",]
combined.SCAN.GSE89540.table.DBAvControl["GATA1",]
combined.SCAN.GSE89540.table.DBAvControl["HBB",]
combined.SCAN.GSE89540.table.DBAvControl["ALAS2",]
combined.SCAN.GSE89540.table.DBAvControl["GATA1",]
combined.SCAN.GSE89540.table.DBAvControl["HBB",]
combined.SCAN.GSE89540.table.DBA_GATA1vControl["ALAS2",]
combined.SCAN.GSE89540.table.DBA_GATA1vControl["GATA1",]
combined.SCAN.GSE89540.table.DBA_GATA1vControl["HBB",]
combined.SCAN.GSE89540.table.DBA_GATA1vDBA["ALAS2",]
combined.SCAN.GSE89540.table.DBA_GATA1vDBA["GATA1",]
combined.SCAN.GSE89540.table.DBA_GATA1vDBA["HBB",]

#' Create synthetic normal and normalize O'Brien et al. expression set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.GSE89540 <- readRDS("../processed/combined_RMA_GSE89540.rds")
CIBERSORT.GSE22552 <- readRDS("../processed/CIBERSORT.GSE22552.mixture.rds")
combined.GSE89540 <- combined.GSE89540 %>%
  dplyr::select(grep("44",names(.)))
combined.RMA.GSE89540.sn <- combined.GSE89540
CFU_E <- combined.GSE22552 %>%
  dplyr::select(grep("CFU_E",names(.))) %>%
  rowMeans()
PRO_E <- combined.GSE22552 %>%
  dplyr::select(grep("Pro_E",names(.))) %>%
  rowMeans()
INT_E <- combined.GSE22552 %>%
  dplyr::select(grep("Int_E",names(.))) %>%
  rowMeans()
LATE_E <- combined.GSE22552 %>%
  dplyr::select(grep("Late_E",names(.))) %>%
  rowMeans()
seqto <- dim(CIBERSORT.GSE22552)[1]
for (i in seq(1, seqto, 1)) {
  combined.RMA.GSE89540.sn[,as.vector(CIBERSORT.GSE22552[i,]$Row.names)] <-
    CIBERSORT.GSE22552[i,"CFU_E"] * CFU_E +
    CIBERSORT.GSE22552[i,"PRO_E"] *  PRO_E +
    CIBERSORT.GSE22552[i,"INT_E"] *  INT_E +
    CIBERSORT.GSE22552[i,"LATE_E"] * LATE_E
}
combined.RMA.GSE89540.sn.norm <- combined.GSE89540 - combined.GSE89540.sn

# GSE89540 synthetic normal 
combined.RMA.GSE89540.sn.temp <- combined.RMA.GSE89540.sn
combined.RMA.GSE89540.sn.melt <- melt(as.matrix(combined.RMA.GSE89540.sn.temp))
combined.RMA.GSE89540.sn.melt <- list(combined.RMA.GSE89540.sn.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.RMA.GSE89540.sn.temp)) %>%
  filter(group2 == 44) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA","DBA_GATA1"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA","DBA_GATA1")
cont_matrix <- makeContrasts(DBA - Control,
                             DBA_GATA1 - Control,
                             DBA_GATA1 - DBA,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.RMA.GSE89540.sn.temp))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.RMA.GSE89540.sn.table.DBAvControl <- topTable(fit, coef="DBA - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl <- topTable(fit, coef="DBA_GATA1 - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA <- topTable(fit, coef="DBA_GATA1 - DBA",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.RMA.GSE89540.sn.table.DBAvControl,"../processed/combined.RMA.GSE89540.sn.table.DBAvControl.rds")
saveRDS(combined.RMA.GSE89540.sn.table.DBA_GATA1vControl,"../processed/combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.rds")
saveRDS(combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA,"../processed/combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA.rds")

# GSE89540 normalized (residuals after removing synthetic normal)
combined.RMA.GSE89540.sn.norm.temp <- combined.RMA.GSE89540.sn.norm
combined.RMA.GSE89540.sn.norm.melt <- melt(as.matrix(combined.RMA.GSE89540.sn.norm.temp))
combined.RMA.GSE89540.sn.norm.melt <- list(combined.RMA.GSE89540.sn.norm.melt, samples) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("Var2"="name")), .)
# Make one vs all comparisons for each
groups <- samples %>%
  filter(name %in% names(combined.RMA.GSE89540.sn.norm.temp)) %>%
  filter(group2 == 44) %>%
  .$group1 %>%
  factor(levels=c("Control","DBA","DBA_GATA1"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("Control","DBA","DBA_GATA1")
cont_matrix <- makeContrasts(DBA - Control,
                             DBA_GATA1 - Control,
                             DBA_GATA1 - DBA,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.RMA.GSE89540.sn.norm.temp))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
combined.RMA.GSE89540.sn.norm.table.DBAvControl <- topTable(fit, coef="DBA - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl <- topTable(fit, coef="DBA_GATA1 - Control",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA <- topTable(fit, coef="DBA_GATA1 - DBA",n = 20000, sort.by="logFC", resort.by="logFC", p.value=1,adjust.method="fdr",lfc=0)
saveRDS(combined.RMA.GSE89540.sn.norm.table.DBAvControl,"../processed/combined.RMA.GSE89540.sn.norm.table.DBAvControl.rds")
saveRDS(combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl,"../processed/combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.rds")
saveRDS(combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA,"../processed/combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA.rds")

#' Check for differential expression of key genes in synthetic normal and after removing synthetic normal from O'Brien et al. expression set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
#Synthetic normal
combined.RMA.GSE89540.sn.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.sn.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.sn.table.DBAvControl["HBB",]
combined.RMA.GSE89540.sn.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.sn.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.sn.table.DBAvControl["HBB",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl["ALAS2",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl["GATA1",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl["HBB",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA["ALAS2",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA["GATA1",]
combined.RMA.GSE89540.sn.table.DBA_GATA1vDBA["HBB",]
#Normalized
combined.RMA.GSE89540.sn.norm.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.sn.norm.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.sn.norm.table.DBAvControl["HBB",]
combined.RMA.GSE89540.sn.norm.table.DBAvControl["ALAS2",]
combined.RMA.GSE89540.sn.norm.table.DBAvControl["GATA1",]
combined.RMA.GSE89540.sn.norm.table.DBAvControl["HBB",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl["ALAS2",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl["GATA1",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl["HBB",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA["ALAS2",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA["GATA1",]
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vDBA["HBB",]

