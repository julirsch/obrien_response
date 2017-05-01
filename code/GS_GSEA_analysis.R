#' ---
#' title: "Create gene sets and run GSEA for DBA by genotype from O'Brien et al. and Ludwig et al."
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Import libraries
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
library(dplyr)
library(reshape2)
library(tidyr)
library(oligo)
library(limma)
library(fgsea)

#' Identify most variable genes between groups in the two reference datasets. Create erythroid signature gene sets.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.GSE22552 <- readRDS("../processed/combined_RMA_GSE22552.rds")
combined.GSE24759 <- readRDS("../processed/combined_RMA_GSE24759.rds")
combined.RMA.GSE41817.table <- readRDS("../processed/combined.RMA.GSE41817.table.rds")
combined.RMA.GSE89540.table.DBAvControl <- readRDS("../processed/combined.RMA.GSE89540.table.DBAvControl.rds")
combined.RMA.GSE89540.table.DBA_GATA1vControl <- readRDS("../processed/combined.RMA.GSE89540.table.DBA_GATA1vControl.rds")
combined.RMA.GSE89540.sn.table.DBAvControl <- readRDS("../processed/combined.RMA.GSE89540.sn.table.DBAvControl.rds")
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl <- readRDS("../processed/combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.rds")
combined.RMA.GSE89540.sn.norm.table.DBAvControl <- readRDS("../processed/combined.RMA.GSE89540.sn.norm.table.DBAvControl.rds")
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl <- readRDS("../processed/combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.rds")
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage","group1","group2")

# GSE22552 RMA
groups <- samples %>%
  filter(name %in% names(combined.GSE22552)) %>%
  .$group1 %>%
  factor(levels=c("CFU_E","PRO_E","INT_E","LATE_E"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("CFU_E","PRO_E","INT_E","LATE_E")
cont_matrix <- makeContrasts(LATEvEARLY = (INT_E+LATE_E)/2 - (CFU_E+PRO_E)/2,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.GSE22552))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
out <- data.frame(fit$coefficients)
out_table <- topTable(fit, n = 20000, sort.by="logFC", resort.by="logFC", p.value=0.05,adjust.method="fdr",lfc=2)
GSE22552.eryth <- row.names(out_table[out_table$logFC>0,]) #row.names(head(out_table,100))
# GSE24759 RMA
groups <- samples %>%
  filter(name %in% names(combined.GSE24759)) %>%
  .$group1 %>%
  factor(levels=c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos"), ordered=TRUE)
condition <- model.matrix(~0+groups)
colnames(condition) <- c("CD34pos_CD71pos_GlyAneg","CD34neg_CD71pos_GlyAneg","CD34neg_CD71pos_GlyApos","CD34neg_CD71lo_GlyApos","CD34neg_CD71neg_GlyApos")
cont_matrix <- makeContrasts(LATEvEARLY = (CD34neg_CD71pos_GlyApos+CD34neg_CD71lo_GlyApos+CD34neg_CD71neg_GlyApos)/3 - (CD34pos_CD71pos_GlyAneg+CD34neg_CD71pos_GlyAneg)/2,
                             levels = condition)
# Convert to eSet and run limma
v <-new("ExpressionSet", exprs=as.matrix(combined.GSE24759))
fit <- lmFit(v, condition)
fit <- contrasts.fit(fit, cont_matrix)
fit <- eBayes(fit)
out_table <- topTable(fit, n = 20000, sort.by="logFC", resort.by="logFC", p.value=0.05,adjust.method="fdr",lfc=2)
GSE24759.eryth <- row.names(out_table[out_table$logFC>0,]) #row.names(head(out_table,100))
#Read in and merge gene sets
geneset <- read.table("../data/GeneSets.txt",header=T,sep="\t",na.strings="",stringsAsFactors=F)
geneset.l <-lapply(geneset, function (x) x[!is.na(x)])
geneset.eryth <- list(GSE22552.eryth=GSE22552.eryth,GSE24759.eryth=GSE24759.eryth)
genesets <- append(geneset.eryth, geneset.l)

#' RMA-normalized, DBA vs. Control from Ludwig et al. 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE41817.table.ranks <- combined.RMA.GSE41817.table$logFC
names(combined.RMA.GSE41817.table.ranks) <- row.names(combined.RMA.GSE41817.table)
combined.RMA.GSE41817.fgseaRes <- fgsea(genesets, combined.RMA.GSE41817.table.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE41817.table.ranks,combined.RMA.GSE41817.fgseaRes,gseaParam=0.5)

#' RMA-normalized, DBA (RP/I) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.table.DBAvControl.ranks <- combined.RMA.GSE89540.table.DBAvControl$logFC
names(combined.RMA.GSE89540.table.DBAvControl.ranks) <- row.names(combined.RMA.GSE89540.table.DBAvControl)
combined.RMA.GSE41817.DBAvControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.table.DBAvControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.table.DBAvControl.ranks,combined.RMA.GSE41817.DBAvControl.fgseaRes,gseaParam=0.5)

#' RMA-normalized, DBA (GATA1) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.table.DBA_GATA1vControl.ranks <- combined.RMA.GSE89540.table.DBA_GATA1vControl$logFC
names(combined.RMA.GSE89540.table.DBA_GATA1vControl.ranks) <- row.names(combined.RMA.GSE89540.table.DBA_GATA1vControl)
combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.table.DBA_GATA1vControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.table.DBA_GATA1vControl.ranks,combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes,gseaParam=0.5)

#' Synthetic norm of RMA-normalized, DBA (RP/I) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.sn.table.DBAvControl.ranks <- combined.RMA.GSE89540.sn.table.DBAvControl$logFC
names(combined.RMA.GSE89540.sn.table.DBAvControl.ranks) <- row.names(combined.RMA.GSE89540.sn.table.DBAvControl)
combined.RMA.GSE41817.DBAvControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.sn.table.DBAvControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.sn.table.DBAvControl.ranks,combined.RMA.GSE41817.DBAvControl.fgseaRes,gseaParam=0.5)

#' Synethic norm of RMA-normalized, DBA (GATA1) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.ranks <- combined.RMA.GSE89540.sn.table.DBA_GATA1vControl$logFC
names(combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.ranks) <- row.names(combined.RMA.GSE89540.sn.table.DBA_GATA1vControl)
combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.sn.table.DBA_GATA1vControl.ranks,combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes,gseaParam=0.5)

#' Synthetic norm-normalization of RMA-normalized, DBA (RP/I) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.sn.norm.table.DBAvControl.ranks <- combined.RMA.GSE89540.sn.norm.table.DBAvControl$logFC
names(combined.RMA.GSE89540.sn.norm.table.DBAvControl.ranks) <- row.names(combined.RMA.GSE89540.sn.norm.table.DBAvControl)
combined.RMA.GSE41817.DBAvControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.sn.norm.table.DBAvControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.sn.norm.table.DBAvControl.ranks,combined.RMA.GSE41817.DBAvControl.fgseaRes,gseaParam=0.5)

#' Synethic norm-normalization of RMA-normalized, DBA (GATA1) vs. Control from O'Brien et al.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=7, fig.height=4, fig.align='center'
combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.ranks <- combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl$logFC
names(combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.ranks) <- row.names(combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl)
combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes <- fgsea(genesets, combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.ranks, minSize=15, maxSize=500, nperm=1000)
plotGseaTable(genesets,combined.RMA.GSE89540.sn.norm.table.DBA_GATA1vControl.ranks,combined.RMA.GSE41817.DBA_GATA1vControl.fgseaRes,gseaParam=0.5)

