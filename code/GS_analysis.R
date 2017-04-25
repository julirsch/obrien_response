#' ---
#' title: "Analyze GSE89540 with reference data using SCAN (gene sets)"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Compare gene sets used in O'Brien et al. to reference data.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries
library(ggplot2)
library(reshape2)

#' Read in processed data and a subset of gene sets.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.SCAN.bn <- readRDS("../processed/Combined_SCAN_BN.rds")
geneset <- read.table("../data/GeneSets.txt",header=T,sep="\t")
combined.SCAN.bn.ALAS2 <- as.data.frame(combined.SCAN.bn["ALAS2",])

#' Compare gene sets across erythroid maturation or for DBA/healthy for each study.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
#temp <- combined.qn
#temp2 <- cbind(CFUE=rowMeans(temp[,1:3]),PROE=rowMeans(temp[,4:6]),INTE=rowMeans(temp[,7:9]),LATEE=rowMeans(temp[,10:12]))
#temp2 <- cbind(Ery1=rowMeans(temp[,17:23]),Ery2=rowMeans(temp[,24:30]),Ery3=rowMeans(temp[,31:36]),Ery4=rowMeans(temp[,37:43]),Ery5=rowMeans(temp[,44:49]))
#temp2 <- temp[,grep("_44",colnames(temp))]; temp2 <- cbind(DBA_GATA1=rowMeans(temp2[,1:3]),Normal=rowMeans(temp2[,9:16]),DBA=rowMeans(temp2[,c(4:8,17:22)]))
temp2 <- cbind(DBA=rowMeans(temp[,98:100]),Normal=rowMeans(temp[,101:106]))
colMeans(temp2)
temp2.melt <- melt(temp2)
ggplot(temp2.melt, aes(x=value, fill = Var2, color = Var2)) +
  geom_density(alpha=0.1) + 
  scale_fill_manual(values = jdb_palette("Zissou")[c(3,1)]) + 
  scale_color_manual(values = jdb_palette("Zissou")[c(3,1)]) +
  theme_bw()

H1_heme <- combined.SCAN[row.names(combined.SCAN.bn) %in% geneset[,"H1_heme"],]


