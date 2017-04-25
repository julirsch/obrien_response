#' ---
#' title: "Analyze GSE89540 with reference data using CIBERSORT"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Here we create both reference and test datasets for CIBERSORT, which has an online submission portal (https://cibersort.stanford.edu). We then read back into R the results and analyze them.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
# Import libraries

#' Read in both SCAN- and RMA-normalized combined datasets with and without BN
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
combined.SCAN <- readRDS("../processed/Combined_SCAN.rds")
combined.SCAN.bn <- readRDS("../processed/Combined_SCAN_BN.rds")

#' Load in sample annotation
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
samples <- read.table("../data/Samples.txt",sep="\t",stringsAsFactors=F)
names(samples) <- c("GSM","name","GSE","stage")
samples$name <- factor(samples$name,levels=names(combined.SCAN),ordered=T)
samples <- samples[order(samples$name),]

