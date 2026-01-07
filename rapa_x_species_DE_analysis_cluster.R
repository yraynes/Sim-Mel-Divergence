params <-list(FDR = 0.05, LFC = 0)

## ----include=FALSE------------------------------------------------------------
FDR.Thresh <- params$FDR
LFC.Thresh <- params$LFC
set.seed(1234)


## ----message=FALSE, warning = FALSE-------------------------------------------
library(IHW)
library(DESeq2)
library(apeglm)
library(vsn)

## ----data acquisition---------------------------------------------------------

cts <- as.matrix(read.csv("rapa_count_table_sim_mel.csv", row.names = 1))
cts <- subset(cts, select = -OOFCH4)

gene_universe_all <- row.names(cts[rowSums(cts)>0, ])
length(gene_universe_all)

gene_universe_a <- row.names(cts[, ])

cts_a <- cts[,substring(colnames(cts), 5, 5)=="A"] 
gene_universe_a <- row.names(cts_a[rowSums(cts_a)>0, ])

cts_t <- cts[,substring(colnames(cts), 5, 5)=="T"] 
gene_universe_t <- row.names(cts_t[rowSums(cts_t)>0, ])

cts_h <- cts[,substring(colnames(cts), 5, 5)=="H"] 
gene_universe_h <- row.names(cts_h[rowSums(cts_h)>0, ])

ifelse(!dir.exists("results/All_Data"), dir.create("results/All_Data", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_all, "results/All_Data/gene_universe_ad.txt")

coldata <- read.csv("rapa_group_data_sim_mel.csv", row.names = 1)
coldata <- coldata[colnames(cts),]

coldata$SpeciesTreatment <- paste0(coldata$Species, coldata$Treatment)

#save experimental features as factors
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Tissue <- as.factor(coldata$Tissue)
coldata$Species <- as.factor(coldata$Species)
coldata$Group <- as.factor(coldata$Group)
coldata$Sex <- as.factor(coldata$Sex)
coldata$SpeciesTreatment <-as.factor(coldata$SpeciesTreatment)

## ----subset count data to female heads----------------------------------------

cts_fh <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="H"] 
coldata_fh <- coldata[colnames(cts_fh),]

gene_universe_fh <- rownames(cts[rowSums(cts_fh)>0, ])
length(gene_universe_fh)

ifelse(!dir.exists("results/Female_Head"), dir.create("results/Female_Head", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_fh, "results/Female_Head/gene_universe_fh.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_fh) %in% colnames(cts_fh)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_fh) == colnames(cts_fh)) #are they in the same order?


## ----make a new DESeqDataSet for female heads---------------------------------
dds_fh <- DESeqDataSetFromMatrix(countData = cts_fh,
                              colData = coldata_fh,
                              design = ~ Species + Treatment + Species:Treatment) 
keep <- rowSums(counts(dds_fh)) >= 1 #minimal pre-filtering
dds_fh <- dds_fh[keep,]


## ----make a new DESeqDataSet for female heads RT interaction------------------

dds_fh <- DESeq(dds_fh)
resultsNames(dds_fh)
fhi <- results(dds_fh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
saveRDS(fhi, "fhi.rds")
summary(fhi)
fhi_Ordered <- fhi[order(fhi$pvalue),]
fhi_Sig <- subset(fhi_Ordered, padj < FDR.Thresh) #507 gene
dim(fhi_Sig)
#write.csv(fhi_Sig, "results/Female_Head/fh_sp_rapa.csv")
gc()

## ----subset count data to female thorax---------------------------------------

cts_ft <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="T"] 
coldata_ft <- coldata[colnames(cts_ft),]

gene_universe_ft <- rownames(cts[rowSums(cts_ft)>0, ])
length(gene_universe_ft)

ifelse(!dir.exists("results/Female_Thorax"), dir.create("results/Female_Thorax", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_ft, "results/Female_Thorax/gene_universe_ft.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_ft) %in% colnames(cts_ft)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_ft) == colnames(cts_ft)) #are they in the same order?


## ----make a new DESeqDataSet for female thorax--------------------------------
dds_ft <- DESeqDataSetFromMatrix(countData = cts_ft,
                              colData = coldata_ft,
                              design = ~ Species + Treatment + Species:Treatment) #design matrix with no intercept

keep <- rowSums(counts(dds_ft)) >= 1 #minimal pre-filtering
dds_ft <- dds_ft[keep,]


## ----make a new DESeqDataSet for female thorax RT interaction-----------------
dds_ft <- DESeq(dds_ft)
resultsNames(dds_ft)
fti <- results(dds_ft, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
summary(fti)
saveRDS(fti, "fti.rds")

fti_Ordered <- fti[order(fti$pvalue),]
fti_Sig <- subset(fti_Ordered, padj < FDR.Thresh) #467 genes
dim(fti_Sig)
#write.csv(fti_Sig, "results/Female_Thorax/ft_sp_rapa.csv")
gc()

## ----subset count data to female abdomen--------------------------------------

cts_fa <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="A"] 
coldata_fa <- coldata[colnames(cts_fa),]

gene_universe_fa <- rownames(cts[rowSums(cts_fa)>0, ])
length(gene_universe_fa)

ifelse(!dir.exists("results/Female_Abdomen"), dir.create("results/Female_Abdomen", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_fa, "results/Female_Abdomen/gene_universe_fa.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_fa) %in% colnames(cts_fa)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_fa) == colnames(cts_fa)) #are they in the same order?


## ----make a new DESeqDataSet for female abdomen-------------------------------
dds_fa <- DESeqDataSetFromMatrix(countData = cts_fa,
                              colData = coldata_fa,
                              design = ~ Species + Treatment + Species:Treatment) #design matrix with no intercept

keep <- rowSums(counts(dds_fa)) >= 1 #minimal pre-filtering
dds_fa <- dds_fa[keep,]


## ----make a new DESeqDataSet for female abdomen RT interaction----------------
dds_fa <- DESeq(dds_fa)
resultsNames(dds_fa)
fai <- results(dds_fa, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
summary(fai)
saveRDS(fai, "fai.rds")

fai_Ordered <- fai[order(fai$pvalue),]
fai_Sig <- subset(fai_Ordered, padj < FDR.Thresh) #430 genes
dim(fai_Sig)
#write.csv(fai_Sig, "results/Female_Abdomen/fa_sp_rapa.csv")
gc()


## ----subset count data to male heads------------------------------------------

cts_mh <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="H"] 
coldata_mh <- coldata[colnames(cts_mh),]

gene_universe_mh <- rownames(cts[rowSums(cts_mh)>0, ])
length(gene_universe_mh)

ifelse(!dir.exists("results/Male_Head"), dir.create("results/Male_Head", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_mh, "results/Male_Head/gene_universe_mh.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_mh) %in% colnames(cts_mh)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_mh) == colnames(cts_mh)) #are they in the same order?


## ----make a new DESeqDataSet for male heads-----------------------------------
dds_mh <- DESeqDataSetFromMatrix(countData = cts_mh,
                              colData = coldata_mh,
                              design = ~ Species + Treatment + Species:Treatment) 
keep <- rowSums(counts(dds_mh)) >= 1 #minimal pre-filtering
dds_mh <- dds_mh[keep,]


## ----make a new DESeqDataSet for male heads RT interaction--------------------
dds_mh <- DESeq(dds_mh)
resultsNames(dds_mh)
mhi <- results(dds_mh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
summary(mhi)
saveRDS(mhi, "mhi.rds")

mhi_Ordered <- mhi[order(mhi$pvalue),]
mhi_Sig <- subset(mhi_Ordered, padj < FDR.Thresh) #555 gene
dim(mhi_Sig)
#write.csv(mhi_Sig, "results/Male_Head/mh_sp_rapa.csv")
gc()


## ----subset count data to male thorax-----------------------------------------

cts_mt <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="T"] 
coldata_mt <- coldata[colnames(cts_mt),]

gene_universe_mt <- rownames(cts[rowSums(cts_mt)>0, ])
length(gene_universe_mt)

ifelse(!dir.exists("results/Male_Thorax"), dir.create("results/Male_Thorax", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_mt, "results/Male_Thorax/gene_universe_mt.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_mt) %in% colnames(cts_mt)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_mt) == colnames(cts_mt)) #are they in the same order?


## ----make a new DESeqDataSet for male thorax----------------------------------
dds_mt <- DESeqDataSetFromMatrix(countData = cts_mt,
                              colData = coldata_mt,
                              design = ~ Species + Treatment + Species:Treatment) #design matrix with no intercept

keep <- rowSums(counts(dds_mt)) >= 1 #minimal pre-filtering
dds_mt <- dds_mt[keep,]


## ----make a new DESeqDataSet for male thorax RT interaction-------------------
dds_mt <- DESeq(dds_mt)
resultsNames(dds_mt)
mti <- results(dds_mt, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
summary(mti)
saveRDS(mti, "mti.rds")


mti_Ordered <- mti[order(mti$pvalue),]
mti_Sig <- subset(mti_Ordered, padj < FDR.Thresh) #59 genes

mtiShrunk <- lfcShrink(dds_mt, coef="SpeciesSS.TreatmentRapa", type="apeglm")

dim(mti_Sig)
#write.csv(mti_Sig, "results/Male_Thorax/mt_sp_rapa.csv")
gc()

## ----subset count data to male abdomen----------------------------------------

cts_ma <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="A"] 
coldata_ma <- coldata[colnames(cts_ma),]

gene_universe_ma <- rownames(cts[rowSums(cts_ma)>0, ])
length(gene_universe_ma)

ifelse(!dir.exists("results/Male_Abdomen"), dir.create("results/Male_Abdomen", recursive=TRUE), print("directory already exists"))
#writeLines(gene_universe_ma, "results/Male_Abdomen/gene_universe_ma.txt")

#columns of the count matrix and the rows of the column data (information about samples) must be the same so check:
all(rownames(coldata_ma) %in% colnames(cts_ma)) #all rownames of coldate match colnames of count matrix
all(rownames(coldata_ma) == colnames(cts_ma)) #are they in the same order?


## ----make a new DESeqDataSet for male abdomen---------------------------------
dds_ma <- DESeqDataSetFromMatrix(countData = cts_ma,
                              colData = coldata_ma,
                              design = ~ Species + Treatment + Species:Treatment) #design matrix with no intercept

keep <- rowSums(counts(dds_ma)) >= 1 #minimal pre-filtering
dds_ma <- dds_ma[keep,]


## ----make a new DESeqDataSet for Male abdomen RT interaction------------------
dds_ma <- DESeq(dds_ma)
resultsNames(dds_ma)
mai <- results(dds_ma, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh, name="SpeciesSS.TreatmentRapa", filterFun=ihw)
summary(mai)
saveRDS(mai, "mai.rds")

mai_Ordered <- mai[order(mai$pvalue),]
mai_Sig <- subset(mai_Ordered, padj < FDR.Thresh) #557 genes
dim(mai_Sig)
#write.csv(mai_Sig, "results/Male_Abdomen/ma_sp_rapa.csv")
gc()

