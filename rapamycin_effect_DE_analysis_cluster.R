library(DESeq2)
library(apeglm)
library(vsn)
library(IHW)
library(ashr)

FDR.Thresh <- 0.05
LFC.Thresh <- 0.585
set.seed(1234)

cts <- as.matrix(read.csv("rapa_count_table_sim_mel.csv", row.names = 1))
cts <- subset(cts, select = -OOFCH4)

coldata <- read.csv("rapa_group_data_sim_mel.csv", row.names = 1)
coldata <- coldata[colnames(cts),]

coldata$SpeciesTreatment <- paste0(coldata$Species, coldata$Treatment)
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Tissue <- as.factor(coldata$Tissue)
coldata$Species <- as.factor(coldata$Species)
coldata$Group <- as.factor(coldata$Group)
coldata$Sex <- as.factor(coldata$Sex)
coldata$Treatment <- relevel(coldata$Treatment, "Rapa")
coldata$SpeciesTreatment <-as.factor(coldata$SpeciesTreatment)


ifelse(!dir.exists("species_treatment_effect"), dir.create("species_treatment_effect", recursive=TRUE), print("directory already exists"))


## Female Head

cts_fh <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="H"] 
coldata_fh <- coldata[colnames(cts_fh),]
dds_fh <- DESeqDataSetFromMatrix(countData = cts_fh,
                                 colData = coldata_fh,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_fh)) >= 1 #minimal pre-filtering
dds_fh <- dds_fh[keep,]
dds_fh <- DESeq(dds_fh)

res_fh_OOvsSS_C <- results(dds_fh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_fh_OOvsSS_C, "species_treatment_effect/res_fh_OOvsSS_C.rds")

res_fh_OOvsSS_R <- results(dds_fh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_fh_OOvsSS_R, "species_treatment_effect/res_fh_OOvsSS_R.rds")

res_fh_RvsC_OO <- results(dds_fh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)

saveRDS(res_fh_RvsC_OO, "species_treatment_effect/res_fh_RvsC_OO.rds")

res_fh_RvsC_SS <- results(dds_fh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_fh_RvsC_SS, "species_treatment_effect/res_fh_RvsC_SS.rds")


res_fh_RvsC_OO <- lfcShrink(dds_fh,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_fh_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_fh_RvsC_OO, "species_treatment_effect/res_fh_RvsC_OO.rds")

res_fh_RvsC_SS <- lfcShrink(dds_fh,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_fh_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_fh_RvsC_SS, "species_treatment_effect/res_fh_RvsC_SS.rds")

## Female Thorax

cts_ft <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="T"] 
coldata_ft <- coldata[colnames(cts_ft),]
dds_ft <- DESeqDataSetFromMatrix(countData = cts_ft,
                                 colData = coldata_ft,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_ft)) >= 1 #minimal pre-filtering
dds_ft <- dds_ft[keep,]
dds_ft <- DESeq(dds_ft)

res_ft_OOvsSS_C <- results(dds_ft, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_ft_OOvsSS_C, "species_treatment_effect/res_ft_OOvsSS_C.rds")

res_ft_OOvsSS_R <- results(dds_ft, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_ft_OOvsSS_R, "species_treatment_effect/res_ft_OOvsSS_R.rds")

res_ft_RvsC_OO <- results(dds_ft, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)
saveRDS(res_ft_RvsC_OO, "species_treatment_effect/res_ft_RvsC_OO.rds")

res_ft_RvsC_SS <- results(dds_ft, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_ft_RvsC_SS, "species_treatment_effect/res_ft_RvsC_SS.rds")

res_ft_RvsC_OO <- lfcShrink(dds_ft,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_ft_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_ft_RvsC_OO, "species_treatment_effect/res_ft_RvsC_OO.rds")

res_ft_RvsC_SS <- lfcShrink(dds_ft,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_ft_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_ft_RvsC_SS, "species_treatment_effect/res_ft_RvsC_SS.rds")


## Female Abdomen

cts_fa <- cts[,substring(colnames(cts), 3, 3)=="F" & substring(colnames(cts), 5, 5)=="A"] 
coldata_fa <- coldata[colnames(cts_fa),]
dds_fa <- DESeqDataSetFromMatrix(countData = cts_fa,
                                 colData = coldata_fa,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_fa)) >= 1 #minimal pre-filtering
dds_fa <- dds_fa[keep,]
dds_fa <- DESeq(dds_fa)

res_fa_OOvsSS_C <- results(dds_fa, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_fa_OOvsSS_C, "species_treatment_effect/res_fa_OOvsSS_C.rds")

res_fa_OOvsSS_R <- results(dds_fa, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_fa_OOvsSS_R, "species_treatment_effect/res_fa_OOvsSS_R.rds")

res_fa_RvsC_OO <- results(dds_fa, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)
saveRDS(res_fa_RvsC_OO, "species_treatment_effect/res_fa_RvsC_OO.rds")

res_fa_RvsC_SS <- results(dds_fa, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_fa_RvsC_SS, "species_treatment_effect/res_fa_RvsC_SS.rds")

res_fa_RvsC_OO <- lfcShrink(dds_fa,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_fa_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_fa_RvsC_OO, "species_treatment_effect/res_fa_RvsC_OO.rds")

res_fa_RvsC_SS <- lfcShrink(dds_fa,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_fa_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_fa_RvsC_SS, "species_treatment_effect/res_fa_RvsC_SS.rds")

## Male Head

cts_mh <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="H"] 
coldata_mh <- coldata[colnames(cts_mh),]
dds_mh <- DESeqDataSetFromMatrix(countData = cts_mh,
                                 colData = coldata_mh,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_mh)) >= 1 #minimal pre-filtering
dds_mh <- dds_mh[keep,]
dds_mh <- DESeq(dds_mh)

res_mh_OOvsSS_C <- results(dds_mh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_mh_OOvsSS_C, "species_treatment_effect/res_mh_OOvsSS_C.rds")

res_mh_OOvsSS_R <- results(dds_mh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_mh_OOvsSS_R, "species_treatment_effect/res_mh_OOvsSS_R.rds")

res_mh_RvsC_OO <- results(dds_mh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)
saveRDS(res_mh_RvsC_OO, "species_treatment_effect/res_mh_RvsC_OO.rds")

res_mh_RvsC_SS <- results(dds_mh, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_mh_RvsC_SS, "species_treatment_effect/res_mh_RvsC_SS.rds")

res_mh_RvsC_OO <- lfcShrink(dds_mh,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_mh_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_mh_RvsC_OO, "species_treatment_effect/res_mh_RvsC_OO.rds")

res_mh_RvsC_SS <- lfcShrink(dds_mh,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_mh_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_mh_RvsC_SS, "species_treatment_effect/res_mh_RvsC_SS.rds")

## Male Thorax

cts_mt <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="T"] 
coldata_mt <- coldata[colnames(cts_mt),]
dds_mt <- DESeqDataSetFromMatrix(countData = cts_mt,
                                 colData = coldata_mt,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_mt)) >= 1 #minimal pre-filtering
dds_mt <- dds_mt[keep,]
dds_mt <- DESeq(dds_mt)

res_mt_OOvsSS_C <- results(dds_mt, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_mt_OOvsSS_C, "species_treatment_effect/res_mt_OOvsSS_C.rds")

res_mt_OOvsSS_R <- results(dds_mt, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_mt_OOvsSS_R, "species_treatment_effect/res_mt_OOvsSS_R.rds")

res_mt_RvsC_OO <- results(dds_mt, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)
saveRDS(res_mt_RvsC_OO, "species_treatment_effect/res_mt_RvsC_OO.rds")

res_mt_RvsC_SS <- results(dds_mt, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_mt_RvsC_SS, "species_treatment_effect/res_mt_RvsC_SS.rds")

res_mt_RvsC_OO <- lfcShrink(dds_mt,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_mt_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_mt_RvsC_OO, "species_treatment_effect/res_mt_RvsC_OO.rds")

res_mt_RvsC_SS <- lfcShrink(dds_mt,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_mt_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_mt_RvsC_SS, "species_treatment_effect/res_mt_RvsC_SS.rds")

## Male Abdomen

cts_ma <- cts[,substring(colnames(cts), 3, 3)=="M" & substring(colnames(cts), 5, 5)=="A"] 
coldata_ma <- coldata[colnames(cts_ma),]
dds_ma <- DESeqDataSetFromMatrix(countData = cts_ma,
                                 colData = coldata_ma,
                                 design = ~ SpeciesTreatment) #design matrix with no intercept
keep <- rowSums(counts(dds_ma)) >= 1 #minimal pre-filtering
dds_ma <- dds_ma[keep,]
dds_ma <- DESeq(dds_ma)

res_ma_OOvsSS_C <- results(dds_ma, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OOControl","SSControl"), filterFun = ihw)
saveRDS(res_ma_OOvsSS_C, "species_treatment_effect/res_ma_OOvsSS_C.rds")

res_ma_OOvsSS_R <- results(dds_ma, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                           contrast=c("SpeciesTreatment","OORapa","SSRapa"), filterFun = ihw)
saveRDS(res_ma_OOvsSS_R, "species_treatment_effect/res_ma_OOvsSS_R.rds")

res_ma_RvsC_OO <- results(dds_ma, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","OORapa","OOControl"), filterFun = ihw)
saveRDS(res_ma_RvsC_OO, "species_treatment_effect/res_ma_RvsC_OO.rds")
res_ma_RvsC_SS <- results(dds_ma, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                          contrast=c("SpeciesTreatment","SSRapa","SSControl"), filterFun = ihw)
saveRDS(res_ma_RvsC_SS, "species_treatment_effect/res_ma_RvsC_SS.rds")

res_ma_RvsC_OO <- lfcShrink(dds_ma,
                                   contrast = c("SpeciesTreatment", "OOControl", "OORapa"),
                                   res = res_ma_RvsC_OO,     
                                   type = "ashr")
saveRDS(res_ma_RvsC_OO, "species_treatment_effect/res_ma_RvsC_OO.rds")

res_ma_RvsC_SS <- lfcShrink(dds_ma,
                                   contrast = c("SpeciesTreatment", "SSControl", "SSRapa"),
                                   res = res_ma_RvsC_SS,     
                                   type = "ashr")
saveRDS(res_ma_RvsC_SS, "species_treatment_effect/res_ma_RvsC_SS.rds")