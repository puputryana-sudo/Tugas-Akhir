## Analisis data microarray_Differential expression analysis with limma
##LOAD LIBRARY
library(GEOquery)
library(limma)
library(affy)
library(Biobase)
library(dplyr)
library(EnhancedVolcano)
library(VennDiagram)
library(R.utils)

##MICROARRAY DATA
#Input Raw Data
untar("/Users/putrianadwiagustin/Kuliah/TA/GSE214051_RAW.tar", exdir = "CEL")
list.files("CEL")
# Jika masih .gz
# gz_files <- list.files("CEL", pattern = "\\.gz$", full.names = TRUE)
# sapply(gz_files, gunzip, overwrite = TRUE)

#Load & Preprocessing (RMA)
rawData <- ReadAffy(celfile.path = "CEL")
mentah <- exprs(rawData)
colnames(mentah) <- sub("_.*", "", colnames(mentah))

#Filter Sampel
gsm_numbers <- 6598267:6598341
selected_gsm <- paste0("GSM", gsm_numbers)
mentah_filt <- mentah[, colnames(mentah) %in% selected_gsm]
dim(mentah_filt)

#Simpan raw data yang terfilter ke dalam csv
write.csv(mentah_filt, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/Upload_Git/rawdata_GSE214051.csv",
          row.names = TRUE)

#melihat data secara keseluruhan melalui visualisasi
sampleNames(rawData) <- sub("_.*", "", sampleNames(rawData)) # menyederhanakan nama kolom
boxplot(mentah_filt,
        col = "skyblue",
        border = "gray40",
        outline = FALSE,
        las = 2,
        cex.axis = 0.6,
        main = "Distribusi Ekspresi Gen Sebelum Pre-processing",
        ylab = "Ekspresi Gen")

#RMA Normalization
eset <- rma(rawData)
expr_rma <- exprs(eset)
expr_rma_filt <- expr_rma[, colnames(expr_rma) %in% selected_gsm]

#lihat hasil preprocessing
exprs(eset_rma)[1:5, 1:5] #hasilnya sama seperti file getGEO
prepos <- exprs(eset_rma)
RMA_filt <- prepos[, colnames(prepos) %in% selected_gsm] #filter bakteri saja untuk visualisasi
dim(RMA_filt)
boxplot(RMA_filt,
        col = "skyblue",
        border = "gray40",
        outline = FALSE,
        las = 2,
        cex.axis = 0.6,
        main = "Distribusi Ekspresi Gen Setelah Pre-processing",
        ylab = "Ekspresi Gen")

#Load GEO Metadata
GSE214051 <- getGEO("GSE214051", GSEMatrix =TRUE)
gset <- GSE214051[[1]]

#mengubah nama label
fvarLabels(gset) 
fvarLabels(gset) <- make.names(fvarLabels(gset))

#melihat data secara keseluruhan
dim(exprs(gset))
#mengeluarkan gen control
# cek jumlah gen control -> ada 64
sum(grepl("^AFFX", rownames(exprs(gset))))

# hapus gen control
gset <- gset[!grepl("^AFFX", rownames(exprs(gset))), ]##ekspor ke csv
dim(exprs(gset)) #yang awalnya 22690 jadi 22626

#ekspor ke csv
write.csv(exprs(gset), 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/Data/prepos_GSE214051.csv",
          row.names = TRUE)

pData(gset) #melihat info sampel
fData(gset) #melihat info gene/probe
annot<-fData(gset)

pData(gset)[,"agent_time combined:ch1"]
group <- pData(gset)[,"agent_time combined:ch1"]
table(group)

#filter hanya ambil control dan bakteri
keep <- group %in% c("control_24", "control_48", "control_72",
                     "Streptococcus_24", "Streptococcus_48", "Streptococcus_72")
gset_group <- gset[, keep]
agent.time_group <- pData(gset_group)[,"agent_time combined:ch1"]
Infection_Status <- ifelse(grepl("Streptococcus", agent.time_group),"Strep", "Control")
Time <- ifelse(grepl("24", agent.time_group), "24",
               ifelse(grepl("48", agent.time_group), "48", "72"))
Infection_Status <- factor(Infection_Status, levels=c("Control","Strep"))
Time <- factor(Time, levels=c("24","48","72"))
table(Time, Infection_Status)

#Desain Matriks untuk analisis longitudinal
SP_long <- group %in% c("Streptococcus_24", "Streptococcus_48", "Streptococcus_72")
gset_group_long <- gset[, SP_long]
time_group_long <- pData(gset_group_long)[,"agent_time combined:ch1"]
table(time_group_long)
time_group_long <- factor(time_group_long,
                          levels = c("Streptococcus_24", "Streptococcus_48", "Streptococcus_72"))
design_SPlong <- model.matrix(~ 0 + time_group_long)
design_SPlong
colnames(design_SPlong) <- c("Strep24","Strep48","Strep72")
colnames(design_SPlong)

#Contrast S.Pneumonia longitudinal
contrast.matrix_SPlong <- makeContrasts(
  Strep24_vs_Strep48 = Strep48 - Strep24,
  Strep48_vs_Strep72 = Strep72 - Strep48,
  levels = design_SPlong)

#LIMMA
fit_SPlong <- lmFit(gset_group_long, design_SPlong)
fit_SPlong2 <- contrasts.fit(fit_SPlong, contrast.matrix_SPlong)
fit_SPlong2 <- eBayes(fit_SPlong2)

#Hasil DEG
#rules
logFC_cutoff <- 1
padj_cutoff <- 0.05

#Rentang 24-48 jam
res_SPlong24 <- topTable(
  fit_SPlong2,
  coef = "Strep24_vs_Strep48",
  number = Inf,
  adjust.method = "BH")
head(res_SPlong24)
#Filter yang signifikan
deg_SPlong1 <- res_SPlong24[res_SPlong24$adj.P.Val < 0.05 & abs(res_SPlong24$logFC) >= 1,]
nrow(deg_SPlong1)

res_SPlong24$regulation <- "NotSig"

res_SPlong24$regulation[
  res_SPlong24$logFC >= logFC_cutoff &
    res_SPlong24$adj.P.Val < padj_cutoff
] <- "Up"

res_SPlong24$regulation[
  res_SPlong24$logFC <= -logFC_cutoff &
    res_SPlong24$adj.P.Val < padj_cutoff
] <- "Down"

table(res_SPlong24$regulation)


lamp <- data.frame(res_SPlong48$Gene.Symbol,res_SPlong48$logFC,
                   res_SPlong48$adj.P.Val, res_SPlong48$regulation)

deg_lamp <- lamp[lamp$res_SPlong48.regulation %in% c("Up","Down"), ]
#Visualisasi
keyvals_SPlong1 <- ifelse(
  res_SPlong24$regulation == "Up", "red",
  ifelse(res_SPlong24$regulation == "Down", "blue",
         "grey70"))

names(keyvals_SPlong1)[keyvals_SPlong1 == "red"]  <- "Up"
names(keyvals_SPlong1)[keyvals_SPlong1 == "blue"] <- "Down"
names(keyvals_SPlong1)[keyvals_SPlong1 == "grey70"] <- "Not sig"

EnhancedVolcano(res_SPlong24,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                colCustom = keyvals_SPlong1,
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Volcano S.Pneumoniae 48 vs S.Pneumoniae 24",
                subtitle = NULL,
                caption = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

#Rentang 48-72 jam
res_SPlong48 <- topTable(
  fit_SPlong2,
  coef = "Strep48_vs_Strep72",
  number = Inf,
  adjust.method = "BH")
head(res_SPlong48)

#Filter yang signifikan
deg_SPlong2 <- res_SPlong48[res_SPlong48$adj.P.Val < 0.05 & abs(res_SPlong48$logFC) >= 1,]
nrow(deg_SPlong2)

res_SPlong48$regulation <- "NotSig"

res_SPlong48$regulation[
  res_SPlong48$logFC >= logFC_cutoff &
    res_SPlong48$adj.P.Val < padj_cutoff
] <- "Up"

res_SPlong48$regulation[
  res_SPlong48$logFC <= -logFC_cutoff &
    res_SPlong48$adj.P.Val < padj_cutoff
] <- "Down"

table(res_SPlong48$regulation)

#Visualisasi
keyvals_SPlong2 <- ifelse(
  res_SPlong48$regulation == "Up", "red",
  ifelse(res_SPlong48$regulation == "Down", "blue",
         "grey70"))

names(keyvals_SPlong2)[keyvals_SPlong2 == "red"]  <- "Up"
names(keyvals_SPlong2)[keyvals_SPlong2 == "blue"] <- "Down"
names(keyvals_SPlong2)[keyvals_SPlong2 == "grey70"] <- "Not sig"

EnhancedVolcano(res_SPlong48,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                colCustom = keyvals_SPlong2,
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Volcano S.Pneumoniae 72 vs S.Pneumoniae 48",
                subtitle = NULL,
                caption = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

#intersect
int_24 <- rownames(deg_SPlong1)
int_48 <- rownames(deg_SPlong2)
intersect_SP <- intersect(int_24, int_48)
length(intersect_SP)

#Visualisasi diagram venn
venn.plot <- venn.diagram(
  x = list(
    "24–48 jam" = int_24,
    "48–72 jam" = int_48
  ),
  filename = NULL,
  fill = c("lightblue", "steelblue"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),   # posisi label biar ga tabrakan
  cat.dist = c(0.05, 0.05),
  ext.text = FALSE
)
grid::grid.draw(venn.plot)

#Prepare untuk interpretasi biologis
gen <- annot[intersect_SP, "Gene.Symbol"]
length(gen)
table_gen <- data.frame(
  Gene = gen,
  logFC_24 = deg_SPlong1[intersect_SP, "logFC"],
  padj_24 = deg_SPlong1[intersect_SP, "adj.P.Val"],
  logFC_48 = deg_SPlong2[intersect_SP, "logFC"],
  padj_48 = deg_SPlong2[intersect_SP, "adj.P.Val"]
)

write.csv(table_gen, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/GenSPlong.csv",
          row.names = TRUE)

#Hasil analisis DEG global (S.pneumoniae vs Control)
design_SPglobal <- model.matrix(~ 0 + Infection_Status)
colnames(design_SPglobal) <- levels(design_SPglobal)
colnames(design_SPglobal) <- c("Control", "Strep")
design_SPglobal

#Contrast S.Pneumonia global
contrast.matrixSP <- makeContrasts(
  Infection_StatusStrep_vs_Infection_StatusControl = Strep - Control,
  levels = design_SPglobal)

#LIMMA
fit_global <- lmFit(gset_group, design_SPglobal)
fit_SPglobal <- contrasts.fit(fit_global, contrast.matrixSP)
fit_SPglobal <- eBayes(fit_SPglobal)

res_SPglobal <- topTable(
  fit_SPglobal,
  coef = "Infection_StatusStrep_vs_Infection_StatusControl",
  number = Inf,
  adjust.method = "BH")

#Filtering gen yang signifikan
deg_SPglobal <- subset(res_SPglobal, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(deg_SPglobal)

deg_SPglobal$regulation[
  deg_SPglobal$logFC >= logFC_cutoff &
    deg_SPglobal$adj.P.Val < padj_cutoff
] <- "Up"

deg_SPglobal$regulation[
  deg_SPglobal$logFC <= -logFC_cutoff &
    deg_SPglobal$adj.P.Val < padj_cutoff
] <- "Down"

write.csv(deg_SPglobal[,c("Gene.Symbol","logFC","AveExpr","t","P.Value","adj.P.Val","B","regulation")], 
        file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/Gen_SPglobal.csv",
        row.names = TRUE)

res_SPglobal$regulation <- "NotSig"

res_SPglobal$regulation[
  res_SPglobal$logFC >= logFC_cutoff &
    res_SPglobal$adj.P.Val < padj_cutoff
] <- "Up"

res_SPglobal$regulation[
  res_SPglobal$logFC <= -logFC_cutoff &
    res_SPglobal$adj.P.Val < padj_cutoff
] <- "Down"

table(res_SPglobal$regulation)

#Visualisasi
keyvals_SPglobal <- ifelse(
  res_SPglobal$regulation == "Up", "red",
  ifelse(res_SPglobal$regulation == "Down", "blue",
         "grey70"))

names(keyvals_SPglobal)[keyvals_SPglobal == "red"]  <- "Up"
names(keyvals_SPglobal)[keyvals_SPglobal == "blue"] <- "Down"
names(keyvals_SPglobal)[keyvals_SPglobal == "grey70"] <- "Not sig"

EnhancedVolcano(res_SPglobal,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                colCustom = keyvals_SPglobal,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                title = "Volcano S.Pneumoniae vs Control",
                subtitle = NULL,
                caption = NULL,
                
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

#Cleaning untuk PPI
gen_prep <- deg_SPglobal$Gene.Symbol
length(gen_prep)
gen_prep <- gen_prep[gen_prep != "" & !is.na(gen_prep)]
length(gen_prep)
gen_split <- unlist(strsplit(gen_prep, " /// "))
length(gen_split)
gen_unique <- unique(gen_split)
length(gen_unique)

#Gen untuk PPI
write.csv(gen_unique, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/Gen_clust.csv",
          row.names = TRUE)

## NGS DATA
#Input Data
expr_NGS <- read.csv("/Users/putrianadwiagustin/Kuliah/TA/GSE206534_Superinfection_NormCounts.csv", header = TRUE, sep = ";", row.names = 1)
gse <- getGEO("GSE206534", GSEMatrix = TRUE)
pDataSet <- pData(gse[[1]]) 

#Metadata
metadata <- data.frame(
  GEO_Accession = pDataSet$geo_accession,
  Infection_Status = pDataSet$`infection:ch1`,
  Time_Post_Infection = pDataSet$`time (days after first infection):ch1`)
table(metadata$Infection_Status)

#Filter Sampel
colnames(expr_NGS) <- metadata$GEO_Accession
filtr <- metadata[c(1:10,18:35),1]
expr_NGS <- expr_NGS[,filtr]
dim(expr_NGS)

#Metadata 2 untuk analisis longitudinal
metadata_2 <- metadata[c(1:10,18:35),]
metadata_2$Infection_Status <- factor(metadata_2$Infection_Status,
                                      levels = c("Mock infection", "IAV infection"))
metadata_2$Fase <- c("Awal","Awal","Pertengahan","Pertengahan","Pertengahan","Akhir","Akhir","Akhir","Akhir","Akhir",
                     "Awal","Awal","Awal","Awal","Pertengahan","Pertengahan",
                     "Pertengahan","Pertengahan","Pertengahan","Pertengahan","Pertengahan",
                     "Akhir","Akhir","Akhir","Akhir","Akhir",
                     "Akhir","Akhir")
metadata_2$Fase <- factor(metadata_2$Fase,
                          levels = c("Awal","Pertengahan","Akhir"))
table(metadata_2$Infection_Status, metadata_2$Fase)

#visualisasi library size & distribusi ekspresi gen
lib_size <- colSums(expr_NGS)
barplot(
  lib_size,
  names.arg = colnames(lib_size),
  las = 2,                      
  col = "skyblue",
  ylim = c(0, max(lib_size) * 1.1),  
  cex.names = 0.7,          
  main = "Library Size Setiap Sample")

boxplot(expr_NGS,
        col = "skyblue",
        border = "gray40",
        outline = FALSE,
        las = 2,
        cex.axis = 0.6,
        main = "Distribusi Ekspresi Gen",
        ylab = "Ekspresi Gen")

#Recheck data & transformasi log
dim(expr_NGS)
keep.NGS <- rowSums(expr_NGS >= 10) >= 2
expr_filt <- expr_NGS[keep.NGS, ]
dim(expr_filt)
expr_filt <- as.matrix(expr_filt)
expr_log <- log2(expr_filt + 1)
dim(expr_log)

#Analisis longitudinal IAV
#Desain matriks
metadata_iav <- metadata_2[metadata_2$Infection_Status == "IAV infection", ]
metadata_iav$Fase <- factor(metadata_iav$Fase, 
                            levels = c("Awal", "Pertengahan", "Akhir"))
design_IAVlong <- model.matrix(~ 0 + metadata_iav$Fase )
colnames(design_IAVlong) <- c("Awal","Tengah","Akhir")

#Contrast
contrast.matrix_IAVlong <- makeContrasts(
  Awal_vs_Tengah = Tengah - Awal,
  Tengah_vs_Akhir = Akhir - Tengah,
  levels = design_IAVlong)

#LIMMA
GSM_IAVlong <- metadata_iav[,1]
long_IAV <- expr_log[,GSM_IAVlong]
fit_IAVlong <- lmFit(long_IAV, design_IAVlong)
fit_IAVlong2 <- contrasts.fit(fit_IAVlong, contrast.matrix_IAVlong)
fit_IAVlong2 <- eBayes(fit_IAVlong2)

#Rentang Fase awal-tengah
res_IAVlongawal <- topTable(
  fit_IAVlong2,
  coef = "Awal_vs_Tengah",
  number = Inf,
  adjust.method = "BH")
head(res_IAVlongawal)

#Filter gen signifikan
deg_IAVlong1 <- res_IAVlongawal[res_IAVlongawal$adj.P.Val < 0.05 & abs(res_IAVlongawal$logFC) >= 1,]
nrow(deg_IAVlong1)

res_IAVlongawal$regulation <- "NotSig"

res_IAVlongawal$regulation[
  res_IAVlongawal$logFC >= logFC_cutoff &
    res_IAVlongawal$adj.P.Val < padj_cutoff
] <- "Up"

res_IAVlongawal$regulation[
  res_IAVlongawal$logFC <= -logFC_cutoff &
    res_IAVlongawal$adj.P.Val < padj_cutoff
] <- "Down"

table(res_IAVlongawal$regulation)

write.csv(res_IAVlongawal, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/IAlong1.csv",
          row.names = TRUE)

#Visualisasi
keyvals_IAVlong1 <- ifelse(
  res_IAVlongawal$regulation == "Up", "red",
  ifelse(res_IAVlongawal$regulation == "Down", "blue",
         "grey70"))

names(keyvals_IAVlong1)[keyvals_IAVlong1 == "red"]  <- "Up"
names(keyvals_IAVlong1)[keyvals_IAVlong1 == "blue"] <- "Down"
names(keyvals_IAVlong1)[keyvals_IAVlong1 == "grey70"] <- "Not sig"

EnhancedVolcano(res_IAVlongawal,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                colCustom = keyvals_IAVlong1,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                title = "Volcano IAV Tengah vs Awal",
                subtitle = NULL,
                caption = NULL,
                
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

#Rentang Fase tengah-akhir
res_IAVlongtengah <- topTable(
  fit_IAVlong2,
  coef = "Tengah_vs_Akhir",
  number = Inf,
  adjust.method = "BH")
head(res_IAVlongtengah)

#Filter gen signifikan
deg_IAVlong2 <- res_IAVlongtengah[res_IAVlongtengah$adj.P.Val < 0.05 & abs(res_IAVlongtengah$logFC) >= 1,]
nrow(deg_IAVlong2)

res_IAVlongtengah$regulation <- "NotSig"

res_IAVlongtengah$regulation[
  res_IAVlongtengah$logFC >= logFC_cutoff &
    res_IAVlongtengah$adj.P.Val < padj_cutoff
] <- "Up"

res_IAVlongtengah$regulation[
  res_IAVlongtengah$logFC <= -logFC_cutoff &
    res_IAVlongtengah$adj.P.Val < padj_cutoff
] <- "Down"

table(res_IAVlongtengah$regulation)

#Visualisasi
keyvals_IAVlong2 <- ifelse(
  res_IAVlongtengah$regulation == "Up", "red",
  ifelse(res_IAVlongtengah$regulation == "Down", "blue",
         "grey70"))

names(keyvals_IAVlong2)[keyvals_IAVlong2 == "red"]  <- "Up"
names(keyvals_IAVlong2)[keyvals_IAVlong2 == "blue"] <- "Down"
names(keyvals_IAVlong2)[keyvals_IAVlong2 == "grey70"] <- "Not sig"

EnhancedVolcano(res_IAVlongtengah,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                colCustom = keyvals_IAVlong2,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                title = "Volcano IAV Akhir vs Tengah",
                subtitle = NULL,
                caption = NULL,
                
                gridlines.major = FALSE,
                gridlines.minor = FALSE)


#Analisis global (IAV vs Control)
Infection_StatusIAV <- metadata_2[,"Infection_Status"]
Infection_StatusIAV <- gsub(" infection", "", Infection_StatusIAV)
Infection_StatusIAV <- factor(Infection_StatusIAV, levels=c("Mock","IAV"))

#Desain Matriks
design_IAVglobal <- model.matrix(~ 0 + Infection_StatusIAV )
colnames(design_IAVglobal)<- c("Mock", "IAV")
design_IAVglobal

#Contrast
contrast.matrixIAV <- makeContrasts(
  IAV_infection_vs_Mock_infection = IAV - Mock,
  levels = design_IAVglobal)

#LIMMA
fit_global2 <- lmFit(expr_log, design_IAVglobal)
fit_IAVglobal <- contrasts.fit(fit_global2, contrast.matrixIAV)
fit_IAVglobal <- eBayes(fit_IAVglobal)
res_IAVglobal <- topTable(
  fit_IAVglobal,
  number = Inf,
  adjust.method = "BH")

#Filtering gen signifikan
sig_IAVglobal <- res_IAVglobal[res_IAVglobal$adj.P.Val < 0.05 & abs(res_IAVglobal$logFC) > 1,]
nrow(sig_IAVglobal)

res_IAVglobal$regulation <- "NotSig"

res_IAVglobal$regulation[
  res_IAVglobal$logFC >= logFC_cutoff &
    res_IAVglobal$adj.P.Val < padj_cutoff
] <- "Up"

res_IAVglobal$regulation[
  res_IAVglobal$logFC <= -logFC_cutoff &
    res_IAVglobal$adj.P.Val < padj_cutoff
] <- "Down"
table(res_IAVglobal$regulation)

write.csv(res_IAVglobal, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/IAglobal1.csv",
          row.names = TRUE)

#Visualisasi
keyvals_IAVglobal <- ifelse(
  res_IAVglobal$regulation == "Up", "red",
  ifelse(res_IAVglobal$regulation == "Down", "blue",
         "grey70"))

names(keyvals_IAVglobal)[keyvals_IAVglobal == "red"] <- "Up"
names(keyvals_IAVglobal)[keyvals_IAVglobal == "blue"] <- "Down"
names(keyvals_IAVglobal)[keyvals_IAVglobal == "grey70"] <- "NotSig"

EnhancedVolcano(
  res_IAVglobal,
  
  lab = NA,  
  selectLab = NULL,
  
  x = "logFC",
  y = "adj.P.Val",
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  colCustom = keyvals_IAVglobal,
  
  title = "Volcano Plot IAV vs Control",
  subtitle = NULL,
  caption = NULL,
  
  gridlines.major = FALSE,
  gridlines.minor = FALSE)

#prepare untuk PPI
gen_IAV <- rownames(sig_IAVglobal)
length(gen_IAV)
gen_IAV <- gen_IAV[gen_IAV != "" & !is.na(gen_IAV)]
length(gen_IAV)
gen_IAV <- unique(gen_IAV)
length(gen_IAV)

write.csv(gen_IAV, 
          file = "/Users/putrianadwiagustin/Library/CloudStorage/OneDrive-UniversitasIslamIndonesia/Kuliah/Bismillah Tugas Akhir/output/PPI_IAV.csv",
          row.names = TRUE)
