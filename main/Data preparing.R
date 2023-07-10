#---------------------------------------------------#
# prepare input data for MOVICS (tumor sample only) #

# set working path 
workdir <- "~/Ruan/KIRC/data"
setwd(workdir)

# load R package
library(TCGAbiolinks)
library(data.table)
library(ChAMPdata)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)

# customized functions
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# load annotation
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#----------------------#
# Clinical information #
clinical <- GDCquery(project = "TCGA-KIRC", 
                     data.category = "Clinical", 
                     file.type = "xml")
GDCdownload(clinical,directory = "GDCdata")
cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient",directory = "GDCdata")
colnames(cliquery)[1] <- "Tumor_Sample_Barcode"
write.table(cliquery,"TCGA_KIRC_Clinical.txt",  quote=F, row.names=F,sep = "\t")

#----------------#
# RNA expression #

# count data
expquery <- GDCquery(project = "TCGA-KIRC", 
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
expMatrix <- TCGAanalyze_Preprocessing(expquery2)
colnames(expMatrix) <- substr(colnames(expMatrix), start = 1,stop = 15)
normsamples <- colnames(expMatrix)[which(substr(colnames(expMatrix),14,15) == "11")] # get normal samples
tumorsamples <- colnames(expMatrix)[which(substr(colnames(expMatrix),14,15) == "01")] # get tumor samples
expMatrix <- expMatrix[,tumorsamples]; expMatrix <- as.data.frame(expMatrix[rowSums(expMatrix) > 0,])
#expMatrix <- expMatrix[apply(expMatrix, 1, function(x) sum(x > 0) > 0.9*length(tumorsamples)),]
comgene <- intersect(rownames(expMatrix),rownames(Ginfo))
count <- as.data.frame(expMatrix)[comgene,]; Ginfo <- Ginfo[comgene,]
count$gene <- Ginfo$genename; count <- count[!duplicated(count$gene),]; Ginfo <- Ginfo[rownames(count),]; rownames(count) <- count$gene; count <- count[,-ncol(count)]
colnames(count) <- substr(colnames(count), start = 1,stop = 15)
write.table(count, "TCGA_KIRC_Count.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

# extract mRNA and lncRNA
Lids <- Ginfo[Ginfo$genetype %in% c("lincRNA","antisense","processed_transcript","sense_intronic","sense_overlapping","3prime_overlapping_ncrna","non_coding"),"genename"]
Mids <- Ginfo[Ginfo$genetype == "protein_coding","genename"]

# FPKM
fpkm <- as.data.frame(round(apply(as.matrix(count), 2, countToFpkm, effLen = Ginfo$unqlen),2))
rownames(fpkm) <- rownames(count)
write.table(fpkm, "TCGA_KIRC_FPKM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

# TPM
tpm <- as.data.frame(round(apply(fpkm,2,fpkmToTpm),2))
write.table(tpm, "TCGA_KIRC_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

#-----------------#
# DNA methylation #

# please download this in XENA (GDC TCGA Lung Cancer (KIRC))
# https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
meth <- fread("TCGA-KIRC.methylation450.tsv.gz",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
meth <- as.data.frame(meth); rownames(meth) <- meth[,1]; meth <- as.data.frame(na.omit(meth[,-1]))
colnames(meth) <- substr(colnames(meth), start = 1,stop = 15)
meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]

# get promoter
data("probe.features")
promoter <- intersect(rownames(probe.features[which(probe.features$cgi == "island" & probe.features$feature %in% c("TSS1500","TSS200")),]),rownames(meth))
meth <- round(meth[promoter,],3)

# map to gene name to take median value if multiple mapping
meth$gene <- as.character(probe.features[promoter,"gene"])
meth <- apply(meth[,setdiff(colnames(meth), "gene")], 2, function(x) tapply(x, INDEX=factor(meth$gene), FUN=median, na.rm=TRUE)) # be patient because this will take a while
meth <- as.data.frame(round(meth,3))

#---------------------#
# Copy number segment #

# please download this in XENA (GDC TCGA Lung Cancer (KIRC))
# https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.masked_cnv.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
segment <- fread("TCGA-KIRC.masked_cnv.tsv.gz",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
colnames(segment) <- c("sample","chrom","start","end","value")
segment$sample <- substr(segment$sample, start = 1,stop = 15)
segment <- segment[which(substr(segment$sample, 14, 15) == "01"),]
write.table(segment, "TCGA_KIRC_segment.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

#-----#
# MAF #

# please download this from cBioPortal under Data Sets archive
# https://www.cbioportal.org/datasets
maf <- read_tsv("data_mutations.txt", comment = "#")
flag <- read.table("Mutation Flags 100.txt",sep = "\t")$V1
label <- c("Tumor_Sample_Barcode",
           "Hugo_Symbol",
           "Chromosome",
           "Start_Position",
           "End_Position",
           "Variant_Classification",
           "Variant_Type",
           "Reference_Allele",
           "Tumor_Seq_Allele1",
           "Tumor_Seq_Allele2")
maf <- maf[-which(maf$Hugo_Symbol %in% flag),label]
maf$Hugo_Symbol <- toupper(maf$Hugo_Symbol) # transfer gene name to capital. (eg, C1orf198 to C1ORF198)
write.table(maf, "TCGA_KIRC_MAF.txt", quote=F, row.names=F,col.names = T,sep = "\t")

#-------------------------#
# binary somatic mutation #

# option 1: please download from XENA (TCGA Lung Cancer (KIRC))
# https://xenabrowser.net/datapages/?dataset=mc3_gene_level%2FLIHC_mc3_gene_level.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# mut <- read.delim("LIHC_mc3_gene_level.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
# mut <- as.data.frame(na.omit(mut))
# rownames(mut) <- mut$sample; mut <- mut[,-1]
# mut <- mut[rowSums(mut) > 0,]
# mut <- mut[,which(substr(colnames(mut), 14, 15) == "01")]
# write.table(mut, "TCGA_LIHC_mut1.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

# option 2: convert from MAF file without silent mutation (may not be the same with above file)
# binary mutation
mut.binary <- matrix(0,nrow = length(unique(maf$Hugo_Symbol)),ncol = length(unique(maf$Tumor_Sample_Barcode)),dimnames = list(unique(maf$Hugo_Symbol),unique(maf$Tumor_Sample_Barcode)))

for (i in colnames(mut.binary)) {
  tmp <- maf[which(maf$Tumor_Sample_Barcode == i),]
  if(is.element("Silent",tmp$Variant_Classification)) {
    tmp <- tmp[-which(tmp$Variant_Classification %in% "Silent"),]
  }
  for (j in tmp$Hugo_Symbol)
    mut.binary[j,i] <- 1
}
mut.binary <- as.data.frame(mut.binary); rownames(mut.binary) <- toupper(rownames(mut.binary))
mut.binary <- mut.binary[rowSums(mut.binary) > 0,]
mut.binary <- mut.binary[,which(substr(colnames(mut.binary), 14, 15) == "01")]

# load survival data
surv.info <- read.table("pancancerSurvivalData_XLu.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv.info <- as.data.frame(na.omit(surv.info[which(surv.info$type == "KIRC"),c("OS","OS.time")]))
surv.info <- surv.info[which(surv.info$OS.time > 0),] # get survival time greater than 0
colnames(surv.info) <- c("fustat","futime")
rownames(surv.info) <- paste0(rownames(surv.info),"-01")

cliquery <- cliquery[!duplicated(cliquery$Tumor_Sample_Barcode),]
rownames(cliquery) <- paste0(cliquery$Tumor_Sample_Barcode,"-01")
clin.info <- cliquery[,c("age_at_initial_pathologic_diagnosis",
                         "gender",
                         "race_list",
                         "stage_event_pathologic_stage",
                         "neoplasm_histologic_grade")]
colnames(clin.info) <- c("Age","Gender","Race","pStage","Grade")
write.table(clin.info,"TCGA_KIRC_Clinical_selected.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
clin.info <- read.table("TCGA_KIRC_Clinical_selected.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clin.info[which(clin.info$Race %in% c("AMERICAN INDIAN OR ALASKA NATIVE","BLACK OR AFRICAN AMERICAN")),"Race"] <- "Others"
# clin.info$Virus <- ifelse(clin.info$Virus != "","Affected","Clean")
# clin.info[which(clin.info$pStage %in% c("Stage IIIA","Stage IIIB","Stage IIIC")),"pStage"] <- "Stage III"
# clin.info[which(clin.info$pStage %in% c("Stage IVA","Stage IVB")),"pStage"] <- "Stage IV"
clin.info[clin.info == ""] <- NA

# get common samples
comsam <- intersect(colnames(tpm),colnames(mut.binary))
comsam <- intersect(comsam, colnames(meth))
comsam <- intersect(comsam, rownames(surv.info))
comsam <- intersect(comsam,rownames(clin.info))
#comsam <- intersect(comsam, colnames(mut))
rm.sam <- c("TCGA-BP-4758-01",
            "TCGA-B0-5705-01",
            "TCGA-A3-3313-01",
            "TCGA-A3-3374-01",
            "TCGA-AK-3433-01",
            "TCGA-AK-3440-01",
            "TCGA-AK-3456-01",
            "TCGA-AK-3465-01",
            "TCGA-AS-3777-01",
            "TCGA-B2-3923-01",
            "TCGA-BP-4334-01",
            "TCGA-B0-4699-01",
            "TCGA-BP-4756-01",
            "TCGA-BP-4760-01",
            "TCGA-BP-4769-01",
            "TCGA-BP-4784-01",
            "TCGA-BP-4795-01",
            "TCGA-BP-4994-01",
            "TCGA-BP-4995-01",
            "TCGA-B0-5083-01",
            "TCGA-B0-5117-01",
            "TCGA-B8-5546-01",
            "TCGA-CJ-5681-01")
comsam <- comsam[!(comsam %in% rm.sam)]

surv.info <- cbind.data.frame(surv.info[comsam,],clin.info[comsam,])
#surv.info <- cbind.data.frame(surv.info,TP53 = ifelse(as.numeric(mut.binary["TP53",comsam]) == 1,"Mutated","Wild"))

# output file
write.table(fpkm[Mids,comsam], "TCGA_KIRC_mRNA_FPKM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(fpkm[Lids,comsam], "TCGA_KIRC_lncRNA_FPKM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(tpm[Mids,comsam], "TCGA_KIRC_mRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(tpm[Lids,comsam], "TCGA_KIRC_lncRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(meth[,comsam], "TCGA_KIRC_Methpromoter.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(mut.binary[,comsam],"TCGA_KIRC_mut2.txt",sep = "\t",row.names = T,col.names = NA)
write.table(surv.info[comsam,],"TCGA_KIRC_surv.txt",sep = "\t",row.names = T,col.names = NA)

# filter elites
library(MOVICS)
mrna.tpm <- read.table("TCGA_KIRC_mRNA_TPM.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
lncrna.tpm <- read.table("TCGA_KIRC_lncRNA_TPM.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
meth <- read.table("TCGA_KIRC_Methpromoter.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
mut <- read.table("TCGA_KIRC_mut2.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
surv <- read.table("TCGA_KIRC_surv.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

mo.mrna <- getElites(dat = mrna.tpm,
                     method = "mad",
                     doLog2 = T,
                     lowpct = 0.1,
                     elite.num = 1500, # get top 1000 genes according mad values
                     scaleFlag = F,
                     centerFlag = F)

# if you want to use survival 
mo.mrna.surv <- getElites(dat = mrna.tpm,
                          surv.info = surv.info,
                          method = "cox",
                          doLog2 = T,
                          lowpct = 0.1,
                          elite.num = 1500, # this argument will be discard
                          p.cutoff = 0.05,
                          scaleFlag = F,
                          centerFlag = F)

mo.lncrna <- getElites(dat = lncrna.tpm,
                       method = "mad",
                       doLog2 = T,
                       lowpct = 0.1,
                       elite.num = 1000, # get top 1000 genes according mad values
                       scaleFlag = F,
                       centerFlag = F)

mo.meth <- getElites(dat = meth,
                     method = "mad",
                     doLog2 = F,
                     elite.num = 1500, # get top 1000 genes according mad values
                     scaleFlag = F,
                     centerFlag = F)

mo.mut <- getElites(dat = mut,
                    method = "freq",
                    doLog2 = F,
                    elite.pct = 0.03, # get mutations that have mutated rate greater than 3%
                    scaleFlag = F,
                    centerFlag = F)

# create input for MOVICS (very IMPORTANT!!!)
KIRC.tcga <- list(mRNA      = mo.mrna$elite.dat,
                  lncRNA    = mo.lncrna$elite.dat,
                  meth      = mo.meth$elite.dat,
                  mut       = mo.mut$elite.dat,
                  count     = count,
                  tpm       = tpm,
                  fpkm      = fpkm,
                  maf       = maf,
                  segment   = segment,
                  surv.info = surv)

# save RData
save(KIRC.tcga, file = "KIRC.tcga.RData")
save(KIRC.tcga, file = "../movics_pipeline/KIRC.tcga.RData")
save.image("data_preparing.RData")
