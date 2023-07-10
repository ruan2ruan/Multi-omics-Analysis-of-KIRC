#-----------------#
# MOVICS pipeline #

# set working path and creat directory
tumor.path <- "~/Ruan/KIRC/movics_pipeline"
setwd(tumor.path) #create dir
res.path    <- file.path(tumor.path, "Results")
fig.path    <- file.path(tumor.path, "Figures")
data.path  <- file.path(tumor.path, "Data")
comAnn.path <- file.path(tumor.path,"Annotation")
script.path <- file.path(tumor.path,"Scripts")
comRFun.path <- file.path(tumor.path,"commonFun")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }

# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")

mutect.dataframe <- function(x){
  # delete rows of Silent
  #cut_id <- x$Variant_Classification %in% c("Silent")
  sel_id <- x$Variant_Type == "SNP" #& x$Variant_Classification %in% c("Silent","Missense_Mutation","","Nonsense_Mutation","Nonstop_Mutation")
  #x <- x[!cut_id,]
  x <- x[sel_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}

indel.dataframe <- function(x){
  # delete rows of Silent
  #cut_id <- x$Variant_Classification %in% c("Silent")
  sel_id <- x$Variant_Type %in% c("DEL","INS")
  #x <- x[!cut_id,]
  x <- x[sel_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}

# load R package
library(MOVICS)
library(ggplot2)
library(RColorBrewer)

# load MOVICS input data
load("KIRC.tcga.RData")

# print name of example data
names(KIRC.tcga)
# [1] "mRNA"      "lncRNA"    "meth"      "mut"       "count"     "tpm"       "fpkm"      "maf"       "segment"  
# [10] "surv.info"       "segment"     "clin.info"

# extract multi-omics data
mo.data   <- KIRC.tcga[1:4]

# extract raw count data for downstream analyses
count     <- KIRC.tcga$count

# extract tpm data for downstream analyses
tpm       <- KIRC.tcga$tpm

# extract maf for downstream analysis
maf       <- KIRC.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- KIRC.tcga$segment

# extract survival information
surv.info <- KIRC.tcga$surv.info

#----------------------------------------------------------#
# 1. identify optimal clustering number (may take a while) #
optk.KIRC <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:6, # try cluster number from 2 to 6
                         fig.path    = fig.path,
                         fig.name    = "CLUSTER NUMBER OF TCGA-KIRC")
# take cluster number as 2 

# 2. perform multi-omics integrative clustering with the 10 algorithms #
moic.res.list <- getMOIC(data        = mo.data,
                             N.clust     = 2,
                             methodslist = list("iClusterBayes","SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"), 
                             type        = c("gaussian","gaussian","gaussian","binomial"))
save(moic.res.list, file = file.path(res.path,"moic.res.list.rda"))

#------------------------------------------------------#
# 3. get consensus clustering from different algorithm #
cmoic.KIRC <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               fig.path      = fig.path,
                               distance      = "pearson",
                               linkage       = "ward.D2")

getSilhouette(sil      = cmoic.KIRC$sil, # a sil object returned by getConsensusMOIC()
              fig.path = fig.path,
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

#------------------------------#
# 4. get comprehensive heatmap #
# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))

# set color for each omics data
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# extract sample annotation
annCol    <- surv.info[,c("Age","Race","Gender","pStage","Grade"), drop = FALSE]
annCol[is.na(annCol)] <- "Missing"
# annCol$Age <- as.numeric(annCol$Age)
# annCol$Smoke <- as.numeric(annCol$Smoke)
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age,na.rm = TRUE),
                                                           median(annCol$Age,na.rm = TRUE),
                                                           max(annCol$Age,na.rm = TRUE)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  pStage = c("Stage I"    = alpha(darkred,0.2),
                             "Stage II"   = alpha(darkred,0.4),
                             "Stage III"  = alpha(darkred,0.7),
                             "Stage IV"   = alpha(darkred,0.9), 
                             "Missing"    = "white"),
                  Grade = c("G1"    = alpha(darkgreen,0.2),
                             "G2"   = alpha(darkgreen,0.4),
                             "G3"  = alpha(darkgreen,0.7),
                             "G4"   = alpha(darkgreen,0.9), 
                             "GX"    = "white"),
                  Race   = c("ASIAN" = yellow,
                             "WHITE" = "grey80",
                             "Others" = brown,
                             "Missing" = "white"),
                  Gender = c("FEMALE" = jco[1],
                             "MALE"   = jco[2]))
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation
# library(impute)
# plotdata$meth <- impute.knn(plotdata$meth)$data
feat   <- moic.res.list[["CIMLR"]][["feat.res"]]
feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TPM","lncRNA.TPM","M value","Mutated"),
             clust.res     = cmoic.KIRC$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.path      = fig.path,
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")

#---------------------#
# 5. compare survival #
surv.KIRC <- compSurv(moic.res         = cmoic.KIRC,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      surv.cut         = 120,
                      p.adjust.method  = "none",
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

#------------------------------#
# 6. compare clinical features #
clust <- cmoic.KIRC$clust.res
clust$Subtype <- ifelse(clust$clust==1,"CS1","CS2")
tmp <- clust[order(clust$Subtype),]
tmp <- tmp[,"Subtype",drop=FALSE]
surv.info <- cbind(tmp,surv.info[rownames(tmp),])

clin.KIRC <- compClinvar(moic.res      = cmoic.KIRC,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("Race","pStage","fustat","Grade","Gender"), # features that are considered categorical variables
                         nonnormalVars = c("futime","Age"), # feature(s) that are considered using nonparametric test
                         exactVars     = c("Race","pStage","fustat","Grade","Gender"), # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         res.path      = res.path,
                         includeNA     = FALSE,
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
print(clin.KIRC$compTab)
write.csv(clin.KIRC$compTab,"clincompTab.csv")

library(tableone)
tab <- CreateTableOne(vars = c("Age","Gender","Race","pStage","Grade"), 
                       strata = "Subtype" , data = surv.info, 
                       factorVars = c("Race","Gender","pStage","Grade"))
tabMat <-print(tab, nonnormal = c("Age"), exact =c("Race","Gender","pStage","Grade"),showAllLevels = TRUE,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tabMat, file = "myTable.csv")

# customize the factor level for pstage
surv.info$pStage <- factor(surv.info$pStage, levels = c("Stage I","Stage II","Stage III","Stage IV"))
surv.info$Grade <- factor(surv.info$Grade, levels = c("G1","G2","G3","G4","GX"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.KIRC <- compAgree(moic.res  = cmoic.KIRC,
                        subt2comp = surv.info[,c("pStage","Grade")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")

#------------------------------------#
# 7. mutational frequency comparison #
mut.KIRC <- compMut(moic.res     = cmoic.KIRC,
                    mut.matrix   = KIRC.tcga$mut, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    #                    p.cutoff     = 0.05, # keep those genes with nominal p value < 0.25 to draw OncoPrint
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 1 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 8, 
                    height       = 10,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION",
                    res.path     = res.path,
                    fig.path     = fig.path)
print(mut.KIRC)

#----------------------------------#
# 7. compare total mutation burden #
tmb.KIRC <- compTMB(moic.res     = cmoic.KIRC,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

#------------------------------------#
# 8. compare fraction genome altered #
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")

# compare FGA, FGG, and FGL
fga.KIRC <- compFGA(moic.res     = cmoic.KIRC,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "BARPLOT OF FGA")

# compare agreement with other subtypes
subtype <- read.csv("subtype.csv")
rownames(subtype)<- subtype$PID
subtype$Gu <- ifelse(subtype$Gu=="KIRC_SdA_G2","G2","G1")
subtype$Gu <- factor(subtype$Gu, levels = c("G1","G2"))
subtype$Gu.new <- ifelse(subtype$Gu=="G2",2,1)

agree.Gu <- compAgree(moic.res  = cmoic.KIRC,
                      subt2comp = subtype[, "Gu",drop=FALSE],
                      doPlot    = TRUE,
                      fig.name  = "Agreement between the CMOIC subtype and Gu classification",
                      fig.path     = fig.path)

clust <- cmoic.KIRC$clust.res
clust$PID <- rownames(clust)
tmp <- merge(clust,subtype,by="PID",all=FALSE)
runKappa(subt1     = tmp$clust,
         subt2     = tmp$Gu.new,
         subt1.lab = "CMOIC",
         subt2.lab = "Gu",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and Gu",
         fig.path     = fig.path)

#--------------------------------#
# 9. drug sensitivity comparison #
drug.KIRC <- compDrugsen(moic.res    = cmoic.KIRC,
                         norm.expr   = tpm[,cmoic.KIRC$clust.res$samID], # double guarantee sample order
                         drugs       = c("Sorafenib","Sunitinib", "Axitinib", "Pazopanib"), # a vector of names of drug in GDSC (I get this from https://www.cancer.org/cancer/liver-cancer/treating/chemotherapy.html)
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         fig.path    = fig.path,
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 

drug.KIRC <- compDrugsen(moic.res    = cmoic.KIRC,
                         norm.expr   = tpm[,cmoic.KIRC$clust.res$samID], # double guarantee sample order
                         drugs       = c("Ruxolitinib"), # a vector of names of drug in GDSC (I get this from https://www.cancer.org/cancer/liver-cancer/treating/chemotherapy.html)
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         fig.path    = fig.path,
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50")

# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count, # raw count data
       moic.res   = cmoic.KIRC,
       res.path   = res.path,
       prefix     = "TCGA-KIRC") # prefix of figure name

#--------------------------------------------#
# 12. run biomarker identification procedure #
# choose limma result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.KIRC,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                       dat.path      = res.path, # path of DEA files
                       res.path      = res.path, # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.path      = fig.path,
                       width         = 12,
                       height        = 12,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.KIRC,
                       dea.method    = "deseq2",
                       dat.path      = res.path, # path of DEA files
                       res.path      = res.path, # path to save marker files
                       prefix        = "TCGA-KIRC",
                       dirct         = "down",
                       n.marker      = 100, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = tpm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.path      = fig.path,
                       width         = 12,
                       height        = 12,
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")

#################### GSEA ####################
library(clusterProfiler)
library(enrichplot)

MSigDB=read.gmt(file.path(comAnn.path,"h.all.v7.4.symbols.gmt"))
MSigDB$term <- str_replace(MSigDB$term, "HALLMARK_","")
MSigDB$term <- str_replace_all(MSigDB$term, "_"," ")

CS12 <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- CS12$log2fc
names(geneList) <- rownames(CS12)
geneList <- sort(geneList,decreasing = T)
GSEA.CS12 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F,pvalueCutoff = 0.25)
pdf(file.path(fig.path,"GSEA.CS12.h_bubble.pdf"),width =8,height = 9)
dotplot(GSEA.CS12, showCategory=15)
dev.off()

pdf(file.path(fig.path,"GSEA enrichplot.pdf"),width = 8.5,height = 5)
gseaplot2(GSEA.CS12, geneSetID = c(1,4,8,13,14),base_size=12,
          color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
          rel_heights = c(1.7, 0.5, 0.8))
dev.off()
# pdf(file.path(fig.path,"GSEA KRAS enrichplot.pdf"),width = 8,height = 6)
# gseaplot2(GSEA.CS12, geneSetID = c(13),base_size=9,
#           color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
#           rel_heights = c(1.7, 0.5, 0.8))
# dev.off()

#################### GO ####################
tmp <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
gene_up <- rownames(tmp[which(tmp$padj < 0.05 & tmp$pvalue < 0.05  & tmp$log2fc > log2(2)),])
#gene <- rownames(tmp[which(tmp$padj < 0.25 & tmp$pvalue < 0.05  & tmp$log2FoldChange > log2(1.5)),])
tmp <- bitr(as.character(gene_up), fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
ego.up <-     enrichGO(gene          = as.character(tmp$ENSEMBL),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENSEMBL",
                       ont           = "ALL",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.25,
                       minGSSize     = 10,
                       readable      = T)
#ego2.up <- simplify(ego.up)
write.table(as.data.frame(summary(ego.up)),file.path(res.path,"Upregualted genes in CS1_vs_CS2 GO all clusterprofiler.txt"),sep = "\t",row.names = F)

tmp <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
gene_dn <- rownames(tmp[which(tmp$padj < 0.05 & tmp$pvalue < 0.05  & tmp$log2fc < -log2(2)),])
#gene <- rownames(tmp[which(tmp$padj < 0.25 & tmp$pvalue < 0.05  & tmp$log2FoldChange < -log2(1.5)),])
tmp <- bitr(as.character(gene_dn), fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
ego.dn <-     enrichGO(gene          = as.character(tmp$ENSEMBL),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENSEMBL",
                       ont           = "ALL",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.25,
                       minGSSize     = 10,
                       readable      = T)
#ego2 <- simplify(ego)
write.table(data.frame(summary(ego.dn)),file.path(res.path,"Dnregualted genes in CS1_vs_CS2 GO all clusterprofiler.txt"),sep = "\t",row.names = F)

library(cowplot)
pdf(file.path(fig.path,"CS1_vs_CS2.GO_bubble_top15.pdf"),width = 15,height = 10)
p1 <- dotplot(ego.up, showCategory=15) + ggtitle("CS1 vs CS2 dotplot for GO UP")
p2 <- dotplot(ego.dn, showCategory=15) + ggtitle("CS1 vs CS2 dotplot for GO DOWN")
plot_grid(p1, p2, ncol=2)
dev.off()

#################### KEGG ######################
tmp <- data.frame(SYMBOL=rownames(CS12),logFC=CS12$log2fc)
gsym.id <- bitr(tmp$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gsym.fc.id <- merge(tmp, gsym.id, by="SYMBOL", all=F)
head(gsym.fc.id)
GSEA_input <- gsym.fc.id$logFC
names(GSEA_input) = gsym.fc.id$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
head(GSEA_input)
kk.polyps <- gseKEGG(GSEA_input, organism = 'hsa',pvalueCutoff = 0.25)
sortkk <- kk.polyps[order(kk.polyps$enrichmentScore, decreasing = T),]
colnames(sortkk)[1] <- "PATHID"
write.table(as.data.frame(sortkk),file.path(res.path,"KEGG results for CS1 VS CS2.txt"),sep = "\t",row.names = F,quote = F)

pdf(file.path(fig.path," CS1 VS CS2 KEGG_bubble_top15.pdf"),width = 8,height = 8)
dotplot(kk.polyps, showCategory=15)
dev.off()

pdf(file.path(fig.path,"KEGG enrichplot.pdf"),width = 8.5,height = 5)
gseaplot2(kk.polyps, geneSetID = c(12,14,24,32),base_size=13,
          color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
          rel_heights = c(1.7, 0.5, 0.8))
dev.off()

pdf(file.path(fig.path,"KEGG enrichplot (PPAR).pdf"),width = 8.5,height = 5)
gseaplot2(kk.polyps, geneSetID = c(16),base_size=12,
          color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
          rel_heights = c(1.7, 0.5, 0.8))
dev.off()


TPM <- read.table(file.path(res.path,"TCGA_KIRC_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
gene <- c("PGAM1","ENO1","LDHA","TPI1","P4HA1","MRPS17","ADM",
                       "NDRG1","TUBB6","ALDOA","MIF","SLC2A1","CDKN3","ACOT7")
Hypoxia_Signature <- TPM[gene,rownames(annCol_new)]
choose_matrix = standarize.fun(indata=Hypoxia_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 

annCol1 <- annCol_new[,"Subtype",drop=FALSE]
annColors1 <- list(Subtype = c("CS1" = "#2EC4B6", "CS2" = "#E71D36"))

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = T,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_KIRC_15gene_Hypoxia_Signatures.pdf"),
         border_color = NA,fontsize_row = 6, fontsize_col = 5,cellwidth = 1.3,cellheight = 8)

###############Volcano##############
library( "ggplot2" )

x <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
x$label<- rownames(x)
head(x)
#colnames(x) <- c("baseMean","logFC","lfcSE","stat","pvalue","padj","Genesymbol","label")
colnames(x) <- c("baseMean","logFC","pvalue","padj","label")
#plot_mode <- "classic" 
plot_mode <- "advanced"

logFCcut <- log2(2) 
logFCcut2 <- log2(4) #for advanced mode

adjPcut <- 0.25
adjPcut2 <- 0.05 #for advanced mode

DEG_up <- x[x$log2FoldChange > log2(2) & x$padj<0.05,]
DEG_DN <- x[x$log2FoldChange < -log2(2) & x$padj<0.05,]

xmin <- -4
xmax <- 4
ymin <- 0
ymax <- 6

if (plot_mode == "classic"){
  # setting for color
  x$color_transparent <- ifelse((x$padj < adjPcut & x$logFC > logFCcut), "red", ifelse((x$padj < adjPcut & x$logFC < -logFCcut), "blue","grey30"))
  # setting for size
  size <- ifelse((x$padj < adjPcut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  cols[x$padj < adjPcut & x$logFC >logFCcut]<- "#EBE645"
  cols[x$padj < adjPcut2 & x$logFC > logFCcut2]<- "#EBE645"
  cols[x$padj < adjPcut & x$logFC < -logFCcut]<- "#577BC1"
  cols[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- "#577BC1"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  size[x$padj < adjPcut & x$logFC > logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC > logFCcut2]<- 5
  size[x$padj < adjPcut & x$logFC < -logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- 5
  
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(padj), label = label, color = pathway)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~log[2]~"(FoldChange)"), y=bquote(~-log[10]~"FDR"), title="") + 
  scale_y_continuous(
    breaks = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5), 
    labels = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5),
    limits = c(ymin,ymax)
  ) +
  scale_x_continuous(
    breaks = c(-4,-3, -round(logFCcut2,1),-round(logFCcut,1), 0, round(logFCcut,1),round(logFCcut2,1),3, 4),
    labels = c(-4,-3, -round(logFCcut2,1),"-1", 0, "1",round(logFCcut2,1),3, 4),
    limits = c(-4, 4) 
  ) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = -log10(adjPcut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times"
  ) +
  theme(panel.grid=element_blank()) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=11),
        axis.text.y = element_text(face="bold",  color="black", size=11),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11))

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(adjPcut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}

ggsave(file.path(fig.path,"volcano plot for DEGs between CS1 and CS2.pdf"),width = 8,height = 6)

#################### GSVA ######################
library( "pheatmap" )
library(GSVA)
library("gplots")
library(stringr)

TPM.HUGO <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
MSigDB=read.gmt(file.path(comAnn.path,"h.all.v7.4.symbols.gmt"))
met_sig <- MSigDB
met_sig$term <- str_replace_all(string = met_sig$term,pattern = "HALLMARK_",replacement = "")
met_sig$term <- str_replace_all(string = met_sig$term,pattern = "_",replacement = " ")
label <- unique(met_sig$term)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$term == i),"gene"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for hallmark by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap
Hallmark_Signature <- read.table(file.path(res.path,"gsva results for hallmark by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Hallmark_Signature <- Hallmark_Signature[,rownames(surv.info)]
choose_matrix = standarize.fun(indata=Hallmark_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 

annCol1 <- surv.info[,"Subtype",drop=FALSE]
annColors1 <- list(Subtype = c("CS1" = "#2EC4B6", "CS2" = "#E71D36"))

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = T,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_KIRC_gsva_Hallmark_Signatures.pdf"),
         border_color = NA,fontsize_row = 6, fontsize_col = 5,cellwidth = 1.3,cellheight = 8)

#################### Immune genesets ####################
### CCR Curated Immune Cell Signature ###
TPM.HUGO <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`CellType`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`CellType` == i),"Symbol"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for CCR_Curated_Immune_Cell_Signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")
# heatmap
Immune_Cell_Signature <- read.table(file.path(res.path,"gsva results for CCR_Curated_Immune_Cell_Signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Immune_Cell_Signature <- Immune_Cell_Signature[,rownames(annCol1)]
choose_matrix = standarize.fun(indata=Immune_Cell_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_KIRC_gsva_CCR_Curated_Immune_Cell_Signature.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.2,cellheight = 12)
#dev.off()

############## Oncogenetic signature
# TPM.HUGO <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
# met_sig <- read.table(file.path(comAnn.path,"10 Oncogenetic_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
# label <- unique(met_sig$`Pathway`)
# meta_sig <- list()
# for (i in label) {
#   meta_sig[[i]] <- met_sig[which(met_sig$`Pathway` == i),"Symbol"]  
# }
# meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
# write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for 10 Oncogenetic_signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")
# # heatmap
# Oncogenetic_Signature <- read.table(file.path(res.path,"gsva results for 10 Oncogenetic_signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
# Oncogenetic_Signature <- Oncogenetic_Signature[,rownames(annCol1)]
# choose_matrix = standarize.fun(indata=Oncogenetic_Signature,halfwidth = 1.5)
# mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 
# 
# pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
#          annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
#          annotation_legend = T, filename = file.path(fig.path,"heatmap_KIRC_gsva_10 Oncogenetic_Signature.pdf"),
#          border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.5,cellheight = 20)
# #dev.off()

############## Some signatures (JAK-STAT)
TPM.HUGO <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"some signatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`Pathway`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`Pathway` == i),"Symbol"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for some signatures by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap
JAK_Signature <- read.table(file.path(res.path,"gsva results for some signatures by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
JAK_Signature <- JAK_Signature[,rownames(annCol1)]
choose_matrix = standarize.fun(indata=JAK_Signature,halfwidth = 1.5)
rownames(choose_matrix) <- str_replace_all(string = rownames(choose_matrix),pattern = "_",replacement = " ")
mycol <- colorpanel(256,low="#3f88c5",mid = "black",high="#ffd424") 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1[,1,drop=FALSE], annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_KIRC_gsva_JAK_signatures.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 2,cellheight = 16)
#dev.off()
p.JAK <- c()
for (i in 1:9) {
  tmp1 <- as.numeric(JAK_Signature[i,CS1])
  tmp2 <- as.numeric(JAK_Signature[i,CS2])
  p.JAK <- c(p.JAK,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.JAK) <- rownames(JAK_Signature)

########################## C1 VS C2 Immune #############################
### CIBERSORT ###
CS1 <- rownames(annCol1[which(annCol1$Subtype=="CS1"),,drop=FALSE])
CS2 <- rownames(annCol1[which(annCol1$Subtype=="CS2"),,drop=FALSE])

TPM <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp <- log2(TPM + 1)
tmp <- sweep(tmp,1,apply(tmp, 1, median))
write.table(tmp,file.path(res.path,"CIBERSOT.TPM.logTransMedCent.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- read.table(file.path(res.path,"CIBERSORT.Output_Job3.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
tmp <- as.data.frame(t(tmp[,1:22]))
cibersort <- tmp[,rownames(annCol)]
p <- c()
for (i in 1:22) {
  tmp1 <- as.numeric(tmp[i,CS1])
  tmp2 <- as.numeric(tmp[i,CS2])
  p <- c(p,wilcox.test(tmp1,tmp2)$p.value)
}
names(p) <- rownames(tmp)
cibersort.p<- as.data.frame(p)

### TIDE ###
tmp <- log2(TPM + 1)
tmp <- sweep(tmp,2,apply(tmp, 2, median))
tmp <- sweep(tmp,1,apply(tmp, 1, median))
write.table(tmp,file.path(res.path,"TIDE.TPM.logTransMedCentBoth.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

TIDE <- read.csv(file.path(res.path,"KIRC_TIDE_output.csv"),check.names = F, header = T,stringsAsFactors = F,row.names = 1)
TIDE <- TIDE[rownames(annCol1),]
TIDE$Subtype <- rep(c("CS1","CS2"),
                    c(length(CS1),length(CS2)))
chisq.test(table(TIDE$Responder,TIDE$Subtype))  # X-squared = 20.307, df = 1, p-value = 6.595e-06
table(TIDE$Responder,TIDE$Subtype)

TIDE_result <- data.frame(table(TIDE$Responder,TIDE$Subtype)) 
TIDE_result<- TIDE_result %>% 
  group_by(Var2) %>% 
  mutate(sumVal = sum(Freq)) %>%
  ungroup() %>%
  mutate(per=Freq/sumVal) %>%
  mutate(label=paste0(round(per*100,2),"%"))%>%
  mutate(Var1=ifelse(TIDE_result$Var1=="False","Non-Response","Response"))

pdf(file.path(fig.path,"TIDE plot.pdf"),width = 8,height = 6)
ggplot(aes(x=Var2, y=per, fill=Var1), data=TIDE_result)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#119da4","#ffc857"))+
  xlab("")+ 
  ylab("Percentage")+
  geom_text(aes(label = label), 
            color="black", size=5,position=position_fill(0.5))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw() +
  theme(axis.title.y=element_text(size=18),
        axis.text=element_text(size=13),
        legend.text=element_text(size=13))
dev.off()

### MCPCOUNTER and ESITIMATE ###
tmp <- log2(TPM + 1)

# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
filterCommonGenes(input.f=file.path(res.path, "TCGA_KIRC_mRNA_TPM.txt") , output.f=file.path(res.path,"KIRC_HUGO_symbol_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"KIRC_HUGO_symbol_ESTIMATE.txt"), file.path(res.path,"KIRC_HUGO_symbol_estimate_score.txt"), platform="affymetrix")
est <- read.table(file = file.path(res.path,"KIRC_HUGO_symbol_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est) <- est[,2]; colnames(est) <- est[1,]; est <- est[-1,c(-1,-2)];
est <- sapply(est, as.numeric); rownames(est) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")
est.raw <- est; colnames(est.raw) <- colnames(TPM)
source(file.path(script.path,"annTrackScale.R"))
tmp <- annTrackScale(indata = est, halfwidth = 2, poolsd = F); tmp <- as.data.frame(t(tmp))
rownames(tmp) <- colnames(TPM)
annCol1 <- cbind.data.frame(annCol1,tmp[rownames(annCol1),])

library(MCPcounter)
MCPscore <- MCPcounter.estimate(expression = tmp,featuresType = "HUGO_symbols")
write.table(MCPscore,"KIRC_TPM_MCPscore_HUGO_symbols.txt",row.names=T, col.names=NA, sep="\t", quote=F)
MCPscore <- read.table(file.path(res.path,"KIRC_TPM_MCPscore_HUGO_symbols.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)

tmp <- as.data.frame(MCPscore[,rownames(annCol1)])
MCPscore.raw <- tmp
tmp <- annTrackScale(indata = tmp, halfwidth = 2, poolsd = F)
tmp <- as.data.frame(t(tmp))
annCol1 <- cbind.data.frame(annCol1,tmp)

p.mcp <- c()
for (i in 1:10) {
  tmp1 <- as.numeric(MCPscore.raw[i,CS1])
  tmp2 <- as.numeric(MCPscore.raw[i,CS2])
  p.mcp <- c(p.mcp,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.mcp) <- rownames(MCPscore.raw)

p.est <- c()
for (i in 1:4) {
  tmp1 <- as.numeric(est.raw[i,CS1])
  tmp2 <- as.numeric(est.raw[i,CS2])
  p.est <- c(p.est,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.est) <- rownames(est.raw)

# boxplot for scores
require(tidyr)
annCol_new <- annCol1[which(annCol1$Subtype=="CS1"|annCol1$Subtype=="CS2"),]

dd <- as.data.frame(est.raw[,c(CS1,CS2)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:252)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS1","CS2"),c(length(CS1),length(CS2))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(est.raw))
d2$Subtype <- factor(d2$Subtype,levels = c("CS1","CS2"))
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
d2_1 <- d2[d2$cell %in% c("StromalScore","ImmuneScore","ESTIMATEScore"),]
d2_2 <- d2[!(d2$cell %in% c("StromalScore","ImmuneScore","ESTIMATEScore")),]
pv_1 <- pv[pv$cell %in% c("StromalScore","ImmuneScore","ESTIMATEScore"),]
pv_2 <- pv[!(pv$cell %in% c("StromalScore","ImmuneScore","ESTIMATEScore")),]
p1 <- ggplot(d2_1, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2_1$ES) * 1.1, 
                label=pv_1$sigcode),
            data=pv_1, inherit.aes=F) + 
  xlab(NULL)+ylab("ESITIMATE Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank()) 
ggsave(file.path(fig.path,"Boxplot for ESITIMATE score of CS1 and CS2(3).pdf"),width = 8,height = 3)
p2 <- ggplot(d2_2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2_2$ES) * 1.1, 
                label=pv_2$sigcode),
            data=pv_2, inherit.aes=F) + 
  xlab(NULL)+ylab("ESITIMATE Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank()) 
ggsave(file.path(fig.path,"Boxplot for ESITIMATE score of CS1 and CS2(1).pdf"),width = 5,height = 3)

dd <- as.data.frame(MCPscore.raw[,c(CS1,CS2)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:252)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS1","CS2"),c(length(CS1),length(CS2))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(MCPscore.raw))
d2$Subtype <- factor(d2$Subtype,levels = c("CS1","CS2"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
p1 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("MCPcounter Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for MCPcounter score of CS1 and CS2.pdf"),width = 8,height = 3)

dd <- as.data.frame(cibersort[,c(CS1,CS2)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:252)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS1","CS2"),
                                                          c(length(CS1),length(CS2))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(cibersort))
d2$Subtype <- factor(d2$Subtype,levels = c("CS1","CS2"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
library(ggplot2)
p2 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("CIBERSORT Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for CIBERSORT of CS1 and CS2.pdf"),width = 8,height = 3)

########################## CIBERSORTx
tmp <- read.table(file.path(res.path,"CIBERSORTx_Job2_Adjusted.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
tmp <-tmp %>%
  dplyr::select("B cells","CD4+ Effector","CD4+ Naive","CD4+ Proliferating","CD4+ Activated IEG","CD4+ Treg",
         "CD8A+ Exhausted","CD8A+ NK-like","CD8A+ Proliferating","CD8A+ Tissue-resident",
         "NK cells","Macrophages","Endothelial","Monocytes","Dendritic cells","Mast")
tmp <- as.data.frame(t(tmp))

cibersort <- tmp[,rownames(annCol_new)]
  
p <- c()
for (i in 1:16) {
  tmp1 <- as.numeric(tmp[i,CS1])
  tmp2 <- as.numeric(tmp[i,CS2])
  p <- c(p,wilcox.test(tmp1,tmp2)$p.value)
}
names(p) <- rownames(tmp)
cibersort.p<- as.data.frame(p)

library(tidyr)
dd <- as.data.frame(cibersort[,c(CS1,CS2)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:252)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS1","CS2"),
                                                       c(length(CS1),length(CS2))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(cibersort))
d2$Subtype <- factor(d2$Subtype,levels = c("CS1","CS2"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
library(ggplot2)
ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("CIBERSORT Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for CIBERSORTx of CS1 and CS2.pdf"),width = 8,height = 3)


tmp <- read.table(file.path(res.path,"ImmuCellAI_abundance_result.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
tmp <-tmp[,1:24]
tmp <- as.data.frame(t(tmp))

cibersort <- tmp[,rownames(annCol_new)]

p <- c()
for (i in 1:24) {
  tmp1 <- as.numeric(tmp[i,CS1])
  tmp2 <- as.numeric(tmp[i,CS2])
  p <- c(p,wilcox.test(tmp1,tmp2)$p.value)
}
names(p) <- rownames(tmp)
cibersort.p<- as.data.frame(p)

library(tidyr)
dd <- as.data.frame(cibersort[,c(CS1,CS2)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:252)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS1","CS2"),
                                                          c(length(CS1),length(CS2))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(cibersort))
d2$Subtype <- factor(d2$Subtype,levels = c("CS1","CS2"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
library(ggplot2)
ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("CIBERSORT Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for CIBERSORTx of CS1 and CS2.pdf"),width = 8,height = 3)

########################## GSVA #############################
### EMT(NTP) ###
library(CMScaller)
EMT_signature <- read.table(file.path(comAnn.path,"EMT_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
emat <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
ntp.matrix <- ntp(emat=emat, templates=EMT_signature, nPerm=1000,distance = "pearson", doPlot =  FALSE)
ntp <- ntp.matrix[rownames(annCol1),]
ntp$Cluster <- annCol1$Subtype
#write.table(ntp,file.path(res.path,"KIRC_EMT_result.txt"),row.names = T,col.names = NA,sep = "\t")
table(ntp$prediction,ntp$Cluster)
fisher.test(table(ntp$prediction,ntp$Cluster))   #p-value = 0.1448

################## SubMap ###################
source(file.path(comRFun.path,"generateInputFileForSubMap.R"))

gct_file=file.path(res.path,'Submap_gct.gct')
cls_file=file.path(res.path,'Submap_cls.cls')

sam_info=data.frame('rank'=annCol_new$Subtype)
rownames(sam_info)=rownames(annCol_new)
sam_info$type=sam_info$rank
data <- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
in_gct=log2(data+1)
in_gct <- in_gct[,rownames(sam_info)]
generateInputFileForSubMap(in_gct, gct_file, cls_file, sam_info, type_name = "type")

submap <- read.table(file.path(res.path,"SubMapResult.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F)
rownames(submap) <- submap$Subtype
library( "pheatmap" )
choose_matrix <- submap[,2:5]
annotation_row = data.frame(Pvalue = factor(rep(c("Nominal.p", "Bonferroni"), c(2,2))))
rownames(annotation_row) = rownames(choose_matrix)
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)
mycol <- c("#e5361f","#d16c69","#5587ec","#7968cc","#423293")
library("RColorBrewer")
#brewer.pal(8,"Set2")
ann_colors = list( Pvalue = c(Nominal.p = "#61bb30", Bonferroni = gold))
pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F, annotation_row=annotation_row,
         annotation_colors=ann_colors , show_rownames = T,show_colnames = T,
         gaps_row = 2, annotation_legend = T, legend_breaks=c(0,0.25,0.5,0.75,1),
         filename = file.path(fig.path,"heatmap_submap.pdf"),
         border_color = "white", height=4)
dev.off()

############################ Step Logistic ############################
DE <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
DE <- DE %>%
  filter(pvalue<0.05 & (log2fc>5|log2fc<(-3)))
TPM.norm <- log2(TPM.HUGO + 1)
TPM.norm <- data.frame(t(na.omit(TPM.norm[rownames(DE),])))

tmp <-TPM.norm[rownames(surv.info),]
survival_info_df <- cbind(tmp,surv.info)[,1:19]
survival_info_df$Subtype <- as.factor(survival_info_df$Subtype)

set.seed(12345)
nn=0.7
sub<-sample(1:nrow(survival_info_df),round(nrow(survival_info_df)*nn))
data_train<-survival_info_df[sub,]
data_test<-survival_info_df[-sub,]

model <- glm(Subtype ~., data = data_train, family = binomial)
#summary(model)
step.model=step(model,direction = "both") 
summary(step.model)

# HP+KCNS1+ANGPTL8+SCGN+TMEM174
model.new <- glm(Subtype ~ HP+KCNS1+ANGPTL8+SCGN+TMEM174,
                 data = data_train, family = binomial)
summary(model.new)
model.HP <- glm(Subtype ~ HP,data = data_train, family = binomial)
model.KCNS1 <- glm(Subtype ~ KCNS1, data = data_train, family = binomial)
model.ANGPTL8 <- glm(Subtype ~ ANGPTL8, data = data_train, family = binomial)
model.SCGN <- glm(Subtype ~ SCGN, data = data_train, family = binomial)
model.TMEM174 <- glm(Subtype ~ TMEM174, data = data_train, family = binomial)

genes <- c("HP","KCNS1","ANGPTL8","SCGN","TMEM174")
data_test<-data_test[,genes]

library(pROC)
pred <- predict(model.new,data_test, type="response")
pred.HP <- predict(model.HP,data_test[,"HP",drop=FALSE], type="response")
pred.KCNS1 <- predict(model.KCNS1,data_test[,"KCNS1",drop=FALSE], type="response")
pred.ANGPTL8 <- predict(model.ANGPTL8,data_test[,"ANGPTL8",drop=FALSE], type="response")
pred.SCGN <- predict(model.SCGN,data_test[,"SCGN",drop=FALSE], type="response")
pred.TMEM174 <- predict(model.TMEM174,data_test[,"TMEM174",drop=FALSE], type="response")
pred_class <- ifelse(pred > 0.5, "CS2","CS1")
pred_class<-factor(pred_class,levels = c("CS1","CS2"),order=TRUE)
# confusion <- table(tmp$Subtype, pred_class)
# confusion
# library(caret)
# confusionMatrix(tmp$Subtype, pred_class)

tmp <-survival_info_df[-sub,]

pdf(file.path(fig.path,"5 Sig ROC.pdf"),width = 8,height = 8)
rocobj1 <- plot.roc(tmp$Subtype, pred,
                    main="Statistical comparison",
                    percent=TRUE,print.auc=TRUE,
                    col="#f94144",lwd=1)
rocobj2 <- lines.roc(tmp$Subtype, pred.HP,
                     percent=TRUE,
                     col="#f9844a",lwd=1)
rocobj3 <- lines.roc(tmp$Subtype, pred.KCNS1,
                     percent=TRUE,
                     col="#f9c74f",lwd=1)
rocobj4 <- lines.roc(tmp$Subtype, pred.ANGPTL8,
                     percent=TRUE,
                     col="#90be6d",lwd=1)
rocobj5 <- lines.roc(tmp$Subtype, pred.SCGN,
                     percent=TRUE,
                     col="#577590",lwd=1)
rocobj6 <- lines.roc(tmp$Subtype, pred.TMEM174,
                     percent=TRUE,
                     col="#277da1",lwd=1)
legend("bottomright", legend=c("5 Sig", "HP","KCNS1","ANGPTL8","SCGN","TMEM174"),
       col=c("#f94144","#f9844a","#f9c74f","#90be6d","#577590","#277da1"), lwd=2,cex=0.8)
dev.off()

# KCNS1+SCGN+TMEM174
model.new <- glm(Subtype ~ KCNS1+SCGN+TMEM174, 
                 data = data_train, family = binomial)
summary(model.new)
model.KCNS1 <- glm(Subtype ~ KCNS1, data = data_train, family = binomial)
model.SCGN <- glm(Subtype ~ SCGN, data = data_train, family = binomial)
model.TMEM174 <- glm(Subtype ~ TMEM174, data = data_train, family = binomial)

genes <- c("KCNS1","SCGN","TMEM174")
data_test<-data_test[,genes]

library(pROC)
pred <- predict(model.new,data_test, type="response")
pred.KCNS1 <- predict(model.KCNS1,data_test[,"KCNS1",drop=FALSE], type="response")
pred.SCGN <- predict(model.SCGN,data_test[,"SCGN",drop=FALSE], type="response")
pred.TMEM174 <- predict(model.TMEM174,data_test[,"TMEM174",drop=FALSE], type="response")
pred_class <- ifelse(pred > 0.5, "CS2","CS1")
pred_class<-factor(pred_class,levels = c("CS1","CS2"),order=TRUE)
# confusion <- table(tmp$Subtype, pred_class)
# confusion
# library(caret)
# confusionMatrix(tmp$Subtype, pred_class)

tmp <-survival_info_df[-sub,]

pdf(file.path(fig.path,"3 Sig ROC.pdf"),width = 8,height = 8)
rocobj1 <- plot.roc(tmp$Subtype, pred,
                    main="Statistical comparison",
                    percent=TRUE,print.auc=TRUE,
                    col="#f94144",lwd=1)
rocobj2 <- lines.roc(tmp$Subtype, pred.KCNS1, 
                     percent=TRUE, 
                     col="#f9c74f",lwd=1)
rocobj3 <- lines.roc(tmp$Subtype, pred.SCGN, 
                     percent=TRUE, 
                     col="#577590",lwd=1)
rocobj4 <- lines.roc(tmp$Subtype, pred.TMEM174, 
                     percent=TRUE, 
                     col="#277da1",lwd=1)
legend("bottomright", legend=c("3 gene Sig","KCNS1","SCGN","TMEM174"),
       col=c("#f94144","#f9c74f","#577590","#277da1"), lwd=2,cex=0.8)
dev.off()

############################ COX ############################
library('survival')
library('survminer')
library('survivalROC')
library('plotROC')

DE <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
DE <- DE %>%
  filter(pvalue<0.05 & log2fc>3.5)
TPM.norm <- log2(TPM.HUGO + 1)
TPM.norm <- data.frame(t(na.omit(TPM.norm[rownames(DE),])))

tmp <-TPM.norm[rownames(surv.info),]
survival_info_df <- cbind(tmp,surv.info)
#rowname <-gsub(rownames(DE),pattern='-',replacement ='_')

# filter potential useful sig genes using univariate cox regression.
uni_cox_in_bulk<-function(gene_list, survival_info_df){
  library('survival')
#  gene_list<-gsub(gene_list, pattern='-', replacement ='_')
  uni_cox<-function(single_gene){
    formula<-as.formula(paste0('Surv(futime, fustat)~',single_gene))
    surv_uni_cox<-summary(coxph(formula, data=survival_info_df))
    ph_hypothesis_p<-cox.zph(coxph(formula, data=survival_info_df))$table[1,3]
    if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){ # get the pvalue 
      single_cox_report<-data.frame('uni_cox_sig_genes'=single_gene,
                                    'beta'=surv_uni_cox$coefficients[,1], 
                                    'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]), 
                                    'z_pvalue'=surv_uni_cox$coefficients[,5], 
                                    'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                    'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report 
    }
  }
  uni_cox_list<-lapply(gene_list, uni_cox)
  do.call(rbind, uni_cox_list)
}
uni_cox_df<-uni_cox_in_bulk(gene_list=colnames(TPM.norm), survival_info_df=survival_info_df)

mySurv<-survival_info_df[,c("futime","fustat",uni_cox_df$uni_cox_sig_genes)]
# sig_gene_multi_cox <- uni_cox_df$uni_cox_sig_genes
# formula_for_multivariate <- as.formula(paste0('Surv(futime, fustat)~',paste(sig_gene_multi_cox,sep='',collapse = '+')))
# multi_variate_cox <- coxph(formula_for_multivariate,data=survival_info_df)
# ph_hypo_multi <- cox.zph(multi_variate_cox)
# ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
# formula_for_multivariate <- as.formula(paste0('Surv(futime, fustat)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep='',collapse = '+')))
# multi_variate_cox2 <- coxph(formula_for_multivariate,data=survival_info_df)
# multi_variate_cox

multi_COX<-coxph(Surv(futime, fustat) ~ ., data=mySurv)  
summary(multi_COX)#看一下结果COX回归分析的结果
step.multi_COX=step(multi_COX,direction = "both")  #方法可以选择"both"（两边，这是默认的方法）,"backward"（向后逐步回归）,"forward"（向前逐步回归）
step.multi_COX

write.csv(uni_cox_df,"uni_cox_df.csv")
#write.csv(step.multi_COX,"step.multi_COX.csv")

pdf(file.path(fig.path,"Multi COX Hazard ratios of candidate genes.pdf"),width = 12,height = 8)
ggforest(model = step.multi_COX,data=mySurv,main = "Hazard ratios of candidate genes",fontsize = 1)
dev.off()

# genes <- c("ANGPTL8","PKP3","TRPM3","PI3","KRTAP2.3","SLC5A1","HOXB13","SLC6A19","CYP4A11","SLC22A6","PLG")
# 
# TPM_RiskScore <- TPM.norm[,genes]
# TPM_RiskScore<-TPM_RiskScore %>%
#   mutate(RiskScore=0.15980*ANGPTL8+0.46967*PKP3+0.70215*TRPM3-0.25178*PI3+
#            0.46607*KRTAP2.3-0.14873*SLC5A1+0.31296*HOXB13-0.20905*SLC6A19+0.21071*CYP4A11-
#            0.33183*SLC22A6-0.17567*PLG) %>%
#   mutate(risk_group=ifelse(RiskScore>=median(RiskScore),'HRisk','LRisk'))
#   
# TPM_RiskScore <-TPM_RiskScore[rownames(surv.info),c("RiskScore","risk_group")]
# survival_risk <- cbind(TPM_RiskScore,surv.info)
# 
# fit<- survfit(Surv(round(futime/30.5,4), fustat) ~ risk_group, data = survival_risk,type = "kaplan-meier", error = "greenwood", conf.type = "plain", 
#               na.action = na.exclude)
# names(fit$strata) <- gsub("Subtype=", "", names(fit$strata))
# ggsurvplot(fit, data = survival_risk,
#            surv.median.line = "h", # Add medians survival
#            
#            # Change legends: title & labels
#            legend.title = "RiskGroup",
#            legend.labs = c("HRisk", "LRisk"),
#            # Add p-value and tervals
#            pval = TRUE,
#            
#            conf.int = FALSE,
#            # Add risk table
#            risk.table = TRUE,
#            tables.height = 0.2,
#            risk.table.col = "strata",
#            risk.table.y.text = FALSE,
#  #          tables.theme = theme_cleantable(),
#            
#            # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
#            # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
#            palette = c("#E7B800", "#2E9FDF"),
#            ggtheme = theme_classic() # Change ggplot2 theme
# )

###### TCGA(N=255)
TPM.ALL <- read.table(file.path(res.path,"TCGA_KIRC_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
colnames <- colnames(TPM.ALL)
colnames <- colnames[!(colnames %in% rm.sam)]
colnames <- colnames[!(colnames %in% colnames(TPM.HUGO))]
TPM.ALL <- TPM.ALL[,colnames]
  
genes <- c("ANGPTL8","IL20RB","SAA1","WISP2","PI3","CCNA1","SYT8","TNNT1","CIDEC",
           "ZPLD1","PPP1R1A","SAA4","TF","FDCSP","ARL14","HOXB13","SMPX","CCDC185","KRT75")

TPM.norm <- log2(TPM.ALL + 1)
TPM.norm <- data.frame(t(TPM.norm))
TPM_RiskScore <- TPM.norm[,genes]
TPM_RiskScore<-TPM_RiskScore %>%
  mutate(RiskScore=0.24973*ANGPTL8+0.18295*IL20RB+0.24684*SAA1-0.39887*WISP2-
           0.21868*PI3+ 0.43355*CCNA1+0.44356*SYT8+0.34130*TNNT1-0.32117*CIDEC+
           0.47614*ZPLD1-0.25511*PPP1R1A-0.34612*SAA4+0.22565*TF+0.18836*FDCSP-
           0.51783*ARL14+0.25297*HOXB13+0.45957*SMPX+1.34280*CCDC185+2.10814*KRT75) %>%
  mutate(risk_group=ifelse(RiskScore>=median(RiskScore),'HRisk','LRisk'))

surv <- read.table("pancancerSurvivalData_XLu.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv <- as.data.frame(na.omit(surv[which(surv$type == "KIRC"),c("OS","OS.time")]))
surv <- surv[which(surv$OS.time > 0),] # get survival time greater than 0
colnames(surv) <- c("fustat","futime")
rownames(surv) <- paste0(rownames(surv),"-01")
surv <- surv[rownames(TPM_RiskScore),]

TPM_RiskScore <-TPM_RiskScore[,c("RiskScore","risk_group")]
survival_risk <- cbind(TPM_RiskScore,surv)

fit<- survfit(Surv(round(futime/30.5,4), fustat) ~ risk_group, data = survival_risk,type = "kaplan-meier", error = "greenwood", conf.type = "plain", 
              na.action = na.exclude)
names(fit$strata) <- gsub("Subtype=", "", names(fit$strata))
ggsurvplot(fit, data = survival_risk,
           surv.median.line = "h", # Add medians survival
           
           # Change legends: title & labels
           legend.title = "RiskGroup",
           legend.labs = c("HRisk", "LRisk"),
           # Add p-value and tervals
           pval = TRUE,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           #          tables.theme = theme_cleantable(),
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_classic() # Change ggplot2 theme
)


### Copy Number Variation ###
# compare copy number variation
seg <- read.table(file.path(comAnn.path,"KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
seg <- seg[grepl("-01A-",seg$Sample),]
seg$Sample <- substr(seg$Sample,1,15)

CS1.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol1$Subtype == "CS1"),]))
CS2.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol1$Subtype == "CS2"),]))

seg.CS1 <- seg[which(seg$Sample %in% CS1.cnv),]
seg.CS2 <- seg[which(seg$Sample %in% CS2.cnv),]

seg <- seg[which(seg$Sample %in% c(CS1.cnv,CS2.cnv)),]
seg$Sample <- as.character(seg$Sample)
rownames(seg) <- 1:nrow(seg)
write.table(as.data.frame(seg),file.path(res.path,"KIRC_all_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS1),file.path(res.path,"KIRC_CS1_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS2),file.path(res.path,"KIRC_CS2_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)

marker <- seg[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"KIRC_all_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS1[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"KIRC_CS1_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS2[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"KIRC_CS2_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# generate gistic plot
# Create a chromosomes reference objects function
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#str(chrom)
col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#013a63", alpha.f = .7)

### CS1 ###
pdf(file.path(fig.path,"CS1.cnv.scores.gistic.pdf"),14,5)

scores <- read.table("./GISTIC/CS1.scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("KIRC CS1 copy number"," ","n=",length(CS1.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, col = adjustcolor("#800f2f", alpha.f = .7),
     xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.08,ylim[2]-0.1)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#0466c8", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), bty="n", fill=c(col1,col2))
dev.off()

### CS2 ###
pdf(file.path(fig.path,"CS2.cnv.scores.gistic.pdf"),14,5)

scores <- read.table("./GISTIC/CS2.scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.4)
title=paste0("KIRC CS2 copy number"," ","n=",length(CS2.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = c(-0.3,0.7),
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, col = adjustcolor("#800f2f", alpha.f = .7),
     xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.08)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2)+2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#0466c8", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"),bty="n", fill=c(col1,col2))
dev.off()

############## Heatmap of CNA ##############
arm1 <- read.table("./GISTIC/CS1.broad_values_by_arm.txt",sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
arm2 <- read.table("./GISTIC/CS2.broad_values_by_arm.txt",sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
arm1$chr <- rownames(arm1)
arm2$chr <- rownames(arm2)
arm <- join(arm1,arm2,by="chr")
rownames(arm) <- arm$chr
arm <- arm[,-grep("chr",colnames(arm))]

library( "pheatmap" )
library("gplots")
arm[arm < 0] <- -1
arm[arm > 0] <- 1

annotation <- annCol1[,1,drop=FALSE]
annotation$sample <- rownames(annotation)

choose_matrix <- arm[,which(colnames(arm) %in% rownames(annotation))]
annotation <- annotation[which(rownames(annotation) %in% colnames(choose_matrix)),]
choose_matrix <- choose_matrix[,rownames(annotation)]

annotation_col = annotation[,1,drop=FALSE]
mycol <- colorRampPalette(c("#225ea8","white","#b5042a"))(200)
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)
library("RColorBrewer")
#brewer.pal(8,"Set2")
ann_colors = list(Subtype = c("CS1" = "#2EC4B6", "CS2" = "#E71D36"))
pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annotation_col, annotation_colors=ann_colors , show_rownames = T,show_colnames = F, 
         annotation_legend = T, legend=FALSE, 
         filename = file.path(fig.path,"heatmap_Subtype_CNA.pdf"),
         border_color = NA, fontsize_row = 12, fontsize_col = 5,cellwidth = 2.5, height=8)

############## fisher test ##############
arm.del <- arm
arm.del[arm.del>0] <- 0
arm.amp <- arm
arm.amp[arm.amp<0] <- 0

CS1.group <- intersect(colnames(arm),CS1)
CS2.group <- intersect(colnames(arm),CS2)
ans <- rep(c("CS1","CS2"),c(112,136)); names(ans) <- c(CS1.group,CS2.group)

genelist <- rownames(arm.amp[rowSums(arm.amp) > 0.05*ncol(arm.amp),])
arm.amp1 <- arm.amp[genelist,]
arm.amp1 <- as.data.frame(t(arm.amp1))
out <- matrix(0, nrow=length(genelist), ncol=3)
colnames(out) <- c(levels(factor(ans)), "pvalue")
rownames(out) <- paste(genelist, "gain", sep=" ")
for (k in 1:length(genelist)) {
  chrk <- genelist[k]
  x <- ans
  y <- arm.amp1[names(x), chrk, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["1", ]
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out1 <- as.data.frame(out)
write.table(out1, file.path(res.path,"arm.amp.fisher.txt"), quote=F, row.names = T, col.names = NA ,sep = "\t")

genelist <- rownames(arm.del[abs(rowSums(arm.del)) > 0.05*ncol(arm.del),])
arm.del1 <- arm.del[genelist,]
arm.del1 <- as.data.frame(t(arm.del1))
out <- matrix(0, nrow=length(genelist), ncol=3)
colnames(out) <- c(levels(factor(ans)), "pvalue")
rownames(out) <- paste(genelist, "del", sep=" ")
for (k in 1:length(genelist)) {
  chrk <- genelist[k]
  x <- ans
  y <- arm.del1[names(x), chrk, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["-1", ]
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out2 <- as.data.frame(out)
write.table(out2, file.path(res.path,"arm.del.fisher.txt"), quote=F, row.names = T, col.names = NA ,sep = "\t")

################ Venne ################
require(maftools)

CS1_amp <- CS1.laml.gistic@gene.summary[which(CS1.laml.gistic@gene.summary$Amp !=0 ), ]
CS1_del <- CS1.laml.gistic@gene.summary[which(CS1.laml.gistic@gene.summary$Del !=0 ), ]
CS2_amp <- CS2.laml.gistic@gene.summary[which(CS2.laml.gistic@gene.summary$Amp !=0 ), ]
CS2_del <- CS2.laml.gistic@gene.summary[which(CS2.laml.gistic@gene.summary$Del !=0 ), ]
DEG <- read.table(file.path(res.path,"consensusMOIC_TCGA-KIRC_deseq2_test_result.CS1_vs_Others.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
tmp = DEG[DEG$padj < 0.05,]
DEG_up = tmp[tmp$log2fc>1, ]
DEG_dn = tmp[tmp$log2fc< (-1), ]
gene <- DEG_up[which(rownames(DEG_up) %in% CS1_amp$Hugo_Symbol),]

library(VennDiagram)
venn <- venn.diagram(list(CS1_up = rownames(DEG_up) ,
                          CS1_dn = rownames(DEG_dn) , 
                          CS1_amp= CS1_amp$Hugo_Symbol , 
                          CS1_del= CS1_del$Hugo_Symbol ), 
                     file.path(fig.path,"CS1.venn.pdf"),  
                     col = c("cornflowerblue", "#B0DB43", "#FFC532", "#EAAFDB"),
                     fill = c("cornflowerblue", "#B0DB43", "#FFC532", "#EAAFDB"),
                     alpha = 0.6, cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), 
                     cat.cex = 1.5,rotation.degree = 0.3)

################################
library(tidyverse)
library(maftools)

CS1.laml.gistic = readGistic(gisticAllLesionsFile = './GISTIC/CS1.all_lesions.conf_95.txt',
                         gisticAmpGenesFile = './GISTIC/CS1.amp_genes.conf_95.txt', 
                         gisticDelGenesFile = './GISTIC/CS1.del_genes.conf_95.txt', 
                         gisticScoresFile = './GISTIC/CS1.scores.gistic', 
                         isTCGA = T)
CS1.laml.gistic

pdf(file.path(fig.path,"CS1.gisticChromPlot.pdf"),7,6)
gisticChromPlot(gistic = CS1.laml.gistic, markBands = "all",
                y_lims=c(-0.2,0.8))
dev.off()

pdf(file.path(fig.path,"CS1.gisticBubblePlot.pdf"),6,6)
gisticBubblePlot(gistic = CS1.laml.gistic)
dev.off()

CS2.laml.gistic = readGistic(gisticAllLesionsFile = './GISTIC/CS2.all_lesions.conf_95.txt',
                             gisticAmpGenesFile = './GISTIC/CS2.amp_genes.conf_95.txt', 
                             gisticDelGenesFile = './GISTIC/CS2.del_genes.conf_95.txt', 
                             gisticScoresFile = './GISTIC/CS2.scores.gistic', 
                             isTCGA = T)
CS2.laml.gistic

pdf(file.path(fig.path,"CS2.gisticChromPlot.pdf"),7,6)
gisticChromPlot(gistic = CS2.laml.gistic, markBands = "all", y_lims=c(-0.2,0.8))
dev.off()

pdf(file.path(fig.path,"CS2.gisticBubblePlot.pdf"),6,6)
gisticBubblePlot(gistic = CS2.laml.gistic)
dev.off()

# 挑选q<0.25的基因
sig_cytoband <- CS1.laml.gistic@cytoband.summary %>% filter(qvalues<0.25) %>% .$Unique_Name
table(grepl(pattern = "AP",sig_cytoband))
## 挑选q值小于0.25的Amplification cytoband
sig_AP_cytoband <- CS1.laml.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "AP",.)]
length(sig_AP_cytoband)
## 将Amplification 的基因挑选出来
sig_AP_gene <- CS1.laml.gistic@data %>% filter(Cytoband %in% sig_AP_cytoband) %>% .$Hugo_Symbol 
head(unique(sig_AP_gene))
## 挑选q值小于0.25的Deletion cytoband
sig_DP_cytoband <- CS1.laml.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "DP",.)]
length(sig_DP_cytoband)
## 将Deletion cytoband 的基因挑选出来
sig_DP_gene <- CS1.laml.gistic@data %>% filter(Cytoband %in% sig_DP_cytoband) %>% .$Hugo_Symbol 
head(unique(sig_DP_gene))


########### Mutation ###########
################################
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)

maf <- read_tsv(file.path(comAnn.path,"mc3.v0.2.8.PUBLIC.maf"), comment = "#")
label <- c("Tumor_Sample_Barcode","Hugo_Symbol","NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
           "Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
           "HGVSc","HGVSp","HGVSp_Short","BIOTYPE")

maf$bcr_patient_barcode <- substr(maf$Tumor_Sample_Barcode,start = 1,stop = 12)

tmp <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- tmp[which(tmp$type %in% c("KIRC")),"bcr_patient_barcode"]
KIRC.mut.maf <- maf[which(maf$bcr_patient_barcode %in% tmp),label]
KIRC.mut.maf$Tumor_Sample_Barcode <- substr(KIRC.mut.maf$Tumor_Sample_Barcode,start = 1,stop = 15)

mut.maf <- as.data.frame(KIRC.mut.maf)
mut.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% rownames(annCol1)),] 
write.table(as.data.frame(mut.maf),file.path(res.path,"KIRC_MAF.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# TMB with only SNP
mutload.KIRC <- as.data.frame(mutect.dataframe(mut.maf)); rownames(mutload.KIRC) <- mutload.KIRC$Tumor_Sample_Barcode
indel.KIRC <- as.data.frame(indel.dataframe(mut.maf)); rownames(indel.KIRC) <- indel.KIRC$Tumor_Sample_Barcode

dim(mutload.KIRC)
head(mutload.KIRC)

CS1.mut <- intersect(rownames(mutload.KIRC),CS1) #81
CS2.mut <- intersect(rownames(mutload.KIRC),CS2) #120

# tmp1 <- mutload.KIRC[CS1.mut,"TCGA_sum"]
# tmp2 <- mutload.KIRC[CS2.mut,"TCGA_sum"]
# tmp3 <- mutload.KIRC[CS3.mut,"TCGA_sum"]
# tmp4 <- mutload.KIRC[CS4.mut,"TCGA_sum"]
# p <- round(wilcox.test(tmp1,tmp2,tmp3,tmp4,alternative = "less")$p.value,3) #0.28
# df <- data.frame("TMB"=c(tmp1,tmp2),"myasthenia_gravis"=rep(c("T","F"),c(length(tmp1),length(tmp2)))); df$myasthenia_gravis <- factor(df$myasthenia_gravis,levels = c("T","F"))
# p <- ggplot(df,aes(x=myasthenia_gravis,y=TMB,fill=myasthenia_gravis)) + 
#   geom_boxplot() + scale_fill_manual(values = jco[2:1]) + 
#   ggtitle("") + annotate("text",size=3,x=1.5,y=400,label=("pValue = 0.28")) + 
#   theme(legend.position="none") + theme(axis.text.x = element_text( hjust = 1,size = 12),
#                                         plot.title = element_text(size = 12)) + xlab("") + 
#   scale_y_continuous(limits = c(0, 500))
# ggsave(file.path(fig.path,"manuscript boxplot for TMB in myasthenia_gravis TvsF.pdf"),width = 3,height = 4)
# 
# tmp1 <- indel.KIRC[CS1.mut,"TCGA_sum"]
# tmp2 <- indel.KIRC[CS2.mut,"TCGA_sum"]
# p <- round(t.test(tmp1,tmp2,na.rm=T)$p.value,3)
# 
# rm(maf); gc()

# Independent test
tmp <- as.data.frame(mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% c(CS1.mut,CS2.mut)),
                             c("Variant_Classification","Tumor_Sample_Barcode","Hugo_Symbol","Variant_Type")])
mut.binary <- matrix(0,nrow = length(unique(tmp$Hugo_Symbol)),ncol = length(c(CS1.mut,CS2.mut)),
                     dimnames = list(c(unique(tmp$Hugo_Symbol)),c(CS1.mut,CS2.mut)))

for (i in colnames(mut.binary)) {
  tmp1 <- tmp[which(tmp$Tumor_Sample_Barcode == i),]
  tmp1 <- tmp1[which(tmp1$Variant_Type == "SNP"),]
  # if(is.element("Silent",tmp$Variant_Classification)) {
  #   tmp1 <- tmp1[-which(tmp1$Variant_Classification %in% "Silent"),]
  # }
  for (j in tmp1$Hugo_Symbol)
    mut.binary[j,i] <- 1
}
mut.binary <- as.data.frame(mut.binary); #rownames(mut.binary) <- toupper(rownames(mut.binary))
write.table(mut.binary,file.path(res.path,"mut binary SNP.txt"),sep = "\t",row.names = T,col.names = NA)

source(file.path(script.path,"create.anntrack.R")) 
source(file.path(script.path,"createMutSubtype.R")) 
ans <- rep(c("CS1","CS2"),c(92,118))
names(ans) <- c(CS1.mut,CS2.mut)
genelist <- rownames(mut.binary[rowSums(mut.binary) > 0.05*ncol(mut.binary),])
binarymut <- as.data.frame(matrix(0, nrow=210, ncol=length(genelist)))
rownames(binarymut) <- names(ans)
colnames(binarymut) <- genelist
for (k in 1:length(genelist)) {
  res <- create.anntrack(samples=names(ans), subtype=createMutSubtype(mut.binary, names(ans), genelist[k]))
  binarymut[, genelist[k]] <- res$subtype
}

out <- matrix(0, nrow=length(genelist), ncol=5)
colnames(out) <- c(levels(factor(ans)), "pvalue")
rownames(out) <- paste(genelist, "Mutated", sep="_")
for (k in 1:length(genelist)) {
  genek <- genelist[k]
  x <- ans
  y <- binarymut[names(x), genek, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["Mutated", ]
  
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out <- as.data.frame(out)
#out$FDR <- p.adjust(out$pvalue)
write.table(out, file.path(res.path, "Independence test between 5% cut gene mutation and CS 1234 with SNP.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

# MutSigCV
CS1.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS1.mut),]
CS1.maf$Tumor_Seq_Allele1 <- CS1.maf$Tumor_Seq_Allele2
all(CS1.maf$Tumor_Seq_Allele1 == CS1.maf$Tumor_Seq_Allele2)
all(CS1.maf$Tumor_Seq_Allele1 == CS1.maf$Reference_Allele)
outTable <- CS1.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS1.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS2.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS2.mut),]
CS2.maf$Tumor_Seq_Allele1 <- CS2.maf$Tumor_Seq_Allele2
all(CS2.maf$Tumor_Seq_Allele1 == CS2.maf$Tumor_Seq_Allele2)
all(CS2.maf$Tumor_Seq_Allele1 == CS2.maf$Reference_Allele)
outTable <- CS2.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS2.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)),]
CS.maf$Tumor_Seq_Allele1 <- CS.maf$Tumor_Seq_Allele2
all(CS.maf$Tumor_Seq_Allele1 == CS.maf$Tumor_Seq_Allele2)
all(CS.maf$Tumor_Seq_Allele1 == CS.maf$Reference_Allele)
outTable <- CS.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS1234 mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

# mutation signature
library(deconstructSigs)
maf <- as.data.frame(mut.maf[which(mut.maf$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent") & mut.maf$Variant_Type == "SNP"),])
maf$Chromosome <- paste0("chr",maf$Chromosome)
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2")
write.table(sigs.input,file.path(res.path,"mutation.sig.input.snp.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

cut.off <- 0.05
mut.wt <- data.frame()
sigs.out.list <- list()
for (sample in rownames(sigs.input)) {
  tmp <- whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  #Plot output
  # pdf(file.path(fig.path,paste0(sample,"_plotSignatures.pdf")))
  # plotSignatures(tmp)
  # invisible(dev.off())
  # 
  # pdf(file.path(fig.path,paste0(sample,"_weightPie.pdf")))
  # makePie(tmp)
  # invisible(dev.off())
  # 
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt,tmp)
}
write.table(mut.wt,file.path(res.path,"mutation.snp.signature.weightMatrix.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

p.mut <- c()
for (i in 1:30) {
  tmp1 <- mut.wt[CS3.mut,i]
  tmp2 <- mut.wt[CS4.mut,i]
  p.mut <- c(p.mut,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.mut) <- colnames(mut.wt)[1:30]
p.mut

#####################  maftools #####################
require(maftools)

flags <- read.table(file.path(comAnn.path,"Mutation Flags 100.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- as.data.frame(mut.maf)
tmp1 <- tmp[-which(tmp$Hugo_Symbol %in% flags$FLAGS),]

CS1.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS1),]
CS2.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS2),]

clin <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
clin <- as.data.frame(na.omit(clin[which(clin$type == "KIRC"),c("OS","OS.time")]))
clin <- clin[which(clin$OS.time > 0),] # get survival time greater than 0
colnames(clin) <- c("fustat","futime")
rownames(clin) <- paste0(rownames(clin),"-01")
clin$Tumor_Sample_Barcode <- rownames(clin)
clin <- clin[rownames(annCol1),]
clin$Subtype <- annCol1$Subtype

clin1 <- clin[which(rownames(clin) %in% tmp1$Tumor_Sample_Barcode),]
clin1$Tumor_Sample_Barcode <- rownames(clin1)
laml <- read.maf(maf = tmp1, clinicalData= clin1)

library("RColorBrewer")
col = c("#1F78B4","#B2DF8A", "#33A02C" ,"#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
names(col) = c('Missense_Mutation', 'Nonsense_Mutation', "Frame_Shift_Del" ,'Splice_Site','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation')
titvCol = RColorBrewer::brewer.pal(n = 6, name = 'Set2')
names(titvCol) = c('T>G','T>A','T>C','C>T','C>G','C>A')
pdf(file.path(fig.path,"KIRC_mutation_plotmafSummary.pdf"),width=8, height=6)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE , color = col, titvColor = titvCol, fs=0.9,top=8)
invisible(dev.off())
#We will draw oncoplots for top ten mutated genes.
#oncoplot(maf = laml, top = 10, fontSize = 12)

#Changing colors for variant classifications
#brewer.pal(10,"Paired")
col =  RColorBrewer::brewer.pal(n = 8, name = 'Paired')
col =  c("#A6CEE3","#1599ba","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
names(col) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
#Color coding for prediction,myasthenia_gravis
cincolors = list(Subtype = c(CS1 = "#2EC4B6", CS2 = "#E71D36", 
              CS3 = "#FF9F1C", CS4 = "#BDD5EA"))
pdf(file.path(fig.path,"KIRC_mutation_oncoplot.pdf"),width=10, height=8)
oncoplot(maf = laml,top = 10, removeNonMutated=F,clinicalFeatures = c("Subtype"), 
         colors = col, sortByAnnotation = TRUE, groupAnnotationBySize=FALSE,annotationColor = cincolors,
         legendFontSize = 1.5 ,titleFontSize = 1.5)
invisible(dev.off())

pdf(file.path(fig.path,"KIRC_mutation_exclusive and Co-occurrence.pdf"),width=6.5, height=6)
somaticInteractions(maf = laml, top = 20, fontSize=0.7,colPal="PiYG",
                    showSum=FALSE, pvalue = c(0.05, 0.1))
invisible(dev.off())

################### CS2 VS CS3 ###################
clinCS1 <- clin[which(clin$Subtype == "CS1"),]
clinCS2 <- clin[which(clin$Subtype == "CS2"),]

mafCS1 = read.maf(maf = CS1.mut.maf, clinicalData= clinCS1 )
mafCS2 = read.maf(maf = CS2.mut.maf, clinicalData= clinCS2 )
pt.vs.rt <- mafCompare(m1 = mafCS1, m2 = mafCS2, m1Name = 'CS1', m2Name = 'CS2', minMut = 5)
pdf(file.path(fig.path,"forestPlot for CS1 vs CS2 mutation.pdf"),width=8, height=8)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'))
invisible(dev.off())

col =  c("#1F78B4","#B2DF8A", "#FDBF6F","#33A02C")
names(col) = c('Missense_Mutation', 'Nonsense_Mutation', "In_Frame_Ins" ,'Frame_Shift_Del')
pdf(file.path(fig.path,"KIRC_mutation_oncoplot for CS1 vs CS2.pdf"),width=6, height=5)
coBarplot(m1 = mafCS1, m2 = mafCS2, m1Name = 'CS1', m2Name = 'CS2')
invisible(dev.off())

################### ITH ################### 
library("mclust")
maf <- read_tsv(file.path(data.path,"mc3.v0.2.8.PUBLIC.maf"), comment = "#")
label <- c("Tumor_Sample_Barcode","Hugo_Symbol","NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
           "Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
           "HGVSc","HGVSp","HGVSp_Short","BIOTYPE","t_ref_count","t_alt_count","n_alt_count","n_ref_count")
maf$bcr_patient_barcode <- substr(maf$Tumor_Sample_Barcode,start = 1,stop = 12)
tmp <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- tmp[which(tmp$type %in% c("KIRC")),"bcr_patient_barcode"]
KIRC.mut.maf <- maf[which(maf$bcr_patient_barcode %in% tmp),label]
KIRC.mut.maf$Tumor_Sample_Barcode <- substr(KIRC.mut.maf$Tumor_Sample_Barcode,start = 1,stop = 15)
mut.maf <- as.data.frame(KIRC.mut.maf)
mut.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% rownames(annCol1)),]

mut.maf$vaf <- (mut.maf$t_alt_count+mut.maf$n_alt_count)/(mut.maf$t_alt_count+mut.maf$n_alt_count+mut.maf$t_ref_count+mut.maf$n_ref_count)
flags <- read.table(file.path(comAnn.path,"Mutation Flags 100.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- as.data.frame(mut.maf)
#tmp1 <- tmp[-which(tmp$Hugo_Symbol %in% flags$FLAGS),]
tmp1 <- tmp
laml <- read.maf(maf = tmp1, clinicalData= clin )

ITH = inferHeterogeneity(maf = laml, tsb = rownames(annCol1) , vafCol = 'vaf')
ITH1 <- unique(ITH$clusterData[,c("Tumor_Sample_Barcode","MATH","MedianAbsoluteDeviation")])
ITH.CS1 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS1),]
ITH.CS1$Subtype <- "CS1"
ITH.CS2 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS2),]
ITH.CS2$Subtype <- "CS2"
ITH2 <- rbind(ITH.CS1,ITH.CS2)

# shapiro.test(ITH.CS1$MATH) #p=0.7093
# shapiro.test(ITH.CS2$MATH) #p=0.0824
# shapiro.test(ITH.CS3$MATH) #p=0.7093
# shapiro.test(ITH.CS4$MATH) #p=0.0824
# 
# bartlett.test(MATH~MG,data=ITH2) #p=0.8798

kruskal.test(MATH~Subtype, ITH2)
t.test(MATH~Subtype, ITH2)

library(doBy)  
library(ggplot2)
dat <- summaryBy(MATH~Subtype, ITH2, FUN = c(mean, sd))
p <- ggplot(dat, aes(Subtype, MATH.mean, fill = Subtype)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = MATH.mean - MATH.sd, ymax = MATH.mean + MATH.sd), width = 0.15, size = 0.5) +
  scale_fill_manual(values=c("CS2" = "#E71D36", 
                             "CS3" = "#FF9F1C")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Cluster', y = 'MATH', title = 't-test: p-value = 0.058')
ggsave(file.path(fig.path,"plot for MATH between 2Subtype.pdf"),width = 5,height = 4)

###################
### methylation ###
library(missMethyl)
library(ChAMP)
library(data.table)

tmp <- fread(file.path(data.path,"TCGA-KIRC.methylation450.tsv"),sep = "\t",check.names = F,stringsAsFactors = F)
tmp <- as.data.frame(tmp); rownames(tmp) <- tmp[,1]; tmp <- tmp[,-1]
colnames(tmp) <- substr(colnames(tmp),1,15)

meth.tumor <- colnames(tmp)[which(substr(colnames(tmp),14,14)=="0")]
meth.normal <- colnames(tmp)[which(substr(colnames(tmp),14,14)=="1")]

CS1.samples.meth <- intersect(meth.tumor,CS1)
CS2.samples.meth <- intersect(meth.tumor,CS2)

write.table(tmp[,c(CS1.samples.meth,CS2.samples.meth)],file.path(data.path,"KIRC.CS1-CS2.Meth450.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
pd.meth <- data.frame(Sample_Name=c(CS1.samples.meth,CS2.samples.meth),
                      Sample_Group=rep(c("CS1","CS2"),c(length(CS1.samples.meth),length(CS2.samples.meth))),
                      row.names = c(CS1.samples.meth,CS2.samples.meth),
                      stringsAsFactors = F)
write.table(pd.meth,file.path(data.path,"pd.meth.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
rm(tmp); gc()

library(ChAMP)
myLoad <- ChAMP:::champ.read(betaFile = file.path(data.path,"KIRC.CS1-CS2.Meth450.txt"),
                             sampleSheet = file.path(data.path,"pd.meth.txt"), 
                             resultsDir = file.path(res.path,"champ_resultsDir"))

myLoad$pd$Sample_Group <- factor(myLoad$pd$Sample_Group,levels = c("CS1","CS2"))
# myLoad$beta <- as.matrix(na.omit(myLoad$beta))
dim(myLoad$beta)
rownames(myLoad$beta)<-myLoad$beta$X
myLoad$beta <- myLoad$beta[,-1]
myLoad$beta<-as.matrix(myLoad$beta)
colnames(myLoad$beta) <- myLoad$pd$Sample_Name

myFilter <- champ.filter(beta = myLoad$beta,
                         pd = myLoad$pd,
                         arraytype = "450K")
myImpute <- champ.impute(beta=myFilter$beta,pd=myFilter$pd,SampleCutoff = 1)
myNorm <- champ.norm(beta = myImpute$beta)
myDMP <- champ.DMP(beta = myNorm,
                   pheno = myImpute$pd$Sample_Group)
write.table(myDMP[[1]],file.path(res.path,"champ.myDMP.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
#DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)

# rm(myLoad); gc()
# rm(myFilter); gc()
myDMR <- champ.DMR(beta = myNorm,
                   pheno = myImpute$pd$Sample_Group,
                   method="Bumphunter")
write.table(myDMR,file.path(res.path,"champ.myDMR.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

comRFun.path <- file.path(tumor.path,"commonFun") 
source(file.path(comRFun.path,"champ.GSEA.R")) 
myGSEA <- champ.GSEA(beta = myNorm,
                     DMP = myDMP[[1]],
                     DMR = myDMR,
                     method = "fisher")
#myGSEA <- ChAMP:::champ.GSEA(beta=myFilter$beta,DMP=myDMP[[1]],
#                     DMR=myDMR,
#                     arraytype="450K",
#                     adjPval=0.1,
#                     method="fisher")

Glist <- list()
Glist[["hyperC1"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta > 0)]
Glist[["hyperC2"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta < 0)]
Glist[["hyperC1_promCGI"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta > 0 & myDMP[[1]]$cgi == "island" & myDMP[[1]]$feature %in% c("TSS1500","TSS200"))]
Glist[["hyperC2_promCGI"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta < 0 & myDMP[[1]]$cgi == "island" & myDMP[[1]]$feature %in% c("TSS1500","TSS200"))]

# Glist2 <- list()
# Glist2[["hyperC1_promCGI"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta > 0 & myDMP[[1]]$cgi == "island" & myDMP[[1]]$feature %in% c("TSS1500","TSS200"))]
# Glist2[["hyperC2_promCGI"]] <- rownames(myDMP[[1]])[which(myDMP[[1]]$deltaBeta < 0 & myDMP[[1]]$cgi == "island" & myDMP[[1]]$feature %in% c("TSS1500","TSS200"))]

myGSEA.add <- champ.GSEA(beta=myNorm,
                         DMP=myDMP[[1]],DMR=myDMR,CpGlist = Glist,
                         pheno = myImpute$pd$Sample_Group,
                         method = "fisher",
                         arraytype = "450K")

write.table(myGSEA.add[[1]],file.path(res.path,"champ.myDMP.gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(myGSEA.add[[2]],file.path(res.path,"champ.myDMR.gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)

# myGSEA.add2 <- champ.GSEA(beta=myNorm,
#                           DMP=NULL,DMR=NULL,CpGlist =  Glist2,
#                           pheno = myImpute$pd$Sample_Group,
#                           method = "gometh",
#                           arraytype = "450K")

for (i in 1:length(myGSEA.add)) {
  label <- names(myGSEA.add)[i]
  if(label %in% c("DMP","DMR")) {
    next()}
  tmp <- myGSEA.add[[i]]
  write.table(tmp,file.path(res.path,paste0("ChAMP_TvsF_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}

write.table(probe.features,file.path(comAnn.path,"Meth450Anno.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

meth450 <- read.table(file.path(data.path,"KIRC.CS1-CS2.Meth450.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
meth <- log(meth450/(1-meth450))
DMP <- read.table(file.path(res.path,"champ.myDMP.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill=TRUE)
group <- read.table(file.path(data.path,"pd.meth.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
#DMP_up = DMP[DMP$myasthenia_gravis_T_AVG>0.3 & DMP$myasthenia_gravis_F_AVG<0.3 & DMP$deltaBeta<(-0.2), ]
#DMP_down = DMP[DMP$myasthenia_gravis_T_AVG<0.3 & DMP$myasthenia_gravis_F_AVG>0.3 & DMP$deltaBeta>0.2, ]
#DMP1 <- DMP[c( rownames( DMP_up ), rownames( DMP_down )), ]
DMP <- myDMP[[1]]
DMP_up = DMP[DMP$CS1_AVG>0.3 & DMP$CS2_AVG<0.3 & DMP$deltaBeta<(-0.2), ]
DMP_down = DMP[DMP$CS1_AVG<0.3 & DMP$CS2_AVG>0.3 & DMP$deltaBeta>0.2, ]
DMP1 <- DMP[c( rownames( DMP_up ), rownames( DMP_down )), ]
DEG <- read.table(file.path(res.path,"mRNA_deseq2_test_result.CS3_vs_CS4.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
DEG = DEG[DEG$padj < 0.25,]
DEG_Z = DEG[DEG$log2FoldChange>1.5, ]
DEG_F = DEG[DEG$log2FoldChange< (-1.5), ]
DEG=DEG[c( rownames( DEG_Z ), rownames( DEG_F )),]
DMP1 <- DMP1[which(DMP1$gene %in% rownames(DEG)),]
DEG1 <- DEG[which(rownames(DEG) %in% DMP1$gene),]

library( "pheatmap" )
library("gplots")
choose_matrix = meth[c( rownames( DMP_up ), rownames( DMP_down )), ]
#choose_matrix <- sweep(choose_matrix,1,apply(choose_matrix,1,median))
choose_matrix = standarize.fun(indata=choose_matrix,halfwidth = 3)

ann_colors = list( myasthenia_gravis = c(myasthenia_gravis_T = cyan, myasthenia_gravis_F = peach),
                   EMT = c(Epi = blue, Mes = yellow))
rownames( annotation_col ) = colnames( choose_matrix )
mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)
pheatmap(choose_matrix, color = mycol,cluster_col=F, cluster_rows = F,
         gaps_col = 34 , gaps_row = 174,
         annotation_col = annotation_col, annotation_colors=ann_colors , show_rownames = F,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_THYM_Meth450.pdf"), 
         width = 6, height = 5,border_color = NA, fontsize = 6)

####################### Drug ###########################
library(oncoPredict)
library(pRRophetic)

library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

set.seed(12345)
trainingExprData = readRDS(file =paste0 (data.path,"/GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
dim(trainingExprData) #51847 829
trainingPtype = readRDS(file = paste0 (data.path,"/GDSC2_Res.rds"))
dim(trainingPtype) #829 545 
trainingPtype <- trainingPtype[,c("Ruxolitinib_1507","Cisplatin_1005")]

testExpr<- read.table(file.path(res.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
testExpr<- as.matrix(testExpr)
dim(testExpr)  

calcPhenotype(trainingExprData = trainingExprData,
              trainingPtype = trainingPtype,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'homogenizeData' )


prefix      = "BOXVIOLIN OF ESTIMATED IC50 (5methods)"
DrugPredictions <- read.csv(file.path("~/Ruan/KIRC/movics_pipeline/calcPhenotype_Output/DrugPredictions.csv"),check.names = F, header = T,stringsAsFactors = F,row.names = 1)
DrugPredictions <- DrugPredictions[rownames(annCol_new),]
Drug <- cbind(annCol_new,DrugPredictions)
statistic = "kruskal.test"

ic50.test  <- kruskal.test(Drug$Ruxolitinib_1507 ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$Ruxolitinib_1507,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("Ruxolitinib_1507",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- c("#2EC4B6", "#E71D36")
names(colvec) <- c("CS1","CS2")
n.moic=2
p <- ggplot(data = Drug,
            aes(x = Subtype, y = Ruxolitinib_1507, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"Ruxolitinib_1507")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = min(Drug$Ruxolitinib_1507))
p
outFig <- paste0(prefix, " for Ruxolitinib_1507", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)


ic50.test  <- kruskal.test(Drug$ML210 ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$ML210,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("ML210",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
p <- ggplot(data = Drug,
            aes(x = Subtype, y = ML210, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"ML210")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = 1.2)
p
outFig <- paste0(prefix, " for ML210", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)


ic50.test  <- kruskal.test(Drug$erastin ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$erastin,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("erastin",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
p <- ggplot(data = Drug,
            aes(x = Subtype, y = erastin, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"erastin")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = 6)
p
outFig <- paste0(prefix, " for erastin", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)


######################## Regulon ########################
# load R package
library(RTN)

# load TFs
tfs <- read.table(file.path(data.path,"1639TFs.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
TPM <- read.table(file.path(data.path,"TCGA_KIRC_mRNA_TPM.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

# get intersection of TFs
regulatoryElements <- intersect(tfs$GENE, rownames(TPM))

# Run the TNI constructor
rtni_KIRC <- tni.constructor(expData = as.matrix(TPM), 
                             regulatoryElements = regulatoryElements)

# Compute the reference regulatory network by permutation and bootstrap analyses.
# Please set 'spec' according to your available hardware
options(cluster=snow::makeCluster(spec=12, "SOCK")) # 12 cores
rtni_KIRC <- tni.permutation(rtni_KIRC, pValueCutoff = 0.05, nPermutations = 1000)# 1e-5
rtni_KIRC <- tni.bootstrap(rtni_KIRC, nBootstraps = 1000)
snow::stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network
rtni_KIRC <- tni.dpi.filter(rtni_KIRC, eps = 0, sizeThreshold = TRUE, minRegulonSize = 5)# 15
tni.regulon.summary(rtni_KIRC)

# Save the TNI object for subsequent analyses
saveRDS(rtni_KIRC, file=file.path("rtni_KIRC.rds"))

# Compute regulon activity for individual samples
rtnigsea_KIRC <- tni.gsea2(rtni_KIRC, regulatoryElements = regulatoryElements)
KIRC_regact <- tni.get(rtnigsea_KIRC, what = "regulonActivity")
saveRDS(KIRC_regact,file = file.path( "KIRC_regact.rds")) # 



# Get sample attributes 
matrix <- t(KIRC_regact$dif)
matrix <- matrix[,rownames(annCol1)]
try <- matrix[c("NFE2L1","NFE2L2"),]

# Plot regulon activity profiles
pheatmap(try, 
         main="TCGA-KIRC cohort subset (n=437 samples)",
         annotation_col = annCol1, annotation_colors=annColors1, cluster_col = F,
         show_colnames = FALSE, annotation_legend = FALSE, 
         clustering_method = "ward.D2", fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")


raw_count <- read.table(file.path(data.path,"TCGA_KIRC_Count.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
CS1.raw_count <- raw_count[rownames(TPM.HUGO),CS1]
CS2.raw_count <- raw_count[rownames(TPM.HUGO),CS2]
write.table(CS1.raw_count,file.path(data.path,"CS1.raw_count.txt"),sep = "\t",row.names = T)
write.table(CS2.raw_count,file.path(data.path,"CS2.raw_count.txt"),sep = "\t",row.names = T)
