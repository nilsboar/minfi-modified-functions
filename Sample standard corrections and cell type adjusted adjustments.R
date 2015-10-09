library(lumi)
library(minfi)
library(FlowSorted.Blood.450k)
library(GenomicRanges)
library("IlluminaHumanMethylation450k.db")
library("FDb.InfiniumMethylation.hg19")
library(FDb.UCSC.snp137common.hg19)

#################################################################
# Read the targets file and data
#################################################################
# Set the path to the directory with the files 
baseDir <-"~/Documents/big data/methylation_data"
resultDir<-"~/Documents/analysis/AIDS/methylation study/results_directory"
setwd(baseDir)
# Read the targets files from the *.csv file
targets <- read.450k.sheet(baseDir, pattern = )
head(targets); dim(targets)
sub(baseDir, "", targets$Basename) #list the files
# Read in the data
RGset.example <- read.450k.exp(base = baseDir, targets = targets)
rm(targets)

# # # # save.image(file = "~/Documents/big data/methylation big data/RG data read in.Rdat")
# # # # load(file = "~/Documents/big data/methylation big data/RG data read in.Rdat")

###################################################################
# Get reference blood data
###################################################################
referencePkg <- "FlowSorted.Blood.450k"
data(list = referencePkg)
referenceRGset <- get(referencePkg)

###################################################################
# Find bad probes, in both example and reference blood data
###################################################################
detP<-detectionP(RGset.example) # this gets the p-values
failed.01<- detP > 0.01  #set detection cutoff
### LIST Bad Probes
failedProbes.example <-rownames(failed.01)[rowMeans(failed.01)>0.05]     #list of probes that failed in more than 5% of the sample
sum(rowMeans(failed.01)>0.05) # how many probes failed in more than 5% of samples?  552.
detP<-detectionP(referenceRGset) # this gets the p-values
failed.01<-detP > 0.01  #set detection cutoff
### LIST Bad Probes
failedProbes.bloodref <-rownames(failed.01)[rowMeans(failed.01)>0.05]     #list of probes that failed in more than 5% of the sample
sum(rowMeans(failed.01)>0.05) # how many probes failed in more than 5% of samples?  552.
# probes that failed in either example or blood reference
failedProbes <- unique(c(failedProbes.example, failedProbes.bloodref))

###################################################################
# Combine the example and reference blood data;  Make pData consistent to allow this
###################################################################
pdata.celldat.full <- pData(referenceRGset)
pdata.exampledat.full <- pData(RGset.example)
# look at the pdata
names(pdata.celldat.full)
names(pdata.exampledat.full)
head(pdata.celldat.full)
head(pdata.exampledat.full)
# make Pdata consistent so sets can be joined
Pdata.example <- pdata.exampledat.full[, c("Individual.ID", "Slide", "Array", "Batch", "HIV.infection", "Stage")]
Pdata.celldat <- pdata.celldat.full[, c("Sample_Name", "Slide", "Array", "CellType")]
Pdata.celldat $Slide <- as.character(pdata.celldat.full $Slide)
Pdata.example$CellType <- "example.pbmc"
Pdata.example$Study <- "example.methy"
Pdata.celldat$Batch <- as.integer(4) # our study batches are 1-3
Pdata.celldat$Study <- "bloodcell.methy"
Pdata.celldat$HIV.infection <- "NEG"
Pdata.celldat$Stage <- "PRE"
Pdata.celldat <- Pdata.celldat[,c("Sample_Name", "Slide", "Array",  "Batch", "HIV.infection",  "Stage", "CellType",  "Study")]
names(Pdata.example) <- c("Sample_Name", "Slide", "Array",  "Batch", "HIV.infection",  "Stage", "CellType",  "Study")
rownames(Pdata.celldat) <- paste(Pdata.celldat $Slide, Pdata.celldat $Array, sep = "_")
pData(referenceRGset) <- Pdata.celldat
pData(RGset.example) <- Pdata.example
# pdata.celldat.test <- pData(referenceRGset)
# pdata.exampledat.test <- pData(RGset.example)
# summary(pData(referenceRGset))
# summary(pData(RGset.example))

###################################################################
# Combine the example and reference blood data, preprocess, drop bad and irrelevant loci
###################################################################
combinedRGset.all <- combine(RGset.example, referenceRGset)
# # # # # save.image(file = "~/Documents/big data/methylation big data/working results cmbnd sets")
# # # # # load(file = "~/Documents/big data/methylation big data/working results cmbnd sets")
Mset.raw <- preprocessRaw(combinedRGset.all)
Mset.xmeth <- dropMethylationLoci(Mset.raw, dropRS=TRUE, dropCH=TRUE)
Mset <- Mset.xmeth[!rownames(Mset.xmeth)%in% failedProbes]
rm(Mset.raw, Mset.xmeth, referenceRGset, combinedRGset, combinedRGset.all, detP, failedProbes, failedProbes.example, failedProbes.bloodref, failed.01)
rm(failedProbes, failedProbes.example, failedProbes.bloodref, failed.01)
rm(FlowSorted.Blood.450k)
#getMeth and getUnmeth return methylated and unmethylated intensities 
me <- getMeth(Mset)  
unme <- getUnmeth(Mset)
sum(colnames(me)!=colnames(unme))
sum(rownames(me)!=rownames(unme))
# view data
me[1:20,1:8]
unme[1:20,1:8]
B=getBeta(Mset, type = "Illumina")
M <- log2((me+1)/(unme+1))
# view data
B[1:20,1:8]
M[1:20,1:8]

###################################################################
#  Need shorter and more informative names for the samples, to get useful diagonistic graphs
###################################################################
pheno = pData(Mset)
pheno$ordervar <- 1:nrow(pheno)  # to maintain the order of the pheno data, when we substitute a new pheno file
pheno$slide_array <- row.names(pheno)
mypheno <- read.csv("~/Documents/analysis/AIDS/methylation study/methylation data/combPdatToLoad.csv", row.names = 1)
comb.pheno <- merge(pheno[,c("slide_array", "ordervar")], mypheno, by = "slide_array")
comb.pheno <- comb.pheno[order(comb.pheno$ordervar),]
# test sort
# comb.pheno$slide_array  
# pheno$slide_array
newpheno <- comb.pheno[,names(mypheno)]
# order again to be safe
row.names(newpheno) <- newpheno$slide_array # these identify samples by category (disease state or cell type)
newpheno <- newpheno[colnames(me), ]
# change row.names to useful sample names
row.names(newpheno) <- newpheno$SampleName  # these identify samples by category (disease state or cell type)
# check names again
newpheno$slide_array
colnames(me)
rownames(newpheno)
# reset the pheno data in the various methylation data sets:
# merge above must assure that the names are consistent
colnames(me) <- rownames(newpheno)
colnames(unme) <- rownames(newpheno)
colnames(B) <- rownames(newpheno)
colnames(M) <- rownames(newpheno)

###################################################################
#  bring in annotation data, assemble Mset
###################################################################
anno=read.csv("~/Documents/analysis/AIDS/methylation study/methylation data/GPL13534_HumanMethylation450_15017482_v.1.1.csv",header=T, sep=',', stringsAsFactors =F)
color=anno[,c('Genome_Build','CHR','MAPINFO','Strand','Probe_SNPs','Color_Channel','UCSC_RefGene_Name','UCSC_RefGene_Accession','UCSC_RefGene_Group','UCSC_CpG_Islands_Name','Relation_to_UCSC_CpG_Island')]
colnames(color)=toupper(colnames(color))
rownames(color)=anno$IlmnID 
color=color[rownames(M),]
pd <- new("AnnotatedDataFrame", data = as.data.frame(newpheno)) # replacing the original sample names (column names in the data format) with short meaningful names
# pd <- new("AnnotatedDataFrame", data = as.data.frame(pheno))
fd <- new("AnnotatedDataFrame", data = as.data.frame(color))
AD=new.env()
assign("methylated", me, envir = AD)
assign("unmethylated", unme, envir = AD)
assign("exprs", M, envir = AD)
assign("Betas", B, envir = AD)
## create a new Mset with combined 
## integrate phenoData, featureData and expressionData
Mset.combined <-new("MethyLumiM", assayData=AD, phenoData=pd,featureData=fd)
# # # # save(Mset.combined,file='~/Documents/big data/methylation big data/Msetcombined.tstpipe.rdata')
###################################################################
# next after 5/11 continue from here
###################################################################
# # # # load('~/Documents/big data/methylation big data/Msetcombined.tstpipe.rdata')


Mset <- Mset.combined
###################################################################
# remove chrX and ChrY 
# this is standard, but should be reconsidered here, since all subjects are males 
###################################################################
Yidx=featureData(Mset)$CHR=='Y'
Xidx=featureData(Mset)$CHR=='X'
Mset2=Mset[ !(Yidx|Xidx),]
################ color balance adjustment #################
## ## Color balance adjustment between two color channels
lumiMethy.c.adj <- lumiMethyC(Mset2)
## background correction
lumiMethy.bc.adj <- lumiMethyB(lumiMethy.c.adj, method="bgAdjust2C")
#### quantile normalization
lumiMethy.bc.q<- lumiMethyN(lumiMethy.bc.adj, method='quantile')
 
# save different step outputs to plot changes 
Mset.before.adj <- Mset2
save(Mset.before.adj, lumiMethy.c.adj, lumiMethy.bc.adj, lumiMethy.bc.q,file='~/Documents/big data/methylation big data/methyset adj stages pipetest.rdata')
save(lumiMethy.bc.q,file='~/Documents/big data/methylation big data/Mset adj qnorm newnames pipetest.rdata')
load('~/Documents/big data/methylation big data/Mset adj qnorm newnames pipetest.rdata')
# from density plot, PRE7 and PRE18 appear anomalous, remove them
colnames(lumiMethy.bc.q)
badsubs <- which(colnames(lumiMethy.bc.q)  %in% c("PRE7", "PRE18"))
lumiMethy.bc.q <- lumiMethy.bc.q[, -badsubs]
mdataBetas <- estimateBeta(lumiMethy.bc.q, returnType = "matrix")
mdataPdat <- pData(lumiMethy.bc.q)

###################################################################
#		Batch adjustment with ComBat--doesn't seem to help but could be explored further
###################################################################
library(sva)
modcombat <- model.matrix(~1, data = mdataPdat)
Batch <- mdataPdat$Batch
combat.mdata <- ComBat(dat = mdataBetas, batch = Batch, mod = modcombat, prior.plots = T)
save(combat.mdata,file='~/Documents/big data/methylation big data/betas adjusted combat.rdata')

###################################################################
# Adjustment of methylation profiles for cell composition
# This is done by modifying two minfi functions
###################################################################

insertSource('~/Documents/analysis/R/R code and data/AIDS R/AIDS R code/methylation R code/minfi and lumi pls modified fcns/my cell type fcns.R', package = "minfi")
# untrace(c(estimateCellCounts, pickCompProbes, projectCellType))  # how to restore package with untrace??cct.output <- estimateCellCounts(lumiMethy.bc.q)

cct.output <- estimateCellCounts(lumiMethy.bc.q)
cell.counts <- cct.output[[1]]
cell.cts.tst <- cct.output[[2]]
ref.table <- cct.output[[3]]
cellAdjData <- cct.output[[4]]
cd4_max <- cct.output[[5]]
cd4.min <- cct.output[[6]]
rm(cct.output)


# print histograms of estimated cell composition in sample, and from adjusted data
# in adjusted data all samples should appear to have the same cell compostion; testing this
par(mfrow=c(6,2))
hist(cell.counts[,"CD4T"], breaks = 0.02*(0:40), col = "blue", main = "CD4 T cells", xlab = "")
hist(cell.cts.tst[,"CD4T"], breaks = 0.02*(0:40), col = "blue", main = "CD4 T cells", xlab = "")
hist(cell.counts[,"CD8T"], breaks = 0.02*(0:40), col = "green", main = "CD8 T cells", xlab = "")
hist(cell.cts.tst[,"CD8T"], breaks = 0.02*(0:40), col = "green", main = "CD8 T cells", xlab = "")
hist(cell.counts[,"Bcell"], breaks = 0.02*(0:40), col = "orange", main = "B cells", xlab = "")
hist(cell.cts.tst[,"Bcell"], breaks = 0.02*(0:40), col = "orange", main = "B cells", xlab = "")
hist(cell.counts[,"Gran"], breaks = 0.02*(0:40), col = "yellow", main = "Granulocytes", xlab = "")
hist(cell.cts.tst[,"Gran"], breaks = 0.02*(0:40), col = "yellow", main = "Granulocytes", xlab = "")
hist(cell.counts[,"Mono"], breaks = 0.02*(0:40), col = "purple", main = "Monocytes", xlab = "")
hist(cell.cts.tst[,"Mono"], breaks = 0.02*(0:40), col = "purple", main = "Monocytes", xlab = "")
# hist(cell.counts[,"NK"], breaks = 0.02*(0:40), col = "ochre")
# hist(cell.cts.tst[,"NK"], breaks = 0.02*(0:40), col = "ochre")
hist(cell.counts[,"Eos"], breaks = 0.02*(0:40), col = "turquoise", main = "Eosinophils", xlab = "Frequency distribution over samples, input data")
hist(cell.cts.tst[,"Eos"], breaks = 0.02*(0:40), col = "turquoise", main = "Eosinophils", xlab = "Frequency distribution over samples, methylation profile subtracted")
par(mfrow=c(1,1))

save(cellAdjData, file = "~/Documents/big data/methylation big data/cell adj data example pipeline.Rdat")


