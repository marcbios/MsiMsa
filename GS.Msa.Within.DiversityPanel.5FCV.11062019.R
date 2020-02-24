
rm(list=ls())
setwd("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/Analysis.4.Revisions/")
date="11112019"
#setwd("")
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")
#library(heritability)
library(rrBLUP)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/CDMean.Algorithm.R")
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Within.Diversity.Panel.GS.RRBLUP.Algorithm.CDmean.PCAtrn.Revised.11052019.R")
source("https://raw.githubusercontent.com/lvclark/R_genetics_conv/master/pairwiseJostDnumeric.R")

#data_hmp <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.Geno.5pmaf.90missing.42288snps.594lines.hmp.txt", head=F)
#myG <- data_hmp
#myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
#dat=myGAPIT$GD

#nan=dat[,1]
#nan <- as.vector(as.matrix(nan))
#dat <- dat[,2:ncol(dat)]
#dat=as.matrix(dat)
#dat[1:6,1:6]

#rownames(dat)=nan
#dat[1:6,1:6]


geno.all <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Geno.All.Combined.Msi.Msa.09F2.1373lines.5141snps.11052019.csv", header=T)
geno.all[1:6,1:6]
name.geno <- as.vector(geno.all$X)
df.nam.geno <- data.frame(name.geno)
names(df.nam.geno) <- "Taxa"
head(df.nam.geno)
df.nam.geno$Taxa <- as.character(df.nam.geno$Taxa)


#write.table(com.taxa, "com.taxa.numeric.pheno.geno.Msi", species,"csv", sep=",", quote=F, row.names = F, col.names = T)
phen.taxa <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.Raw.Selection.Index.Phenotypes.594lines.11052019.csv", header=T)
phen.taxa$Taxa <- as.character(phen.taxa$Taxa)
head(phen.taxa)
phen.taxa <- phen.taxa[,-c(6,7,11)] # Remove traits not needed and pop column
head(phen.taxa)

com.taxa <- merge(phen.taxa, df.nam.geno)
dim(com.taxa)
head(com.taxa)
com.taxa$Taxa <- as.character(com.taxa$Taxa)

com.taxa.df <- data.frame(com.taxa$Taxa)
colnames(com.taxa.df) <- "Taxa"

phen<-com.taxa


Msa.polyRAD.geno <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/geno.msa.polyRAD.txt", header = T)
rownames(Msa.polyRAD.geno) <- Msa.polyRAD.geno$X
Msa.polyRAD.geno <- Msa.polyRAD.geno[,-1]
Msa.polyRAD.geno <- as.matrix(Msa.polyRAD.geno)
Msa.polyRAD.geno[1:6,1:6]
Msa.polyRAD.geno.2 <- Msa.polyRAD.geno*2
Msa.polyRAD.geno.2 <- round(Msa.polyRAD.geno.2)
Msa.polyRAD.geno.2[1:6,1:6]
Msa.polyRAD.geno.3 <- Msa.polyRAD.geno.2-1
A.Aher <- A.mat(Msa.polyRAD.geno.3)
#A.Aher <- A.mat(dat2, impute.method="mean",return.imputed = TRUE) # additive relationship matrix

#A.Aher <- A.mat(dat2, impute.method="mean",return.imputed = TRUE) # additive relationship matrix

#h2<-NULL
#for(i in 2:dim(phen)[2]){
#  y<-phen[,i]
#  y[is.na(y)]<-mean(y,na.rm=T)
#  model <- mixed.solve(y,K=A.Aher)
#  h_sq<-model$Vu/(model$Vu+model$Ve)
#  h2<-c(h2,h_sq)
  
#}

#h2<-data.frame(h2,traits=c(names(phen)[-1]))

#write.table(h2, "/Users/omo/Box Sync/Desktop/PostDoc/Manuscript2/Data/Msa.trait.heritablility.polRAD.11112019.csv", sep=",", quote=F, row.names = F, col.names = T)

her.vec <- read.csv('http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/DiversityPanel.GS/Msa/Msa.trait.heritablility.polRAD.11112019.csv', header=T)
her.vec$trait <- as.character(her.vec$traits)

com.geno <- geno.all[match(com.taxa$Taxa, geno.all$X),]
com.geno[1:6,1:6]
dim(com.geno)
com.geno2 <- com.geno[,-1]
rownames(com.geno2) <- as.character(com.geno$X)


polyRAD.Geno <- com.geno2

head(com.taxa)
phenames <- as.character(colnames(com.taxa[,-1]))

her.vec$trait <- as.character(her.vec$traits)

PCA.Msa <- read.table("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.polyRAD.estimated.PCA.txt", header=T)
PCA.Msa[1:6,1:6]

Opt.PCA=5
num.vec <- 200 #,150,200,250,300,350
s <- 0.2  #  enter selection index here (proportion to select)                 
Iter=2500
burn=1000
species <- "Msa" #either Msi or Msa
# m_train1
k = 5
itrn=3000
cycles=10

trn.app <- c( "All", "CDmean", "Random")

number.of.folds=k


Groups.3.Pop3 <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Groups.3.Pop3.txt", header=T)
geno.diplo.tetra <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/geno.diplo.tetra.txt", header=T)


Within.Diversity.Panel.GS.Pipeline(phenames=phenames, polyRAD.Geno=polyRAD.Geno, m_train.pheno=com.taxa, her.vec=her.vec, dat=Msa.polyRAD.geno.2, PCA.info=PCA.Msa, Opt.PCA=Opt.PCA, 
                                   num.vec=num.vec, s=s, Iter=Iter, burn=burn, species=species, k=k, itrn=itrn, cycles = cycles,
                                   trn.app=trn.app, number.of.folds=number.of.folds, geno.diplo.tetra=geno.diplo.tetra, Groups.3.Pop3=Groups.3.Pop3, date=date)




