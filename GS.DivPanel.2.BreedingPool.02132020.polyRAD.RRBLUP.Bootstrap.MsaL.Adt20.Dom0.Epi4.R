rm(list=ls())
#setwd("")
date="02132020"
setwd("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/09F2/simulations.polyRAD/")
hm.dir <- "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/09F2/simulations.polyRAD/"
par.dir <- "Diff.QTNs.Diff.Eff.MsaL/"

par.msi.dir <- "Msi_20_Add_QTN_1_Dom_QTN_2_Epi_QTN_h.2_0.6_reps_1/"
par.msa.dir <- "Msa_20_Add_QTN_1_Dom_QTN_2_Epi_QTN_h.2_0.6_reps_1/"
f2.dir <- "F2.Pop1/"

ga.dir <- "F2.Traits.Adt20.Dom0.Epi4.Traits/"
species <- "MsaL.Adt20.Dom0.Epi4"

phename <- "Adt20.Dom0.Epi4"

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
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Diversity.Panel.To.09F2.GS.Objective2.polyRAD.RRBLUP.Algorithm.02132020.R")



Msi.polyRAD.geno <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/geno.msi.polyRAD.txt", header = T)
rownames(Msi.polyRAD.geno) <- Msi.polyRAD.geno$X
Msi.polyRAD.geno <- Msi.polyRAD.geno[,-1]
Msi.polyRAD.geno <- as.matrix(Msi.polyRAD.geno)
Msi.polyRAD.geno[1:6,1:6]
Msi.polyRAD.geno.2 <- Msi.polyRAD.geno*2
Msi.polyRAD.geno.2 <- round(Msi.polyRAD.geno.2)


nan.msi=rownames(Msi.polyRAD.geno.2)
nan.msi <- data.frame(nan.msi)
colnames(nan.msi)[1] <- "Taxa"

Msa.polyRAD.geno <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/geno.msa.polyRAD.txt", header = T)
rownames(Msa.polyRAD.geno) <- Msa.polyRAD.geno$X
Msa.polyRAD.geno <- Msa.polyRAD.geno[,-1]
Msa.polyRAD.geno <- as.matrix(Msa.polyRAD.geno)
Msa.polyRAD.geno[1:6,1:6]
Msa.polyRAD.geno.2 <- Msa.polyRAD.geno*2
Msa.polyRAD.geno.2 <- round(Msa.polyRAD.geno.2)
Msa.polyRAD.geno.2[1:6,1:6]

nan.msa=rownames(Msa.polyRAD.geno.2)
nan.msa <- data.frame(nan.msa)
colnames(nan.msa)[1] <- "Taxa"


Msi.div <- read.delim(paste(hm.dir, par.dir, ga.dir, par.msi.dir, "Simulated.Data.1.Reps.Herit.0.6.txt", sep=""), header=T)
colnames(Msi.div)[1] <- "Taxa"
colnames(Msi.div)[2] <- phename
Msi.div$Taxa <- as.character(Msi.div$Taxa)

Msi.Div.RawPheno <- Msi.div
msi.div.taxa <- Msi.div$Taxa

Msa.div <- read.delim(paste(hm.dir, par.dir, ga.dir, par.msa.dir, "Simulated.Data.1.Reps.Herit.0.6.txt", sep=""), header=T)
colnames(Msa.div)[1] <- "Taxa"
colnames(Msa.div)[2] <- phename
Msa.div$Taxa <- as.character(Msa.div$Taxa)

Msa.Div.RawPheno <- Msa.div
msa.div.taxa <- Msa.div$Taxa

# Read in Msi data
PCA.Msi <- read.table("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msi.polyRAD.estimated.PCA.txt", header=T)
PCA.Msi[1:6,1:6]
names(PCA.Msi)[1] <- "Taxa" # Select 7 PCs as optimal
PCA.Msi <- PCA.Msi[,c(1:8)]

PCA.Msa <- read.table("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.polyRAD.estimated.PCA.txt", header=T)
PCA.Msa[1:6,1:6]
PCA.Msa <- PCA.Msa[,c(1:6)] # Select 5 PCs

phenames <- phename



# Read in 09F2 Data
geno.09F2 <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/F2.PolyRadGeno.216lines.5141snps.11052019.csv", header=T)
geno.F2 <- geno.09F2[,-1]
rownames(geno.F2) <- geno.09F2$X
geno.F2 <- as.matrix(geno.F2)
valid.geno <- geno.F2

Phenotypes.09F2 <- read.delim(paste(hm.dir, par.dir, ga.dir, f2.dir, "Simulated.Data.1.Reps.Herit.0.37.txt", sep=""), header=T)
colnames(Phenotypes.09F2)[1] <- "Taxa"
colnames(Phenotypes.09F2)[2] <- phename

valid.pheno <- Phenotypes.09F2

phen<-Phenotypes.09F2

her.vec <- read.csv('http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/F2.trait.heritablility.csv', header=T)
her.vec$trait <- as.character(her.vec$traits)


num.vec <- 200 #,150,200,250,300,350
f2.vec <- 216
bootstraps.vec=1000
s <- 0.2  #  enter selection index here (proportion to select)                 
Iter=2500
burn=1000

#either Msi or Msa
# m_train1
k = 5
itrn=3000

trn.app <- c( "Msi.Random", "Msa.Random", "Msi.Div", "Msa.Div", "Msi.Div.CDMean", "Msa.Div.CDMean", "msi.plus.msa")

number.of.folds=k

DiversityPanel.2.BreedingPool.GS.Pipeline(phenames=phenames, valid.geno=valid.geno, f2.vec=f2.vec, bootstraps.vec=bootstraps.vec,
                                          valid.pheno=valid.pheno, Msi.div.taxa=msi.div.taxa, Msa.div.taxa=msa.div.taxa, her.vec=her.vec, PCA.info.Msi=PCA.Msi, PCA.info.Msa=PCA.Msa,
                                          dat.msi=Msi.polyRAD.geno.2, dat.msa=Msa.polyRAD.geno.2, num.vec=num.vec, s=s, Iter=Iter, burn=burn, species=species, k=k, itrn=itrn, 
                                          trn.app=trn.app, number.of.folds=number.of.folds, Msi.Div.RawPheno=Msi.Div.RawPheno, Msa.Div.RawPheno=Msa.Div.RawPheno, date=date, 
                                          Msi.polyRAD.geno=Msi.polyRAD.geno, Msa.polyRAD.geno=Msa.polyRAD.geno)


phenames=phenames; valid.geno=valid.geno; f2.vec=f2.vec;bootstraps.vec=bootstraps.vec;
valid.pheno=valid.pheno; Msi.div.taxa=msi.div.taxa; Msa.div.taxa=msa.div.taxa; her.vec=her.vec; PCA.info.Msi=PCA.Msi; PCA.info.Msa=PCA.Msa;
dat.msi=Msi.polyRAD.geno.2; dat.msa=Msa.polyRAD.geno.2; num.vec=num.vec; s=s; Iter=Iter; burn=burn; species=species; k=k; itrn=itrn; 
trn.app=trn.app; number.of.folds=number.of.folds; Msi.Div.RawPheno=Msi.Div.RawPheno; Msa.Div.RawPheno=Msa.Div.RawPheno; date=date




