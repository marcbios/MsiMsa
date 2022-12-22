rm(list=ls())
#setwd("")
date="02112020"
setwd("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/09F2/PartialQTNs/F2.Traits.Adt20.Dom4.Epi0.Traits/")
home.dir <- "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/09F2/PartialQTNs/F2.Traits.Adt20.Dom4.Epi0.Traits/"
Traitname <- "A20D4E0"
msi.dir <- "Msi_20_Add_QTN_4_Dom_QTN_1_Epi_QTN_h.2_0.6_reps_1/"
msa.dir <- "Msa_20_Add_QTN_4_Dom_QTN_1_Epi_QTN_h.2_0.6_reps_1/"
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
source("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Diversity.Panel.To.09F2.GS.Objective2.RRBLUP.AD.Algorithm.SimulatedF2.12122019.R")

Msi.polyRAD.geno <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/Msi.584Inds.356SNPs.genotypes.impute.genotypes.txt", header = T)
Msi.polyRAD.geno[1:6,1:6]
rownames(Msi.polyRAD.geno) <- Msi.polyRAD.geno$Snp
Msi.polyRAD.geno <- Msi.polyRAD.geno[,-c(1:5)]
Msi.polyRAD.geno <- t(Msi.polyRAD.geno)
Msi.polyRAD.geno <- as.matrix(Msi.polyRAD.geno)
Msi.polyRAD.geno[1:6,1:6]

nan.msi=rownames(Msi.polyRAD.geno)
nan.msi <- data.frame(nan.msi)
colnames(nan.msi)[1] <- "Taxa"

Msa.polyRAD.geno <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Simulated.F2/Msa.643Inds.356SNPs.genotypes.impute.genotypes.txt", header = T)
rownames(Msa.polyRAD.geno) <- Msa.polyRAD.geno$Snp
Msa.polyRAD.geno <- Msa.polyRAD.geno[,-c(1:5)]
Msa.polyRAD.geno <- t(Msa.polyRAD.geno)
Msa.polyRAD.geno <- as.matrix(Msa.polyRAD.geno)
Msa.polyRAD.geno[1:6,1:6]

nan.msa=rownames(Msa.polyRAD.geno)
nan.msa <- data.frame(nan.msa)
colnames(nan.msa)[1] <- "Taxa"

Msi.div <- read.delim(paste(home.dir, msi.dir, "Simulated.Data.1.Reps.Herit.0.6.txt", sep=""), header=T)
colnames(Msi.div)[1] <- "Taxa"
colnames(Msi.div)[2] <- Traitname
Msi.div$Taxa <- as.character(Msi.div$Taxa)
Msi.Div.RawPheno <- Msi.div
msi.div.taxa <- Msi.div$Taxa

Msa.div <- read.delim(paste(home.dir, msa.dir, "Simulated.Data.1.Reps.Herit.0.6.txt", sep=""), header=T)
colnames(Msa.div)[1] <- "Taxa"
colnames(Msa.div)[2] <- Traitname
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

phenames <- Traitname


F2.dirs <- c(paste("F2.Pop", 1:100, "/",sep=""))
Compliled.Results.Pop.NULL=NULL

for(d in 1:length(F2.dirs)){
  
  print(paste("--------- Running Analysis....in simulated F2 Pop ", d, " !!!!!!!----------", sep = ""))
      # Read in 09F2 Data
      geno.09F2 <- read.delim(paste(home.dir, F2.dirs[d], "09F2pop.geno.txt", sep=""), header=T)
      rownames(geno.09F2) <- geno.09F2$Snp
      geno.F2 <- geno.09F2[,-c(1:5)]
      geno.F2 <- t(geno.F2)
      valid.geno <- geno.F2[c(1:216),]
      valid.geno[1:6,1:6]
      rownames(valid.geno) <- gsub("X", "", rownames(valid.geno))
      valid.geno[1:6,1:6]
      
      Phenotypes.09F2 <- read.delim(paste(home.dir, F2.dirs[d], "Simulated.Data.1.Reps.Herit.0.37.txt", sep=""), header=T)
      head(Phenotypes.09F2)
      colnames(Phenotypes.09F2)[1] <- "Taxa"
      colnames(Phenotypes.09F2)[2] <- Traitname
      valid.pheno <- Phenotypes.09F2[c(1:216),]
      
      num.vec <- 200 #,150,200,250,300,350
      f2.vec <- nrow(valid.pheno)
      bootstraps.vec=1000
      s <- 0.2  #  enter selection index here (proportion to select)                 
      species <- paste("F2.Pop", d, sep="_") #either Msi or Msa
      # m_train1
      k = 5
      itrn=3000
      
      trn.app <- c( "Msi.Random", "Msa.Random", "Msi.Div", "Msa.Div", "Msi.Div.CDMean", "Msa.Div.CDMean", "msi.plus.msa")
      #trn.app <- c("Msi.Div", "Msa.Div", "Msi.Div.CDMean", "Msa.Div.CDMean", "msi.plus.msa")
      
      number.of.folds=k
      
      Compliled.Results.Pop <- DiversityPanel.2.BreedingPool.GS.Pipeline(phenames=phenames, Msi.polyRAD.geno=Msi.polyRAD.geno, Msa.polyRAD.geno=Msa.polyRAD.geno, 
                                                Msi.Div.RawPheno=Msi.Div.RawPheno, Msa.Div.RawPheno=Msa.Div.RawPheno, 
                                                valid.geno=valid.geno, valid.pheno=valid.pheno, f2.vec=f2.vec,
                                                Msi.div.taxa=msi.div.taxa, Msa.div.taxa=msa.div.taxa, 
                                                PCA.info.Msi=PCA.Msi, PCA.info.Msa=PCA.Msa,
                                                num.vec=num.vec, s=s, species=species, k=k, itrn=itrn, 
                                                trn.app=trn.app, number.of.folds=number.of.folds, date=date)
      write.table(Compliled.Results.Pop, paste("Compliled.Results.Pop.Simulated.F2.Pop.PartialQTNs",d,Traitname,date,"csv", sep="."), quote=F, sep=",", row.names = F, col.names = T)
    
      dim(Compliled.Results.Pop)
      head(Compliled.Results.Pop)
      Compliled.Results.Pop$F2.Family <- rep(d, nrow(Compliled.Results.Pop))
      Compliled.Results.Pop.NULL <- rbind(Compliled.Results.Pop.NULL, Compliled.Results.Pop)
      
}

write.table(Compliled.Results.Pop.NULL, paste("Compliled.Results.Across.All.Simulated.F2.Pop.PartialQTNs",Traitname,date,"csv", sep="."), quote=F, sep=",", row.names = F, col.names = T)



#phenames=phenames; Msi.polyRAD.geno=Msi.polyRAD.geno; Msa.polyRAD.geno=Msa.polyRAD.geno; Msi.Div.RawPheno=Msi.Div.RawPheno; Msa.Div.RawPheno=Msa.Div.RawPheno; 
#valid.geno=valid.geno; valid.pheno=valid.pheno; f2.vec=f2.vec; bootstraps.vec=bootstraps.vec; Msi.div.taxa=msi.div.taxa; Msa.div.taxa=msa.div.taxa; 
#PCA.info.Msi=PCA.Msi; PCA.info.Msa=PCA.Msa; num.vec=num.vec; s=s; species=species; k=k; itrn=itrn; 
#trn.app=trn.app; number.of.folds=number.of.folds; date=date


