#############################################################################
# Script to optimize the calibration set in genomic selection (maximize the expected reliability).
# Method based on the generalized CD.
# (Rincent et al. 2012)
#############################################################################
# 28/08/2012, author renaud.rincent@inra.fr


###############
#Functions used
###############
rm(list=ls())

setwd("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/DiversityPanel.GS/Msa")

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")
library(heritability)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")
#Load your hapmap data and use GAPT to convert it to numeric format
data_hmp <- read.delim("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.Geno.5pmaf.90pmissing.hmp.txt", head=F)
myG <- data_hmp
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
dat=myGAPIT$GD

nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
dat[1:6,1:6]

rownames(dat)=nan
# This function creates the matrix of contrast between each of the individual not in the calibration set and the mean of the population
contrasteNonPheno=function(NotSampled_f,Nind_f,Nind_in_Sample_f)
{
mat=matrix(-1/Nind_f,Nind_f,Nind_f-Nind_in_Sample_f)
for (i in 1:ncol(mat)) {
mat[NotSampled_f[i],i]=1-1/Nind_f
}
return(mat)
}

##############################
# Data required
##########################
#Msa.Pheno <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.Pheno.csv", header=T)
Msa.Pheno <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/com.taxa.numeric.pheno.geno.Msi.Msa.csv", header=T)
Msa.Pheno$Taxa <- as.character(Msa.Pheno$Taxa)
library(rrBLUP)
phenames <- colnames(Msa.Pheno[,-1])


geno <- dat-1

impute=A.mat(geno,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed

geno.names <- data.frame(rownames(Markers_impute))
colnames(geno.names)[1] <-"Taxa"
geno.names$Taxa <- as.character(geno.names$Taxa)
com.tax <- merge(Msa.Pheno, geno.names, by="Taxa")

dat2 <- Markers_impute[match(com.tax$Taxa, rownames(Markers_impute)),]

set.seed(234)
train=as.matrix(sample(1:nrow(dat2),0.8*nrow(dat2)), replace=F) # number of genotypes selected out of total, here 80% --> chnage 0.8 if wanted
valid<-setdiff(1:nrow(dat2),train)

# taining set
Pheno_train=com.tax[train,] # phenotypes
m_train1=dat2[train,]
train.pop <- as.character(rownames(m_train1))
write.table(train.pop, "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/DiversityPanel.GS/Msa/Msa.Div.CDmean.CalibrationSet.txt", sep="\t", quote=F, row.names = F, col.names = F)


#validation set 
#Pheno_valid=Pheno[valid,] # phenotype
m_valid1=dat2[valid,]# 
test.pop <- as.character(rownames(m_valid1))
write.table(test.pop, "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/DiversityPanel.GS/Msa/Msa.Div.CDmean.TestSet.txt", sep="\t", quote=F, row.names = F, col.names = F)

Geno <- m_train1+1
matA1 <- A.mat(m_train1, min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,shrink=FALSE,return.imputed=FALSE)

#matA1=read.table("matA1Dent_sansPond.csv") #I think this is the GRM #This is the covariance matrix betw the individuals (size Nind x xNind), estimated with the genotypes.
matA1=as.matrix(matA1)

phen<-com.tax
A.A <- A.mat(dat2, impute.method="mean",return.imputed = TRUE) # additive relationship matrix

h2<-NULL
for(i in 2:dim(phen)[2]){
  y<-phen[,i]
  y[is.na(y)]<-mean(y,na.rm=T)
  model <- mixed.solve(y,K=A.A$A)
  h_sq<-model$Vu/(model$Vu+model$Ve)
  h2<-c(h2,h_sq)
  
}
h2<-data.frame(h2,traits=c(names(phen)[-1]))

write.table(h2, "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/DiversityPanel.GS/Msa/Msa.trait.heritablility.csv", sep=",", quote=F, row.names = F, col.names = T)


Nind=nrow(matA1) # total number of individuals
ind.names <- rownames(matA1)# names of individuals
itrn=3000# number of iterations

# Choose a size for your calibration set
num.vec <- c(50,100,150,200,250,300,350)


for(p in 1:length(phenames)){
  
      for(n in 1:length(num.vec)){
        print(paste("------ Calibration set size being evaluated : ", num.vec[n], " !!!!!!!!!!!---------------", sep = ""))
        
            y<-phen[,phenames[p]]
            y[is.na(y)]<-mean(y,na.rm=T)
            model <- mixed.solve(y,K=A.A$A)
            h_sq<-model$Vu/(model$Vu+model$Ve)
            
            
            nindrep=num.vec[n]
            varG = model$Vu
            varE <- model$Ve
            lambda=varE/varG # lambda is needed to estimate the CDmean # How do we estimate lambda in our work?
            
            #invA1=ginv(matA1) # Inverse of the covariance matrix
            invA1=ginv(matA1) #changes solve to ginv because it complained of singular matrix
            
            
            ##############################
            # Optimization algo
            ##############################
            
            Nind_in_Sample=nindrep
            
            #Design matrices
            Ident<-diag(Nind_in_Sample)
            X<-rep(1,Nind_in_Sample)
            M<-Ident- (X%*%ginv(t(X)%*%X) %*% t(X) )
            #dim(M)
            #[1] 100 100
            
            Sample1<-sample(Nind,Nind_in_Sample) #Calibration set initialization
            SaveSample=Sample1
            NotSampled1<-seq(1:Nind)
            NotSampled<-NotSampled1[-Sample1] # Initial validation set
            
            Z=matrix(0,Nind_in_Sample,Nind)
            
            for (i in 1:length(Sample1)) { Z[i,Sample1[i]]=1 } 
            
            T<-contrasteNonPheno(NotSampled,Nind,Nind_in_Sample)   # T matrice des contrastes
            
            # Calculate of CDmean of the initial set
            matCD<-(t(T)%*%(matA1-lambda*ginv(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
            
            dim(matCD)
            #[1] 454 454
            matCD[1:10,1:10]
            
            CD=diag(matCD)
            CDmeanSave=mean(CD)
            
            CDmeanMax1=rep(NA,itrn)
            
            # Exchange algorithm (maximize CDmean)
            cpt2=1
            cpt=0
            while (cpt2<itrn) {
              print(paste("Trait: ", phenames[p]," --- Trait: ",p, " out of: ", length(phenames)," --- Calibration set: ", num.vec[n],"--- Iteration: ", cpt2, " out of ", itrn," ---", sep = ""))            # Make sure that 800 is enough in your case (that you reached a plateau), for this look at CDmeanMax1.
             NotSampled=NotSampled1[-Sample1] 
            cpt2=cpt2+1
            # Remove one individual (randomly choosen) from the sample :
            Sample2=sample(Sample1,1)
            # Select one individual (randomly choosen) from the individuals that are not in the Calibration set :
            Sample3=sample(NotSampled,1)
            # New calibration set :
            Sample4=c(Sample3,Sample1[Sample1!=Sample2])
            # Calculate the mean CD of the new calibration set :
            Z=matrix(0,Nind_in_Sample,Nind)
            for (i in 1:length(Sample4)) { Z[i,Sample4[i]]=1 } 
            NotSampled=NotSampled1[-Sample4] 
            T<-contrasteNonPheno(NotSampled,Nind,Nind_in_Sample)
            
            matCD<-(t(T)%*%(matA1-lambda*ginv(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
            CD=diag(matCD)
            
            if (mean(CD)>CDmeanSave ) { Sample1=Sample4 # Accept the new Calibration set if CDmean is increased, reject otherwise.
            CDmeanSave=mean(CD)  
            cpt=0 } else { cpt=cpt+1 
            }
            CDmeanMax1[cpt2-1]=CDmeanSave
            }  #Fin du while
            
            SampleOptimiz=Sample1 # SampleOptimiz is the optimized calibration set
            write.table(SampleOptimiz, paste("Msa.CDmeanMax1.DivPanel.3000", phenames[p],"size-n:",num.vec[n],"txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
            SampleOptimiz.names <- ind.names[SampleOptimiz]
            write.table(SampleOptimiz.names, paste("Msa.CDmeanMax1.DivPanel.itrn3000.individuals", phenames[p],"size-n:",num.vec[n],"txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
            
            # End
            pdf(paste("Msa.CDmeanMax1.DivPanel.3000", phenames[p],"size-n:",num.vec[n],"pdf", sep="."))
            plot(c(1:itrn), CDmeanMax1)
            dev.off()
            
      
      }

}



# Five fold cross validation
rm(list=ls())
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")
library(heritability)
#library(kernlab)
#library(BGLR)
#install.packages("sommer")
#install.packages("BGLR")
#install.packages("rrBLUP")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")

specie.name <- c("Msa", "Msi")


geno.all <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/combined_genotypes.csv", header=T)
geno.all[1:6,1:6]
name.geno <- as.vector(geno.all$X)
df.nam.geno <- data.frame(name.geno)
names(df.nam.geno) <- "Taxa"
head(df.nam.geno)
df.nam.geno$Taxa <- as.character(df.nam.geno$Taxa)

#write.table(com.taxa, "com.taxa.numeric.pheno.geno.Msi.Msa.csv", sep=",", quote=F, row.names = F, col.names = T)
phen.taxa <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/com.taxa.numeric.pheno.geno.Msi.Msa.csv", header=T)
phen.taxa$Taxa <- as.character(phen.taxa$Taxa)
#phen.nam.taxa <- data.frame(phen.taxa$Taxa)
#colnames(phen.nam.taxa) <- "Taxa"
com.taxa <- merge(phen.taxa, df.nam.geno)
dim(com.taxa)
head(com.taxa)
com.taxa$Taxa <- as.character(com.taxa$Taxa)
#geno.all <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/combined_genotypes.csv", header=T)

com.geno <- geno.all[match(com.taxa$Taxa, geno.all$X),]
com.geno[1:6,1:6]
dim(com.geno)
com.geno2 <- com.geno[,-1]
rownames(com.geno2) <- as.character(com.geno$X)

test.set.CDmean <- read.delim("/Users/omo/Desktop/PostDoc/Manuscript2/Data/Msa.Div.CDmean.TestSet.txt", header=F)
head(test.set.CDmean)
test.set.CDmean$V1 <- as.character(test.set.CDmean$V1)
colnames(test.set.CDmean) <- "Taxa"
head(test.set.CDmean)
test.set.CDmean.Pheno <- merge(test.set.CDmean, com.taxa)

m_valid1.pheno <- test.set.CDmean.Pheno
m_valid1.polyRAD <- com.geno2[match(test.set.CDmean.Pheno$Taxa, rownames(com.geno2)),]


trainingset.CDmean <- read.delim("/Users/omo/Desktop/PostDoc/Manuscript2/Data/Msa.Div.CDmean.CalibrationSet.txt", header=F)
head(trainingset.CDmean)
trainingset.CDmean$V1 <- as.character(trainingset.CDmean$V1)
colnames(trainingset.CDmean)<- "Taxa"
head(trainingset.CDmean)
trainingset.CDmean.Pheno <- merge(trainingset.CDmean, com.taxa)

m_train.pheno <- trainingset.CDmean.Pheno
m_train1.polyRAD <- com.geno2[match(trainingset.CDmean.Pheno$Taxa, rownames(com.geno2)),]

phenames <- as.character(colnames(com.taxa[,-1]))

her.vec <- read.csv('Msa.trait.heritablility.csv', header=T)
her.vec$trait <- as.character(her.vec$trait)

num.vec <- c(50,100,150,200,250,300,350)
s <- 0.2  #  enter selection index here (proportion to select)                 
Iter=2500
burn=1000

# m_train1
cycles = 20
k = 5

#A.Gen <- A.mat()
#out <- marker_h2(Pheno.trt[,phenames[p]], Pheno_train$Taxa, covariates = NULL, A.Gen, alpha = 0.05, eps = 1e-06, max.iter = 100, fix.h2 = FALSE)


GS.Results.MsiMsa.TRSvec=NULL
GS.Results.MsiMsa.TRSvec.Traits=NULL


for(l in 1:length(phenames)){
  
  Pheno1=as.matrix(m_train.pheno[,phenames[l]])
  rownames(Pheno1) <- m_train.pheno$Taxa
  colnames(Pheno1) <- phenames[l]
    
  G <- m_train1.polyRAD
  herit <- her.vec[which(her.vec$trait==phenames[l]),1]
  
  
  for(n in 1:length(num.vec)){
        print(paste("-------------- Trait being analysed: ", phenames[l], " TRS.Size: ", num.vec[n], " !!!!!!!!!!!---------------", sep = ""))
        Rx.Trts <- NULL
        Ro.Trts <- NULL
        Rk.Trts <- NULL
        Sv.Trts <- NULL
        Ba.Trts <- NULL
        
        CDmean.names.Msa <- read.delim(paste(specie.name[1],"CDmeanMax1.itrn3000.individuals", phenames[l],"size-n:",num.vec[n],"txt", sep="."), header=F)
        #CDmean.names.Msi <- read.delim(paste(specie.name[2],"CDmeanMax1.itrn3000.individuals", phenames[l],"size-n:",num.vec[n],"txt", sep="."), header=T)
        
        trt.vec1 <- as.character(CDmean.names.Msa$V1)
        #trt.vec2 <- as.vector(CDmean.names.Msi[,1])
        trt.vec12 <- trt.vec1
        trt.vec.df <- data.frame(trt.vec12)
        colnames(trt.vec.df)[1] <- "Taxa"
        trt.vec.df$Taxa <- as.character(trt.vec.df$Taxa)
        #G2 <- G[match(trt.vec, rownames(G)),]
        G.df <- data.frame(rownames(G))
        colnames(G.df)[1] <- "Taxa"
        G.df$Taxa <- as.character(G.df$Taxa)
        
        com.TRS.geno <- merge(trt.vec.df, G.df, by="Taxa")
        trt.vec <- as.character(com.TRS.geno$Taxa)
        
        G2 <- G[match(com.TRS.geno$Taxa, rownames(G)),]
         #+ 1
        
        Pheno1.TRS <- as.matrix(Pheno1[match(trt.vec, rownames(Pheno1)),])
        rownames(Pheno1.TRS) <- trt.vec
        colnames(Pheno1.TRS) <- phenames[l]
        mean.trs <- mean(Pheno1.TRS[,1], na.rm=T)
        var.trs <- var(Pheno1.TRS[,1], na.rm=T)
    

      for(j in 1:cycles){
        
            Rx.accuracy=matrix(nrow=cycles,ncol=k+20)# pred acc, intercept, slope, ci, and reliability for 5 cross folds in each cycle
            Ro.accuracy=matrix(nrow=cycles,ncol=k+20)
            Rk.accuracy=matrix(nrow=cycles,ncol=k+20)
            Ba.accuracy=matrix(nrow=cycles,ncol=k+20)
            
            r.gy <- NULL
            r.gy.E <- NULL
            r.gy.BA <- NULL
            r.gy.RK <- NULL
            
            
            myCI.RRB <- NULL
            myCI.E <- NULL
            myCI.BA <- NULL
            myCI.RK <- NULL
            
            
            the.coefficients.Ro <- NULL
            the.coefficients.E <- NULL
            the.coefficients.BA <- NULL
            the.coefficients.RK <- NULL
            
            
            Pheno1.TRS <- as.data.frame(Pheno1.TRS)
            Pheno1.TRS$id <- sample(1:k, nrow(Pheno1.TRS), replace = TRUE)
            
            Pheno2.TS <- as.matrix(m_valid1.pheno[,phenames[l]])
            rownames(Pheno2.TS) <- m_valid1.pheno$Taxa
            colnames(Pheno2.TS) <- phenames[l]
            Pheno2.Valid <- Pheno2.TS
            Pheno2.Valid <- as.data.frame(Pheno2.Valid)
            Pheno2.Valid$id <- sample(1:k, nrow(Pheno2.Valid), replace = TRUE)
            
            list <- 1:k
        
        for (i in 1:k){
          
          # remove rows with id i from dataframe to create training set
          # select rows with id i to create test set
          pheno_trainingset <- subset(Pheno1.TRS, id %in% list[-i])
          
          #colnames(pheno_trainingset)[1] <- phenames[l]
          Geno_trainingset <- G2[match(rownames(pheno_trainingset), rownames(G2)),]
          Geno_trainingset <- as.matrix(Geno_trainingset)
          pheno_testset <- subset(Pheno2.Valid, id %in% list[-i])
          #colnames(pheno_testset) <- phenames[l]
          Geno_testset <- m_valid1.polyRAD[match(rownames(pheno_testset), rownames(m_valid1.polyRAD)),]
          Geno_testset <- as.matrix(Geno_testset)
          # run RRBLUP
          Yt=(pheno_trainingset[,1])
          #LPBL_answer<-mixed.solve(LPBL, Z=Geno_trainingset, K=NULL, SE = FALSE, return.Hinv=FALSE)
          
          rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=Geno_trainingset, K=NULL, SE=F, return.Hinv=F)
          mEff.Ro<-rrMod.Ro$u
          e.Ro= as.matrix(mEff.Ro)
          
          predYv.Ro = Geno_testset%*%e.Ro
          
          predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
          
          Y_valid=pheno_testset[,1]
          names(Y_valid) <- rownames(pheno_testset)
          #Ro.accuracy[j,i] <- cor(predYr.Ro,Y_valid,use="complete")
          r.gy <- c(r.gy , cor(predYr.Ro,Y_valid, use="complete" ))
          the.fitted.model.Ro <- lm(Y_valid ~ predYr.Ro)
          the.coefficients.Ro <- c(the.coefficients.Ro, the.fitted.model.Ro$coefficients[1], the.fitted.model.Ro$coefficients[2])
          
          Ro.accuracy[j,i+5] <- the.coefficients.Ro[1]
          Ro.accuracy[j,i+10] <- the.coefficients.Ro[2]
          
          these.observed.and.predicted.phenotypic.values.Ro <- data.frame(names(Y_valid), Y_valid, predYr.Ro)
          rownames(these.observed.and.predicted.phenotypic.values.Ro) <- NULL
          colnames(these.observed.and.predicted.phenotypic.values.Ro) <- c("Taxa", "Observed.Value", "Predicted.Value")
          x.p.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,3)]
          y.o.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,2)]
          Ro.accuracy[j,i+15] <- round(CI(x.p.Ro,y.o.Ro,s=s,top=T),2)
          
          corr.sqr <- round(summary(the.fitted.model.Ro)$adj.r.squared, 2)
          gsr <- round((corr.sqr/herit), 2)
          Ro.accuracy[j,i+20] <- gsr
          
          print(paste("--------- Trait being analysed: ", phenames[l], " TRS.Size: ", num.vec[n], " --cycle-- ", j, " ; GS Model: ADE"," !!!!!!!----------", sep = ""))
          library(sommer)
          
          Pheno.comb <- rbind(pheno_trainingset, pheno_testset)
          train.comb <- 1:nrow(pheno_trainingset)
          valid.comb <-setdiff(1:nrow(Pheno.comb),train.comb)
          
          y.trn <- Pheno.comb # for prediction accuracy
          ww <- pheno_testset[,-2] # delete data for 1/5 of the population
          names(ww) <- rownames(pheno_testset)
          #y.trn[y.trn%in%Pheno2] <- NA
          y.trn[valid.comb,1] <- NA #valid.comb
          Geno.comb <- rbind(Geno_trainingset, Geno_testset)
          Geno.comb <- as.matrix(Geno.comb)
          M <-tcrossprod(Geno.comb)/ncol(Geno.comb)
          X <- Geno.comb
          
          #RKHS 
          library(BGLR)
          print(paste("--------- Trait being analysed: ", phenames[l], " TRS.Size: ", num.vec[n], " --cycle-- ", j, " ; GS Model: RKHS"," !!!!!!!----------", sep = ""))
          ETA<-list(list(K=M,model='RKHS')) 
          fm.RK<-BGLR(y=y.trn[,1],ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
          Rk.accuracy[j,i] <- cor(fm.RK$yHat[valid.comb], ww, use="complete")
          the.fitted.model.RK <- lm(ww ~ fm.RK$yHat[valid.comb])
          the.coefficients.RK <- c(the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
          Rk.accuracy[j,i+5] <- the.coefficients.RK[1]
          Rk.accuracy[j,i+10] <- the.coefficients.RK[2]
          
          these.observed.and.predicted.phenotypic.values.Rk <- data.frame(names(Y_valid), Y_valid, fm.RK$yHat[valid.comb])
          rownames(these.observed.and.predicted.phenotypic.values.Rk) <- NULL
          colnames(these.observed.and.predicted.phenotypic.values.Rk) <- c("Taxa", "Observed.Value", "Predicted.Value")
          x.p.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,3)]
          y.o.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,2)]
          Rk.accuracy[j,i+15] <- round(CI(x.p.Rk,y.o.Rk,s=s,top=T),2)
          these.observed.and.predicted.phenotypic.values.Rk <- these.observed.and.predicted.phenotypic.values.Rk[order(these.observed.and.predicted.phenotypic.values.Rk$Predicted.Value, decreasing = T),]
         # sel.taxa.Rk <- these.observed.and.predicted.phenotypic.values.Rk$Taxa[c(1:(round(length(these.observed.and.predicted.phenotypic.values.Rk$Taxa)*0.1,0)))]
         # sel.taxa.Rk <- as.character(sel.taxa.Rk)
          #sel.taxa.Rk.numvec <- c(sel.taxa.Rk.numvec, sel.taxa.Rk)
          corr.sqr.Rk <- round(summary(the.fitted.model.RK)$adj.r.squared, 2)
          gsr.Rk <- round((corr.sqr.Rk/herit), 2)
          Rk.accuracy[j,i+20] <- gsr.Rk
          
          #Bayes A
          print(paste("--------- Trait being analysed: ", phenames[l], " TRS.Size: ", num.vec[n], " --cycle-- ", j, " ; GS Model: Bayes C"," !!!!!!!----------", sep = ""))
          ETA<-list(list(X=X,model='BayesA')) 
          fm.BA<-BGLR(y=y.trn[,1],ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
          Ba.accuracy[j,i] <- cor(fm.BA$yHat[valid.comb], ww, use="complete")
          the.fitted.model.BA <- lm(ww ~ fm.BA$yHat[valid.comb])
          the.coefficients.BA <- c(the.fitted.model.BA$coefficients[1], the.fitted.model.BA$coefficients[2])
          Ba.accuracy[j,i+5] <- the.coefficients.BA[1]
          Ba.accuracy[j,i+10] <- the.coefficients.BA[2]
          
          these.observed.and.predicted.phenotypic.values.Ba <- data.frame(names(Y_valid), Y_valid, fm.BA$yHat[valid.comb])
          rownames(these.observed.and.predicted.phenotypic.values.Ba) <- NULL
          colnames(these.observed.and.predicted.phenotypic.values.Ba) <- c("Taxa", "Observed.Value", "Predicted.Value")
          x.p.Ba=these.observed.and.predicted.phenotypic.values.Ba[,c(1,3)]
          y.o.Ba=these.observed.and.predicted.phenotypic.values.Ba[,c(1,2)]
          Ba.accuracy[j,i+15] <- round(CI(x.p.Ba,y.o.Ba,s=s,top=T),2)
          #these.observed.and.predicted.phenotypic.values.Ba <- these.observed.and.predicted.phenotypic.values.Ba[order(these.observed.and.predicted.phenotypic.values.Ba$Predicted.Value, decreasing = T),]
          #sel.taxa.Ba <- these.observed.and.predicted.phenotypic.values.Ba$Taxa[c(1:(round(length(these.observed.and.predicted.phenotypic.values.Ba$Taxa)*0.1,0)))]
          #sel.taxa.Ba <- as.character(sel.taxa.Ba)
          #sel.taxa.Ba.numvec <- c(sel.taxa.Ba.numvec, sel.taxa.Ba)
          corr.sqr.Ba <- round(summary(the.fitted.model.BA)$adj.r.squared, 2)
          gsr.Ba <- round((corr.sqr.Ba/herit), 2)
          Ba.accuracy[j,i+20] <- gsr.Ba
          
        }
            


      }
        Ro.accuracy <- data.frame(Ro.accuracy)
        #Ro.accuracy$Trt <- rep(phenames[l], nrow(Ro.accuracy))
        Ro.Trts <- rbind(Ro.Trts, Ro.accuracy)
        
        
        Rk.accuracy <- data.frame(Rk.accuracy)
        #Rk.accuracy$Trt <- rep(phenames[l], nrow(Rk.accuracy))
        Rk.Trts <- rbind(Rk.Trts, Rk.accuracy)
        
        Ba.accuracy <- data.frame(Ba.accuracy)
        #Ba.accuracy$Trt <- rep(phenames[l], nrow(Ba.accuracy))
        Ba.Trts <- rbind(Ba.Trts, Ba.accuracy)
        
        
        
        Ro.Trts$Model <- rep("RRBLUP", nrow(Ro.Trts))
        Rk.Trts$Model <- rep("RKHS", nrow(Rk.Trts))
        Ba.Trts$Model <- rep("BA", nrow(Ba.Trts))
        
        
        GS.Results.MsiMsa <- rbind(Ro.Trts, Rk.Trts, Ba.Trts) #Rx.Trts, 
        GS.Results.MsiMsa$TRS.Size <- rep(num.vec[n], nrow(GS.Results.MsiMsa))
        GS.Results.MsiMsa.TRSvec <- rbind(GS.Results.MsiMsa.TRSvec, GS.Results.MsiMsa)
        
        write.table(GS.Results.MsiMsa , paste("GS.Results.Msa", phenames[l], num.vec[n], ".txt",sep="_"), sep="\t", quote=F, row.names = F, col.names = F)
        
        GS.Results.MsiMsa=NULL
        
   }
    
      print(paste("--------- Summarizing Result for trait: ", phenames[l]," !!!!!!!----------", sep = ""))
      
      GS.Results.MsiMsa.TRSvec$Trt <- rep(phenames[l], nrow(GS.Results.MsiMsa.TRSvec))
      write.table(GS.Results.MsiMsa.TRSvec , paste("GS.Results.Msa.TRSvec", phenames[l], ".txt",sep="_"), sep="\t", quote=F, row.names = F, col.names = F)
      
      GS.Results.MsiMsa.TRSvec.Traits <- rbind(GS.Results.MsiMsa.TRSvec.Traits, GS.Results.MsiMsa.TRSvec)
      GS.Results.MsiMsa.TRSvec=NULL
      
  }

save.image("GS.GS.Results.Msa.TRSvec.rda")
write.table(GS.Results.MsiMsa.TRSvec.Traits, "GS.Results.Msa.TRSvec.Traits.03172019.csv", sep=",", quote=F, row.names=F, col.names=T)






