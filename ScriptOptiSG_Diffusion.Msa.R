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
Msa.Pheno <- read.csv("http://people.beocat.ksu.edu/~omo/effect_simulations/PostDoc/Manuscript2/Data/Msa.2x_4x.Pheno.csv", header=T)

library(rrBLUP)
phenames <- colnames(Msa.Pheno[,-1])

dat[1:6,1:6]
geno <- dat-1

impute=A.mat(geno,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed

geno.names <- data.frame(rownames(Markers_impute))
colnames(geno.names)[1] <-"Taxa"
com.tax <- merge(Msa.Pheno, geno.names, by="Taxa")

dat <- dat[match(com.tax$Taxa, rownames(dat)),]

Geno <- Markers_impute+1
matA1 <- A.mat(Markers_impute, min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,shrink=FALSE,return.imputed=FALSE)

#matA1=read.table("matA1Dent_sansPond.csv") #I think this is the GRM #This is the covariance matrix betw the individuals (size Nind x xNind), estimated with the genotypes.
matA1=as.matrix(matA1)

A <- A.mat(Markers_impute) # additive relationship matrix

Nind=nrow(matA1) # total number of individuals
ind.names <- rownames(matA1)# names of individuals
itrn=3000# number of iterations

# Choose a size for your calibration set
num.vec <- c(50,100,150,200,250,300,350)

for(p in 1:length(phenames)){
  out <- marker_h2(com.tax[,phenames[p]], com.tax$Taxa, covariates = NULL, A, alpha = 0.05, eps = 1e-06, max.iter = 100, fix.h2 = FALSE)
  
      for(n in 1:length(num.vec)){
        print(paste("------ Calibration set size being evaluated : ", num.vec[n], " !!!!!!!!!!!---------------", sep = ""))
        
            nindrep=num.vec[n]
            varG = out$va
            varE <- out$ve
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
              print(paste("Trait: ", phenames[p]," --- Trait: ",p, " out of: ", length(phenames)," --- Calibration set: ", num.vec[n],"--- Iteration: ", cpt2, " out of ", itrn," ---", sep = ""))
            # Make sure that 800 is enough in your case (that you reached a plateau), for this look at CDmeanMax1.
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
            write.table(SampleOptimiz, paste("Msa.CDmeanMax1.3000", phenames[p],"size-n:",num.vec[n],"txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
            SampleOptimiz.names <- ind.names[SampleOptimiz]
            write.table(SampleOptimiz.names, paste("Msa.CDmeanMax1.itrn3000.individuals", phenames[p],"size-n:",num.vec[n],"txt", sep="."), sep="\t", quote=F, row.names=F, col.names=F)
            
            # End
            pdf(paste("Msa.CDmeanMax1.3000", phenames[p],"size-n:",num.vec[n],"pdf", sep="."))
            plot(c(1:itrn), CDmeanMax1)
            dev.off()
            
      
      }

}

