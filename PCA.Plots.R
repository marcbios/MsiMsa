#Use this to install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt", version = "3.8")
################
library(gdsfmt)
library(SNPRelate)

# Msa
vcf.fn <- "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/09F2/Msa/Msa.Geno.5pmaf.vcf"
snpgdsVCF2GDS(vcf.fn, "pm.gds", method="copy.num.of.ref")
snpgdsSummary("pm.gds")

# open the GDS file
genofile <- openfn.gds("pm.gds")
## LD based pruning of SNPs
set.seed(1000)
snpset2 <- snpgdsLDpruning(genofile, method = c("corr"), ld.threshold = 0.5,missing = 0.8, maf = 0.02)
snpset2.id <- unlist(snpset2)

########## PCA Analysis ###############
# get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#PCA analysis

#pm_pca<-snpgdsPCA(genofile,autosome.only=F)
pca.snpset2 <- snpgdsPCA(genofile, snp.id=snpset2.id, num.thread=2)

write.table(pca.snpset2$eigenvect, "Msa.pca.snpset2$eigenvect.txt")
write.table(pca.snpset2$eigenval, "Msa.pca.snpset2$eigenval.txt")
write.table(pca.snpset2$sample.id,"Msa.pca.snpset2$sample.id.txt")

tab <- data.frame(sample.id = pca.snpset2$sample.id,
                  EV1 = pca.snpset2$eigenvect[,1],    # the first eigenvector
                  EV2 = pca.snpset2$eigenvect[,2],    # the second eigenvector
                  EV3 = pca.snpset2$eigenvect[,3],
                  EV4 = pca.snpset2$eigenvect[,4],
                  stringsAsFactors = FALSE)
#tab <- read.table("Msa.pca.snpset2$eigenvect.txt", header=T)
head(tab)

pc.percent <- pca.snpset2$varprop*100
head(round(pc.percent, 2), 10)

sum(head(pc.percent,10)) #17.7109

phenames=c("Ccirc.year3", "Bcirc.year3", "CmL.year3", "DBI.year3", "Yld.year3")
num.vec <- c(50,100,150,200,250,300,350)
col.vec=c("#999999", "#E69F00", "#56B4E9","#D55E00", "#0072B2", "#CC79A7", "#009E73","#9932CC")

    for(i in 1:length(phenames)){
      
      pdf(paste("Msa.SNPrelate", phenames[i], "pdf", sep="."), width=8, height=16)
      par(mfrow=c(4,2))
      
      for(j in 1:length(num.vec)){
        
        CDmean.names.Msa <- read.delim(paste("Msa.CDmeanMax1.itrn3000.individuals", phenames[i],"size-n:",num.vec[j],"txt", sep="."), header=F)
              tab.vec <- tab[which(CDmean.names.Msa$V1%in%tab$sample.id),]
              par(mar=c(4,5,4,4))
              plot(tab$EV1, tab$EV2, xlab="PC1 = 5.4%", ylab="PC2 = 3.4%", pch=21, cex=1.6, 
                   cex.lab=1.5, cex.axis=1.5, main=num.vec[j], col="black")
              points(tab.vec$EV1, tab.vec$EV2, pch=21, cex=1.4, col=adjustcolor("black", alpha=0.8), bg=adjustcolor(col.vec[j], alpha=0.5))
              CDmean.names.Msa=NULL
              tab.vec=NULL
      }
      dev.off()
    }

#####DAPC Groups
Msa.DAPC <- read.csv("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/09F2/Msa/Msa_DAPC.csv", header=T)
colnames(Msa.DAPC)[1] <- "sample.id"
tab.DAPC <- merge(tab, Msa.DAPC, by="sample.id")
DAPC.vec <- as.character(unique(tab.DAPC$DAPC))

pdf("Msa.DAPCgrps.on.PCA.pdf", 7,6)
par(mar=c(4,6,4,4))
    plot(tab$EV1, tab$EV2, xlab="PC1 = 5.4%", ylab="PC2 = 3.4%", pch=21, cex=1.6, 
         cex.lab=1.5, cex.axis=1.5, main="DAPC Groups", col="black")
    for(k in 1:length(DAPC.vec)){
      tab.DAPC.grp <- tab.DAPC[which(tab.DAPC$DAPC==DAPC.vec[k]),]
    points(tab.DAPC.grp$EV1, tab.DAPC.grp$EV2, pch=21, cex=1.4, col=adjustcolor("black", alpha=0.8), bg=adjustcolor(col.vec[k], alpha=0.5))
    tab.DAPC.grp=NULL
    
    }
    legend("topright",legend=DAPC.vec,
           col=col.vec, pch=20, bty="n",cex=1, pt.cex =1.5)
    
dev.off()



rm(list=ls())
library(gdsfmt)
library(SNPRelate)

setwd("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/09F2/Msi")
# Msi
vcf.fn <- "/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/09F2/Msi/Msi.Geno.5pmaf.vcf"
snpgdsVCF2GDS(vcf.fn, "pm.gds", method="copy.num.of.ref")
snpgdsSummary("pm.gds")

# open the GDS file
genofile <- openfn.gds("pm.gds")
## LD based pruning of SNPs
set.seed(1000)
snpset2 <- snpgdsLDpruning(genofile, method = c("corr"), ld.threshold = 0.5,missing = 0.8, maf = 0.02)
snpset2.id <- unlist(snpset2)

########## PCA Analysis ###############
# get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#PCA analysis

#pm_pca<-snpgdsPCA(genofile,autosome.only=F)
pca.snpset2 <- snpgdsPCA(genofile, snp.id=snpset2.id, num.thread=2)

write.table(pca.snpset2$eigenvect, "Msi.pca.snpset2$eigenvect.txt")
write.table(pca.snpset2$eigenval, "Msi.pca.snpset2$eigenval.txt")
write.table(pca.snpset2$sample.id,"Msi.pca.snpset2$sample.id.txt")

tab <- data.frame(sample.id = pca.snpset2$sample.id,
                  EV1 = pca.snpset2$eigenvect[,1],    # the first eigenvector
                  EV2 = pca.snpset2$eigenvect[,2],    # the second eigenvector
                  EV3 = pca.snpset2$eigenvect[,3],
                  EV4 = pca.snpset2$eigenvect[,4],
                  stringsAsFactors = FALSE)
#tab <- read.table("Msi.pca.snpset2$eigenvect.txt", header=T)
head(tab)

pc.percent <- pca.snpset2$varprop*100
head(round(pc.percent, 2), 10)

sum(head(pc.percent,10)) #17.7109

phenames=c("Ccirc.year3", "Bcirc.year3", "CmL.year3", "DBI.year3", "Yld.year3")
num.vec <- c(50,100,150,200,250,300,350)
col.vec=c("#999999", "#E69F00", "#56B4E9","#D55E00", "#0072B2", "#CC79A7", "#009E73","#9932CC")

for(i in 1:length(phenames)){
  
  pdf(paste("Msi.SNPrelate", phenames[i], "pdf", sep="."), width=8, height=16)
  par(mfrow=c(4,2))
  
  for(j in 1:length(num.vec)){
    
    CDmean.names.Msi <- read.delim(paste("Msi.CDmeanMax1.itrn3000.individuals", phenames[i],"size-n:",num.vec[j],"txt", sep="."), header=F)
    tab.vec <- tab[which(CDmean.names.Msi$V1%in%tab$sample.id),]
    par(mar=c(4,5,4,4))
    plot(tab$EV1, tab$EV2, xlab="PC1 = 5.4%", ylab="PC2 = 3.4%", pch=21, cex=1.6, 
         cex.lab=1.5, cex.axis=1.5, main=num.vec[j], col="black")
    points(tab.vec$EV1, tab.vec$EV2, pch=21, cex=1.4, col=adjustcolor("black", alpha=0.8), bg=adjustcolor(col.vec[j], alpha=0.5))
    CDmean.names.Msi=NULL
    tab.vec=NULL
  }
  dev.off()
}

#####DAPC Groups
Msi.DAPC <- read.csv("/homes/omo/public_html/effect_simulations/PostDoc/Manuscript2/Data/09F2/Msi/Msi_DAPC.csv", header=T)
colnames(Msi.DAPC)[1] <- "sample.id"
tab.DAPC <- merge(tab, Msi.DAPC, by="sample.id")
DAPC.vec <- as.character(unique(tab.DAPC$DAPC))

pdf("Msi.DAPCgrps.on.PCA.pdf", 7,6)
par(mar=c(4,6,4,4))
plot(tab$EV1, tab$EV2, xlab="PC1 = 5.4%", ylab="PC2 = 3.4%", pch=21, cex=1.6, 
     cex.lab=1.5, cex.axis=1.5, main="DAPC Groups", col="black")
for(k in 1:length(DAPC.vec)){
  tab.DAPC.grp <- tab.DAPC[which(tab.DAPC$DAPC==DAPC.vec[k]),]
  points(tab.DAPC.grp$EV1, tab.DAPC.grp$EV2, pch=21, cex=1.4, col=adjustcolor("black", alpha=0.8), bg=adjustcolor(col.vec[k], alpha=0.5))
  tab.DAPC.grp=NULL
  
}
legend("topright",legend=DAPC.vec,
       col=col.vec, pch=20, bty="n",cex=1, pt.cex =1.5)

dev.off()













