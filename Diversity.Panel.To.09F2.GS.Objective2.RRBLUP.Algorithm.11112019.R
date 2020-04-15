
DiversityPanel.2.BreedingPool.GS.Pipeline <- function(phenames=NULL, Msi.Msa.Geno=NULL, Msi.Msa.Pheno=NULL, valid.geno=NULL, f2.vec=NULL, bootstraps.vec=NULL,
                                                      valid.pheno=NULL, Msi.div.taxa=NULL, Msa.div.taxa=NULL, her.vec=NULL, PCA.info.Msi=NULL, PCA.info.Msa=NULL, Opt.PCA=NULL,
                                                      dat.msi=NULL, dat.msa=NULL, num.vec=NULL, s=NULL, Iter=NULL, burn=NULL, species=NULL, k=NULL, itrn=NULL, 
                                                      trn.app=NULL, number.of.folds=NULL, Msi.Div.RawPheno=NULL, Msa.Div.RawPheno=NULL, date=NULL){
  
  
  r.gy.output.ADE.Traits <- NULL
  the.coefficients.output.ADE.int.Traits <- NULL
  the.coefficients.output.ADE.slope.Traits <- NULL
  myCI.output.ADE.Traits <- NULL
  myGSR.output.ADE.Traits <- NULL
  
  for(l in 1:length(phenames)){ #1:length(phenames)
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], "!!!!!!!----------", sep = ""))
    
    Msi.Msa.Pheno.trt=as.matrix(Msi.Msa.Pheno[,phenames[l]])
    rownames(Msi.Msa.Pheno.trt) <- Msi.Msa.Pheno$Taxa
    colnames(Msi.Msa.Pheno.trt) <- phenames[l]
    head(Msi.Msa.Pheno.trt)
    
    Pheno.Msi.DivPanel= as.matrix(Msi.Msa.Pheno.trt[match(Msi.div.taxa, rownames(Msi.Msa.Pheno.trt)),])
    colnames(Pheno.Msi.DivPanel) <- phenames[l]
    
    Pheno.Msa.DivPanel= as.matrix(Msi.Msa.Pheno.trt[match(Msa.div.taxa, rownames(Msi.Msa.Pheno.trt)),])
    colnames(Pheno.Msa.DivPanel) <- phenames[l]
    
    Pheno.MsiMsa.DivPanel = Msi.Msa.Pheno.trt
    
    pheno_testset=as.matrix(valid.pheno[,phenames[l]])
    rownames(pheno_testset) <- valid.pheno$Taxa
    colnames(pheno_testset) <- phenames[l]
    head(pheno_testset)
    
    G <- Msi.Msa.Geno
    
    herit <- her.vec[which(her.vec$trait==phenames[l]),1]

    r.gy.ADE <- matrix(NA, length(trn.app), bootstraps.vec)
    myCI.ADE <- matrix(NA, length(trn.app), bootstraps.vec)
    myGSR.ADE <- matrix(NA, length(trn.app), bootstraps.vec)
    the.coefficients.ADE.int <- matrix(NA, length(trn.app), bootstraps.vec)
    the.coefficients.ADE.slope <- matrix(NA, length(trn.app), bootstraps.vec)
    
    r.gy.PCA<- matrix(NA, length(trn.app), bootstraps.vec)
    myCI.PCA <- matrix(NA, length(trn.app), bootstraps.vec)
    myGSR.PCA <- matrix(NA, length(trn.app), bootstraps.vec)
    the.coefficients.PCA.int <- matrix(NA, length(trn.app), bootstraps.vec)
    the.coefficients.PCA.slope <- matrix(NA, length(trn.app), bootstraps.vec)
    
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- CrossFolds "," !!!!!!!----------", sep = ""))
   
      Geno_testset <- valid.geno
      Geno_testset <- as.matrix(Geno_testset)
      
      # Apply CD Mean to select optimal training set
      Dat_train.Msi <- dat.msi-1
      
      A.A.Msi <- A.mat(Dat_train.Msi)
      matA1.Msi <- A.mat(Dat_train.Msi)
      matA1.Msi=as.matrix(matA1.Msi)
      
      # Estimate variance components
      y.Msi<-Msi.Div.RawPheno[,phenames[l]]
      y.Msi[is.na(y.Msi)]<-mean(y.Msi,na.rm=T)
      model.Msi <- mixed.solve(y.Msi,K=A.A.Msi)
      
      varG.Msi = model.Msi$Vu
      varE.Msi <- model.Msi$Ve
      
      print(paste("--------- Initiating CDMean for Msi .....in trait : ", phenames[l], " --------", sep = ""))
      
      # CDmean function to return optimal 200 individuals
      sel.train.CDmean.Msi <- CDMean(matA1 = matA1.Msi, varG = varG.Msi, varE = varE.Msi, nindrep=num.vec, itrn = itrn)
      write.table(sel.train.CDmean.Msi, paste("Msi", "CDmean.Selected.Individuals.200", "Trait: ", phenames[l], "txt", sep="."), sep="\t", quote=F, row.names=F, col.names = F)
      
      #sel.train.CDmean.Msi.df <- data.frame(sel.train.CDmean.Msi)
      #colnames(sel.train.CDmean.Msi.df) <- "Taxa"
      #Subset phenotypic information for Msi CDmean combination
      sel.train.CDmean.Msi.pheno <- as.matrix(Msi.Msa.Pheno.trt[match(sel.train.CDmean.Msi, rownames(Msi.Msa.Pheno.trt)),])
      colnames(sel.train.CDmean.Msi.pheno) <- phenames[l]
      
      random.samp.msi <- sample(nrow(Pheno.Msi.DivPanel), num.vec, replace = FALSE)
      random_train.Msi=Pheno.Msi.DivPanel[random.samp.msi,]
      random_train.Msi <- as.data.frame(random_train.Msi)
      colnames(random_train.Msi) <- colnames(Pheno.Msi.DivPanel)
      random_train.Msi <- as.matrix(random_train.Msi)
      #### Msa
      Dat_train.Msa <- dat.msa-1
      
      A.A.Msa <- A.mat(Dat_train.Msa)
      matA1.Msa <- A.mat(Dat_train.Msa)
      matA1.Msa=as.matrix(matA1.Msa)
      
      # Estimate variance components
      y.Msa<-Msa.Div.RawPheno[,phenames[l]]
      y.Msa[is.na(y.Msa)]<-mean(y.Msa,na.rm=T)
      model.Msa <- mixed.solve(y.Msa,K=A.A.Msa)
      
      varG.Msa = model.Msa$Vu
      varE.Msa <- model.Msa$Ve
      
      print(paste("--------- Initiating CDMean for Msa .....in trait : ", phenames[l], " --------", sep = ""))
      
      # CDmean function to return optimal 200 individuals
      sel.train.CDmean.Msa <- CDMean(matA1 = matA1.Msa, varG = varG.Msa, varE = varE.Msa, nindrep=num.vec, itrn = itrn)
      write.table(sel.train.CDmean.Msa, paste("Msa", "CDmean.Selected.Individuals.200", "Trait: ", phenames[l], "txt", sep="."), sep="\t", quote=F, row.names=F, col.names = F)
      
      #sel.train.CDmean.Msa.df <- data.frame(sel.train.CDmean.Msa)
      #colnames(sel.train.CDmean.Msa.df) <- "Taxa"
      
      #Subset phenotypic information for Msa CDmean combination
      sel.train.CDmean.Msa.pheno <- as.matrix(Msi.Msa.Pheno.trt[match(sel.train.CDmean.Msa, rownames(Msi.Msa.Pheno.trt)),])
      colnames(sel.train.CDmean.Msa.pheno) <- phenames[l]
      
      #Subset phenotypic information for MsiMsa CDmean combination
      sel.train.CDmean.Msi.Msa <- c(sel.train.CDmean.Msi, sel.train.CDmean.Msa)
      sel.train.CDmean.Msi.Msa.pheno <- as.matrix(Msi.Msa.Pheno.trt[match(sel.train.CDmean.Msi.Msa, rownames(Msi.Msa.Pheno.trt)),])
      colnames(sel.train.CDmean.Msi.Msa.pheno) <- phenames[l]
      
      random.samp.msa <- sample(nrow(Pheno.Msa.DivPanel), num.vec, replace = FALSE)
      random_train.Msa=Pheno.Msa.DivPanel[random.samp.msa,]
      random_train.Msa <- as.data.frame(random_train.Msa)
      colnames(random_train.Msa) <- colnames(Pheno.Msa.DivPanel)
      random_train.Msa <- as.matrix(random_train.Msa)
      
      
      # Loop for calling each of the different approaches Raw phenotypes, PCA, CDmean, DAPC
      for(trn in 1:length(trn.app)){
        
        if(trn.app[trn]=="Msi.Random"){ pheno_trainingset = random_train.Msi
        G.PCA <- PCA.info.Msi
        G2.PCA <- G.PCA[,-1]
        rownames(G2.PCA) <- G.PCA$Taxa
        } else {
          if(trn.app[trn]=="Msa.Random"){ pheno_trainingset= random_train.Msa
          G.PCA <- PCA.info.Msa
          G2.PCA <- G.PCA[,-1]
          rownames(G2.PCA) <- G.PCA$Taxa
          } else { 
            if(trn.app[trn]=="Msi.Div"){ pheno_trainingset=Pheno.Msi.DivPanel
            G.PCA <- PCA.info.Msi
            G2.PCA <- G.PCA[,-1]
            rownames(G2.PCA) <- G.PCA$Taxa
            } else { 
              if(trn.app[trn]=="Msa.Div"){ pheno_trainingset = Pheno.Msa.DivPanel
              G.PCA <- PCA.info.Msa
              G2.PCA <- G.PCA[,-1]
              rownames(G2.PCA) <- G.PCA$Taxa
            }else{
              if(trn.app[trn]=="Msi.Div.CDMean"){ pheno_trainingset=sel.train.CDmean.Msi.pheno
              G.PCA <- PCA.info.Msi
              G2.PCA <- G.PCA[,-1]
              rownames(G2.PCA) <- G.PCA$Taxa
            }else{
                  if(trn.app[trn]=="Msa.Div.CDMean"){pheno_trainingset = sel.train.CDmean.Msa.pheno
                  G.PCA <- PCA.info.Msa
                  G2.PCA <- G.PCA[,-1]
                  rownames(G2.PCA) <- G.PCA$Taxa
                   }else{
                        pheno_trainingset.1 = Pheno.Msi.DivPanel
                        pheno_trainingset.2 = Pheno.Msa.DivPanel
                  }
                  #rt.nam <- rownames(pheno_testset); pheno_testset <- pheno_testset[,-2]; names(pheno_testset) <- rt.nam}
                }
                
              }   # 0% VAT 
            }
          }
        
        }
        
          
        print(paste("--- Analyzing ... Trait : ", phenames[l], " --- CrossFolds ", "Training Approach : ",trn.app[trn]," !!!---", sep = ""))
        #####################################
        
        if(trn.app[trn]!="msi.plus.msa"){
          
        head(pheno_trainingset)
        head(pheno_testset)
        
        
        Yt=(pheno_trainingset[,1])
        
        Geno_trainingset <- G[match(rownames(pheno_trainingset), rownames(G)),]
        
        ww <- pheno_testset

        Geno_trainingset.PCA <- G2.PCA[match(rownames(pheno_trainingset), rownames(G2.PCA)),]
        Geno_trainingset.PCA <- as.matrix(Geno_trainingset.PCA)
        
        #Geno_testset.PCA <- G2.PCA[match(rownames(pheno_testset), rownames(G2.PCA)),]
        #Geno_testset.PCA <- as.matrix(Geno_testset.PCA)
        
        rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=Geno_trainingset, K=NULL, SE=F, return.Hinv=F)
        mEff.Ro<-rrMod.Ro$u
        e.Ro= as.matrix(mEff.Ro)
        
        for(b in 1:bootstraps.vec){
          
          print(paste("--- Analyzing ... Trait : ", phenames[l], " --- CrossFolds ", "Training Approach : ",trn.app[trn]," GS RRBLUP ; Bootsrap : ", b," out of ", bootstraps.vec," !!!---", sep = ""))
                  
                
                  Geno_testset <- Geno_testset[match(rownames(pheno_testset), rownames(Geno_testset)),]
                  
                  sam.boot <- sample(nrow(Geno_testset), f2.vec, replace = TRUE)
                  Geno_testset.boot <- Geno_testset[sam.boot,]
                  pheno_testset.boot <- pheno_testset[sam.boot,]
                  predYv.Ro = Geno_testset.boot%*%e.Ro
                  
                  predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
                  
                  Y_valid=pheno_testset.boot
                  #cor.predacc <- cor(predYr.Ro,Y_valid,use="complete")
                  #print(cor.predacc)
      
                  r.gy.ADE[trn,b] <-  cor(predYr.Ro,Y_valid,use="complete")
                  
                  the.fitted.model.ADE <- lm(Y_valid ~ predYr.Ro)
                  the.coefficients.ADE.int[trn,b] <- the.fitted.model.ADE$coefficients[1]
                  the.coefficients.ADE.slope[trn,b] <- the.fitted.model.ADE$coefficients[2]
                  
                  these.observed.and.predicted.phenotypic.values.Ade <- data.frame(rownames(ww), Y_valid, predYr.Ro)
                  rownames(these.observed.and.predicted.phenotypic.values.Ade) <- NULL
                  colnames(these.observed.and.predicted.phenotypic.values.Ade) <- c("Taxa", "Observed.Value", "Predicted.Value")
                  x.p.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,3)]
                  y.o.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,2)]
                  #write.table(these.observed.and.predicted.phenotypic.values.Ade, paste("these.observed.and.predicted.phenotypic.values.RRBLUP", phenames[l], trn.app[trn],"txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
                  myCI.ADE[trn,b] <- round(CI(x.p.Ade,y.o.Ade,s=s,top=T),2)
                  
                  these.observed.and.predicted.phenotypic.values.Ade <- these.observed.and.predicted.phenotypic.values.Ade[order(these.observed.and.predicted.phenotypic.values.Ade$Predicted.Value, decreasing = T),]
                  #sel.taxa.Ade <- these.observed.and.predicted.phenotypic.values.Ade$Taxa[c(1:(round(length(these.observed.and.predicted.phenotypic.values.Ade$Taxa)*0.1,0)))]
                  #sel.taxa.Ade <- as.character(sel.taxa.Ade)
                  #sel.taxa.Ade.numvec <- c(sel.taxa.Ade.numvec, sel.taxa.Ade)
                  corr.sqr.Ade <- round(summary(the.fitted.model.ADE)$adj.r.squared, 2)
                  myGSR.ADE[trn,b] <- round((corr.sqr.Ade/herit), 2)
        
              }# End Bootstrap
        
        } else{
          
          Yt.1=(pheno_trainingset.1[,1])
          Yt.2=(pheno_trainingset.2[,1])
          
          Geno_trainingset.1 <- G[match(rownames(pheno_trainingset.1), rownames(G)),]
          
          Geno_trainingset.2 <- G[match(rownames(pheno_trainingset.2), rownames(G)),]
          
          ww <- pheno_testset
          
          rrMod.Ro.1 <-mixed.solve(Yt.1, X=NULL, Z=Geno_trainingset.1, K=NULL, SE=F, return.Hinv=F)
          mEff.Ro.1<-rrMod.Ro.1$u
          e.Ro.1= as.matrix(mEff.Ro.1)
          
          rrMod.Ro.2<-mixed.solve(Yt.2, X=NULL, Z=Geno_trainingset.2, K=NULL, SE=F, return.Hinv=F)
          mEff.Ro.2 <- rrMod.Ro.2$u
          e.Ro.2 = as.matrix(mEff.Ro.2)
          
          for(b in 1:bootstraps.vec){
            
            print(paste("--- Analyzing ... Trait : ", phenames[l], " --- CrossFolds ", "Training Approach : ",trn.app[trn]," GS RRBLUP ; Bootsrap : ", b," out of ", bootstraps.vec," !!!---", sep = ""))
            
            
            Geno_testset <- Geno_testset[match(rownames(pheno_testset), rownames(Geno_testset)),]
            
            sam.boot <- sample(nrow(Geno_testset), f2.vec, replace = TRUE)
            Geno_testset.boot <- Geno_testset[sam.boot,]
            pheno_testset.boot <- pheno_testset[sam.boot,]
            
            predYv.Ro.1 = Geno_testset.boot%*%e.Ro.1
            predYr.Ro.1 = predYv.Ro.1[,1]+ rrMod.Ro.1$beta
            
            predYv.Ro.2 = Geno_testset.boot%*%e.Ro.2
            predYr.Ro.2 = predYv.Ro.2[,1]+ rrMod.Ro.2$beta
            
            predYr.Ro <- predYr.Ro.1+predYr.Ro.2
            Y_valid=pheno_testset.boot
            
            r.gy.ADE[trn,b] <-  cor(predYr.Ro,Y_valid,use="complete")
            
            the.fitted.model.ADE <- lm(Y_valid ~ predYr.Ro)
            the.coefficients.ADE.int[trn,b] <- the.fitted.model.ADE$coefficients[1]
            the.coefficients.ADE.slope[trn,b] <- the.fitted.model.ADE$coefficients[2]
            
            these.observed.and.predicted.phenotypic.values.Ade <- data.frame(rownames(ww), Y_valid, predYr.Ro)
            rownames(these.observed.and.predicted.phenotypic.values.Ade) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Ade) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,3)]
            y.o.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,2)]
            #write.table(these.observed.and.predicted.phenotypic.values.Ade, paste("these.observed.and.predicted.phenotypic.values.RRBLUP", phenames[l], trn.app[trn],"txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
            myCI.ADE[trn,b] <- round(CI(x.p.Ade,y.o.Ade,s=s,top=T),2)
            
            these.observed.and.predicted.phenotypic.values.Ade <- these.observed.and.predicted.phenotypic.values.Ade[order(these.observed.and.predicted.phenotypic.values.Ade$Predicted.Value, decreasing = T),]
            #sel.taxa.Ade <- these.observed.and.predicted.phenotypic.values.Ade$Taxa[c(1:(round(length(these.observed.and.predicted.phenotypic.values.Ade$Taxa)*0.1,0)))]
            #sel.taxa.Ade <- as.character(sel.taxa.Ade)
            #sel.taxa.Ade.numvec <- c(sel.taxa.Ade.numvec, sel.taxa.Ade)
            corr.sqr.Ade <- round(summary(the.fitted.model.ADE)$adj.r.squared, 2)
            myGSR.ADE[trn,b] <- round((corr.sqr.Ade/herit), 2)
            
          }# End Bootstrap
          
          }
        
        
      }# Raw phenotype, end of CDmean, end of PCA, end of DAPC
      
      pheno_trainingset=NULL
      pheno_testset=NULL
      Geno_trainingset=NULL
      Geno_testset=NULL
      
      
    #} # End of i in k cross valid
      
      
      
      r.gy.output.ADE <- as.data.frame(r.gy.ADE)
      r.gy.output.ADE$Trait <- rep(phenames[l], nrow(r.gy.output.ADE))
      r.gy.output.ADE$Approach <- as.vector(trn.app)
      
      the.coefficients.output.ADE.int <- as.data.frame(the.coefficients.ADE.int)
      the.coefficients.output.ADE.int$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.int))
      the.coefficients.output.ADE.int$Approach <- trn.app
      
      the.coefficients.output.ADE.slope <- as.data.frame(the.coefficients.ADE.slope)
      the.coefficients.output.ADE.slope$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.slope))
      the.coefficients.output.ADE.slope$Approach <- trn.app
      
      myCI.output.ADE<- as.data.frame(myCI.ADE)
      myCI.output.ADE$Trait <- rep(phenames[l], nrow(myCI.output.ADE))
      myCI.output.ADE$Approach <- trn.app
      
      myGSR.output.ADE <- as.data.frame(myGSR.ADE)
      myGSR.output.ADE$Trait <- rep(phenames[l], nrow(myGSR.output.ADE))
      myGSR.output.ADE$Approach <- trn.app
      
      write.table(r.gy.output.ADE, paste(species,phenames[l], "r.gy.output.RRBLUP.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      write.table(the.coefficients.output.ADE.int, paste(species,phenames[l], "the.coefficients.output.RRBLUP.int.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      write.table(the.coefficients.output.ADE.slope, paste(species,phenames[l], "the.coefficients.output.RRBLUP.slope.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      write.table(myCI.output.ADE, paste(species,phenames[l], "myCI.output.RRBLUP.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      write.table(myGSR.output.ADE, paste(species,phenames[l], "myGSR.output.RRBLUP.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      
      r.gy.output.ADE.Traits <- rbind(r.gy.output.ADE.Traits, r.gy.output.ADE)
      the.coefficients.output.ADE.int.Traits <- rbind(the.coefficients.output.ADE.int.Traits, the.coefficients.output.ADE.int)
      the.coefficients.output.ADE.slope.Traits <- rbind(the.coefficients.output.ADE.slope.Traits, the.coefficients.output.ADE.slope)
      myCI.output.ADE.Traits <- rbind(myCI.output.ADE.Traits, myCI.output.ADE)
      myGSR.output.ADE.Traits <- rbind(myGSR.output.ADE.Traits, myGSR.output.ADE)
    
    
  } # End of trait
  
  print("--------- Writing Final Results In Working Directory ------------")
  
  save.image(paste(species, "GS.GS.Results.Div.5FCV.TRSvec.RRBLUP.rda", sep='.'))
  
  write.table(r.gy.output.ADE.Traits, paste(species, "r.gy.output.RRBLUP.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.int.Traits, paste(species, "the.coefficients.output.RRBLUP.int.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.slope.Traits, paste(species, "the.coefficients.output.RRBLUP.slope.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(myCI.output.ADE.Traits, paste(species, "myCI.output.RRBLUP.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(myGSR.output.ADE.Traits, paste(species, "myGSR.output.RRBLUP.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)

  } #
#} # End of library
# Started running 07022019 09:17 AM


