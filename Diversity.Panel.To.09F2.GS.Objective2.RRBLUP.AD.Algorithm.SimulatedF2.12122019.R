
DiversityPanel.2.BreedingPool.GS.Pipeline <- function(phenames=NULL, Msi.polyRAD.geno=NULL, Msa.polyRAD.geno=NULL, 
                                                      Msi.Div.RawPheno=NULL, Msa.Div.RawPheno=NULL, 
                                                      valid.geno=NULL, valid.pheno=NULL, f2.vec=NULL,
                                                      Msi.div.taxa=NULL, Msa.div.taxa=NULL, 
                                                      PCA.info.Msi=NULL, PCA.info.Msa=NULL,
                                                      num.vec=NULL, s=NULL, species=NULL, k=NULL, itrn=NULL, 
                                                      trn.app=NULL, number.of.folds=NULL, date=NULL){
  
  
  r.gy.output.ADE.Traits <- NULL
  the.coefficients.output.ADE.int.Traits <- NULL
  the.coefficients.output.ADE.slope.Traits <- NULL
  myCI.output.ADE.Traits <- NULL
  myGSR.output.ADE.Traits <- NULL
  
  for(l in 1:length(phenames)){ #1:length(phenames)
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], "!!!!!!!----------", sep = ""))
    
    Msi.Pheno.trt=as.matrix(Msi.Div.RawPheno[,phenames[l]])
    rownames(Msi.Pheno.trt) <- Msi.Div.RawPheno$Taxa
    colnames(Msi.Pheno.trt) <- phenames[l]
    
    Msa.Pheno.trt=as.matrix(Msa.Div.RawPheno[,phenames[l]])
    rownames(Msa.Pheno.trt) <- Msa.Div.RawPheno$Taxa
    colnames(Msa.Pheno.trt) <- phenames[l]
    
    Pheno.Msi.DivPanel= Msi.Pheno.trt
    
    Pheno.Msa.DivPanel= Msa.Pheno.trt
    
    pheno_testset=as.matrix(valid.pheno[,phenames[l]])
    rownames(pheno_testset) <- valid.pheno$Taxa
    colnames(pheno_testset) <- phenames[l]
    head(pheno_testset)
    
    #G.msi <- 
    

    r.gy.ADE <- matrix(NA, length(trn.app), 1)
    myCI.ADE <- matrix(NA, length(trn.app), 1)
    myGSR.ADE <- matrix(NA, length(trn.app), 1)
    the.coefficients.ADE.int <- matrix(NA, length(trn.app), 1)
    the.coefficients.ADE.slope <- matrix(NA, length(trn.app), 1)
    
    r.gy.AD <- matrix(NA, length(trn.app), 1)
    myCI.AD <- matrix(NA, length(trn.app), 1)
    myGSR.AD <- matrix(NA, length(trn.app), 1)
    the.coefficients.AD.int <- matrix(NA, length(trn.app), 1)
    the.coefficients.AD.slope <- matrix(NA, length(trn.app), 1)
    
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- CrossFolds "," !!!!!!!----------", sep = ""))
   
      Geno_testset <- valid.geno
      Geno_testset <- as.matrix(Geno_testset)
      
      # Apply CD Mean to select optimal training set
      Dat_train.Msi <- Msi.polyRAD.geno-1
      
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
      
      sel.train.CDmean.Msi.pheno <- as.matrix(Msi.Pheno.trt[match(sel.train.CDmean.Msi, rownames(Msi.Pheno.trt)),])
      colnames(sel.train.CDmean.Msi.pheno) <- phenames[l]
      
      random.samp.msi <- sample(nrow(Pheno.Msi.DivPanel), num.vec, replace = FALSE)
      random_train.Msi=Pheno.Msi.DivPanel[random.samp.msi,]
      random_train.Msi <- as.data.frame(random_train.Msi)
      colnames(random_train.Msi) <- colnames(Pheno.Msi.DivPanel)
      random_train.Msi <- as.matrix(random_train.Msi)
      #### Msa
      Dat_train.Msa <- Msa.polyRAD.geno-1
      
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

      #Subset phenotypic information for Msa CDmean combination
      sel.train.CDmean.Msa.pheno <- as.matrix(Msa.Pheno.trt[match(sel.train.CDmean.Msa, rownames(Msa.Pheno.trt)),])
      colnames(sel.train.CDmean.Msa.pheno) <- phenames[l]
      
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
        G <- Msi.polyRAD.geno
        } else {
          if(trn.app[trn]=="Msa.Random"){ pheno_trainingset= random_train.Msa
          G.PCA <- PCA.info.Msa
          G2.PCA <- G.PCA[,-1]
          rownames(G2.PCA) <- G.PCA$Taxa
          G <- Msa.polyRAD.geno
          } else { 
            if(trn.app[trn]=="Msi.Div"){ pheno_trainingset=Pheno.Msi.DivPanel
            G.PCA <- PCA.info.Msi
            G2.PCA <- G.PCA[,-1]
            rownames(G2.PCA) <- G.PCA$Taxa
            G <- Msi.polyRAD.geno
            } else { 
              if(trn.app[trn]=="Msa.Div"){ pheno_trainingset = Pheno.Msa.DivPanel
              G.PCA <- PCA.info.Msa
              G2.PCA <- G.PCA[,-1]
              rownames(G2.PCA) <- G.PCA$Taxa
              G <- Msa.polyRAD.geno
            }else{
              if(trn.app[trn]=="Msi.Div.CDMean"){ pheno_trainingset=sel.train.CDmean.Msi.pheno
              G.PCA <- PCA.info.Msi
              G2.PCA <- G.PCA[,-1]
              rownames(G2.PCA) <- G.PCA$Taxa
              G <- Msi.polyRAD.geno
            }else{
                  if(trn.app[trn]=="Msa.Div.CDMean"){pheno_trainingset = sel.train.CDmean.Msa.pheno
                  G.PCA <- PCA.info.Msa
                  G2.PCA <- G.PCA[,-1]
                  rownames(G2.PCA) <- G.PCA$Taxa
                  G <- Msa.polyRAD.geno
                   }else{
                        pheno_trainingset.1 = Pheno.Msi.DivPanel
                        G.1 <- Msi.polyRAD.geno
                        pheno_trainingset.2 = Pheno.Msa.DivPanel
                        G.2 <- Msa.polyRAD.geno
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
        
        library(rrBLUP)
        rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=Geno_trainingset, K=NULL, SE=F, return.Hinv=F)
        mEff.Ro<-rrMod.Ro$u
        e.Ro= as.matrix(mEff.Ro)
                  Geno_testset <- Geno_testset[match(rownames(pheno_testset), rownames(Geno_testset)),]
                  
                  predYv.Ro = Geno_testset%*%e.Ro
                  
                  predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
                  
                  Y_valid=pheno_testset
                  r.gy.ADE[trn,1] <-  cor(predYr.Ro,Y_valid,use="complete")
                  
                  the.fitted.model.ADE <- lm(Y_valid ~ predYr.Ro)
                  the.coefficients.ADE.int[trn,1] <- the.fitted.model.ADE$coefficients[1]
                  the.coefficients.ADE.slope[trn,1] <- the.fitted.model.ADE$coefficients[2]
                  
                  these.observed.and.predicted.phenotypic.values.Ade <- data.frame(rownames(ww), Y_valid, predYr.Ro)
                  rownames(these.observed.and.predicted.phenotypic.values.Ade) <- NULL
                  colnames(these.observed.and.predicted.phenotypic.values.Ade) <- c("Taxa", "Observed.Value", "Predicted.Value")
                  x.p.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,3)]
                  y.o.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,2)]
                  #write.table(these.observed.and.predicted.phenotypic.values.Ade, paste("these.observed.and.predicted.phenotypic.values.RRBLUP", phenames[l], trn.app[trn],"txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
                  myCI.ADE[trn,1] <- round(CI(x.p.Ade,y.o.Ade,s=s,top=T),2)
                  
                  these.observed.and.predicted.phenotypic.values.Ade <- these.observed.and.predicted.phenotypic.values.Ade[order(these.observed.and.predicted.phenotypic.values.Ade$Predicted.Value, decreasing = T),]
                  #corr.sqr.Ade <- round(summary(the.fitted.model.ADE)$adj.r.squared, 2)
                  #myGSR.ADE[trn,b] <- round((corr.sqr.Ade/herit), 2)
                  
                  
                  
                  library(sommer)
                  
                  Pheno.comb <- rbind(pheno_trainingset, pheno_testset)
                  train.comb <- 1:nrow(pheno_trainingset)
                  valid.comb <-setdiff(1:nrow(Pheno.comb),train.comb)
                  
                  y.trn <- Pheno.comb # for prediction accuracy
                  ww <- pheno_testset # delete data for 1/5 of the population
                  
                  y.trn[valid.comb,1] <- NA #valid.comb
                  
                  Geno.comb <- rbind(Geno_trainingset, Geno_testset)
                  Geno.comb <- as.matrix(Geno.comb)
                  A1 <- A.mat(Geno.comb,shrink=TRUE) # additive relationship matrix 
                  # D.mat also take -1,0,and 1
                  D1 <- D.mat(Geno.comb,shrink=TRUE)
                  
                  Za <- diag(length(Pheno.comb[,1])) 
                  Zd <- diag(length(Pheno.comb[,1]))
                  rownames(A1) <- 1:nrow(A1)
                  rownames(D1) <- 1:nrow(D1)
                  
                  Y_valid=pheno_testset

                  ETA.AD <- list(add=list(Z=Za,K=A1), dom=list(Z=Zd,K=D1))
                  ans.AD <- MEMMA(Y=y.trn, ZETA=ETA.AD) 
                  r.gy.AD[trn,1] <- cor(ans.AD$fitted.y[valid.comb], ww, use="complete")
                  
                  the.fitted.model.AD <- lm(ww ~ ans.AD$fitted.y[valid.comb])
                  the.coefficients.AD <- c(the.fitted.model.AD$coefficients[1], the.fitted.model.AD$coefficients[2])
                  the.coefficients.AD.int[trn,1] <- the.coefficients.AD[1]
                  the.coefficients.AD.slope[trn,1] <- the.coefficients.AD[2]
                  
                  these.observed.and.predicted.phenotypic.values.Ad <- data.frame(rownames(Y_valid), Y_valid, ans.AD$fitted.y[valid.comb])
                  rownames(these.observed.and.predicted.phenotypic.values.Ad) <- NULL
                  colnames(these.observed.and.predicted.phenotypic.values.Ad) <- c("Taxa", "Observed.Value", "Predicted.Value")
                  x.p.Ad=these.observed.and.predicted.phenotypic.values.Ad[,c(1,3)]
                  y.o.Ad=these.observed.and.predicted.phenotypic.values.Ad[,c(1,2)]
                  myCI.AD[trn,1] <- round(CI(x.p.Ad,y.o.Ad,s=s,top=T),2)
                  
                  #corr.sqr.Ad <- round(summary(the.fitted.model.AD)$adj.r.squared, 2)
                  #Ad.accuracy[trn,b] <- round((corr.sqr.Ad/herit), 2)
        
              
        
        } else{
          
          Yt.1=(pheno_trainingset.1[,1])
          Yt.2=(pheno_trainingset.2[,1])
          
          Geno_trainingset.1 <- G.1[match(rownames(pheno_trainingset.1), rownames(G.1)),]
          
          Geno_trainingset.2 <- G.2[match(rownames(pheno_trainingset.2), rownames(G.2)),]
          
          ww <- pheno_testset
          
          library(rrBLUP)
          
          rrMod.Ro.1 <-mixed.solve(Yt.1, X=NULL, Z=Geno_trainingset.1, K=NULL, SE=F, return.Hinv=F)
          mEff.Ro.1<-rrMod.Ro.1$u
          e.Ro.1= as.matrix(mEff.Ro.1)
          
          rrMod.Ro.2<-mixed.solve(Yt.2, X=NULL, Z=Geno_trainingset.2, K=NULL, SE=F, return.Hinv=F)
          mEff.Ro.2 <- rrMod.Ro.2$u
          e.Ro.2 = as.matrix(mEff.Ro.2)
            Geno_testset <- Geno_testset[match(rownames(pheno_testset), rownames(Geno_testset)),]
            
            predYv.Ro.1 = Geno_testset%*%e.Ro.1
            predYr.Ro.1 = predYv.Ro.1[,1]+ rrMod.Ro.1$beta
            
            predYv.Ro.2 = Geno_testset%*%e.Ro.2
            predYr.Ro.2 = predYv.Ro.2[,1]+ rrMod.Ro.2$beta
            
            predYr.Ro <- predYr.Ro.1+predYr.Ro.2 # Combine Msi and Msa predicted values together
            Y_valid=pheno_testset
            
            r.gy.ADE[trn,1] <-  cor(predYr.Ro,Y_valid,use="complete")
            
            the.fitted.model.ADE <- lm(Y_valid ~ predYr.Ro)
            the.coefficients.ADE.int[trn,1] <- the.fitted.model.ADE$coefficients[1]
            the.coefficients.ADE.slope[trn,1] <- the.fitted.model.ADE$coefficients[2]
            
            these.observed.and.predicted.phenotypic.values.Ade <- data.frame(rownames(ww), Y_valid, predYr.Ro)
            rownames(these.observed.and.predicted.phenotypic.values.Ade) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Ade) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,3)]
            y.o.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,2)]
            #write.table(these.observed.and.predicted.phenotypic.values.Ade, paste("these.observed.and.predicted.phenotypic.values.RRBLUP", phenames[l], trn.app[trn],"txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
            myCI.ADE[trn,1] <- round(CI(x.p.Ade,y.o.Ade,s=s,top=T),2)
            
            these.observed.and.predicted.phenotypic.values.Ade <- these.observed.and.predicted.phenotypic.values.Ade[order(these.observed.and.predicted.phenotypic.values.Ade$Predicted.Value, decreasing = T),]
           # corr.sqr.Ade <- round(summary(the.fitted.model.ADE)$adj.r.squared, 2)
            #myGSR.ADE[trn,b] <- round((corr.sqr.Ade/herit), 2)
            
            
            library(sommer)
            
            Pheno.comb.1 <- rbind(pheno_trainingset.1, pheno_testset)
            train.comb.1 <- 1:nrow(pheno_trainingset.1)
            valid.comb.1 <-setdiff(1:nrow(Pheno.comb.1),train.comb.1)
            
            y.trn.1 <- Pheno.comb.1 # for prediction accuracy
            
            Pheno.comb.2 <- rbind(pheno_trainingset.2, pheno_testset)
            train.comb.2 <- 1:nrow(pheno_trainingset.2)
            valid.comb.2 <-setdiff(1:nrow(Pheno.comb.2),train.comb.2)
            
            y.trn.2 <- Pheno.comb.2 # for prediction accuracy
            
            ww <- pheno_testset # delete data for 1/5 of the population
            
            y.trn.1[valid.comb.1,1] <- NA #valid.comb
            
            Geno.comb.1 <- rbind(Geno_trainingset.1, Geno_testset)
            Geno.comb.1 <- as.matrix(Geno.comb.1)
            A1.1 <- A.mat(Geno.comb.1,shrink=TRUE) # additive relationship matrix 
            # D.mat also take -1,0,and 1
            D1.1 <- D.mat(Geno.comb.1,shrink=TRUE)
            
            Za.1 <- diag(length(Pheno.comb.1[,1])) 
            Zd.1 <- diag(length(Pheno.comb.1[,1]))
            rownames(A1.1) <- 1:nrow(A1.1)
            rownames(D1.1) <- 1:nrow(D1.1)
            
            Y_valid=pheno_testset
            
            ETA.AD.1 <- list(add=list(Z=Za.1,K=A1.1), dom=list(Z=Zd.1,K=D1.1))
            ans.AD.1 <- MEMMA(Y=y.trn.1, ZETA=ETA.AD.1) 
            
            
            y.trn.2[valid.comb.2,1] <- NA #valid.comb
            
            Geno.comb.2 <- rbind(Geno_trainingset.2, Geno_testset)
            Geno.comb.2 <- as.matrix(Geno.comb.2)
            A1.2 <- A.mat(Geno.comb.2,shrink=TRUE) # additive relationship matrix 
            # D.mat also take -1,0,and 1
            D1.2 <- D.mat(Geno.comb.2,shrink=TRUE)
            
            Za.2 <- diag(length(Pheno.comb.2[,1])) 
            Zd.2 <- diag(length(Pheno.comb.2[,1]))
            rownames(A1.2) <- 1:nrow(A1.2)
            rownames(D1.2) <- 1:nrow(D1.2)
            
            Y_valid=pheno_testset
            
            ETA.AD.2 <- list(add=list(Z=Za.2,K=A1.2), dom=list(Z=Zd.2,K=D1.2))
            ans.AD.2 <- MEMMA(Y=y.trn.2, ZETA=ETA.AD.2) 
            
            
            predicted.1 <- ans.AD.1$fitted.y[valid.comb.1]
            predicted.2 <- ans.AD.2$fitted.y[valid.comb.2]
            predicted.total <- predicted.1+predicted.2
            r.gy.AD[trn,1] <- cor(predicted.total, ww, use="complete")
            
            the.fitted.model.AD <- lm(ww ~ predicted.total)
            the.coefficients.AD <- c(the.fitted.model.AD$coefficients[1], the.fitted.model.AD$coefficients[2])
            the.coefficients.AD.int[trn,1] <- the.coefficients.AD[1]
            the.coefficients.AD.slope[trn,1] <- the.coefficients.AD[2]
            
            these.observed.and.predicted.phenotypic.values.Ad <- data.frame(rownames(Y_valid), Y_valid, predicted.total)
            rownames(these.observed.and.predicted.phenotypic.values.Ad) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Ad) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Ad=these.observed.and.predicted.phenotypic.values.Ad[,c(1,3)]
            y.o.Ad=these.observed.and.predicted.phenotypic.values.Ad[,c(1,2)]
            myCI.AD[trn,1] <- round(CI(x.p.Ad,y.o.Ad,s=s,top=T),2)
            
            #corr.sqr.Ad <- round(summary(the.fitted.model.AD)$adj.r.squared, 2)
            #Ad.accuracy[trn,b] <- round((corr.sqr.Ad/herit), 2)
            
       
          
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
      r.gy.output.ADE$Model <- rep("RRBLUP", nrow(r.gy.output.ADE))
      r.gy.output.ADE$Parameter <- rep("P.A", nrow(r.gy.output.ADE))
      
      the.coefficients.output.ADE.int <- as.data.frame(the.coefficients.ADE.int)
      the.coefficients.output.ADE.int$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.int))
      the.coefficients.output.ADE.int$Approach <- trn.app
      the.coefficients.output.ADE.int$Model <- rep("RRBLUP", nrow(the.coefficients.output.ADE.int))
      the.coefficients.output.ADE.int$Parameter <- rep("Intercept", nrow(the.coefficients.output.ADE.int))
      
      the.coefficients.output.ADE.slope <- as.data.frame(the.coefficients.ADE.slope)
      the.coefficients.output.ADE.slope$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.slope))
      the.coefficients.output.ADE.slope$Approach <- trn.app
      the.coefficients.output.ADE.slope$Model <- rep("RRBLUP", nrow(the.coefficients.output.ADE.slope))
      the.coefficients.output.ADE.slope$Parameter <- rep("Slope", nrow(the.coefficients.output.ADE.slope))
      
      myCI.output.ADE<- as.data.frame(myCI.ADE)
      myCI.output.ADE$Trait <- rep(phenames[l], nrow(myCI.output.ADE))
      myCI.output.ADE$Approach <- trn.app
      myCI.output.ADE$Model <- rep("RRBLUP", nrow(myCI.output.ADE))
      myCI.output.ADE$Parameter <- rep("CI", nrow(myCI.output.ADE))
      

      
      #write.table(r.gy.output.ADE, paste(species,phenames[l], "r.gy.output.RRBLUP.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(the.coefficients.output.ADE.int, paste(species,phenames[l], "the.coefficients.output.RRBLUP.int.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(the.coefficients.output.ADE.slope, paste(species,phenames[l], "the.coefficients.output.RRBLUP.slope.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(myCI.output.ADE, paste(species,phenames[l], "myCI.output.RRBLUP.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      
      #r.gy.output.ADE.Traits <- rbind(r.gy.output.ADE.Traits, r.gy.output.ADE)
      #the.coefficients.output.ADE.int.Traits <- rbind(the.coefficients.output.ADE.int.Traits, the.coefficients.output.ADE.int)
      #the.coefficients.output.ADE.slope.Traits <- rbind(the.coefficients.output.ADE.slope.Traits, the.coefficients.output.ADE.slope)
      #myCI.output.ADE.Traits <- rbind(myCI.output.ADE.Traits, myCI.output.ADE)
      
      
      r.gy.output.AD <- as.data.frame(r.gy.AD)
      r.gy.output.AD$Trait <- rep(phenames[l], nrow(r.gy.output.AD))
      r.gy.output.AD$Approach <- as.vector(trn.app)
      r.gy.output.AD$Model <- rep("AD", nrow(r.gy.output.AD))
      r.gy.output.AD$Parameter <- rep("P.A", nrow(r.gy.output.AD))
      
      the.coefficients.output.AD.int <- as.data.frame(the.coefficients.AD.int)
      the.coefficients.output.AD.int$Trait <- rep(phenames[l], nrow(the.coefficients.output.AD.int))
      the.coefficients.output.AD.int$Approach <- trn.app
      the.coefficients.output.AD.int$Model <- rep("AD", nrow(the.coefficients.output.AD.int))
      the.coefficients.output.AD.int$Parameter <- rep("Intercept", nrow(the.coefficients.output.AD.int))
      
      the.coefficients.output.AD.slope <- as.data.frame(the.coefficients.AD.slope)
      the.coefficients.output.AD.slope$Trait <- rep(phenames[l], nrow(the.coefficients.output.AD.slope))
      the.coefficients.output.AD.slope$Approach <- trn.app
      the.coefficients.output.AD.slope$Model <- rep("AD", nrow(the.coefficients.output.AD.slope))
      the.coefficients.output.AD.slope$Parameter <- rep("Slope", nrow(the.coefficients.output.AD.slope))
      
      myCI.output.AD<- as.data.frame(myCI.AD)
      myCI.output.AD$Trait <- rep(phenames[l], nrow(myCI.output.AD))
      myCI.output.AD$Approach <- trn.app
      myCI.output.AD$Model <- rep("AD", nrow(myCI.output.AD))
      myCI.output.AD$Parameter <- rep("CI", nrow(myCI.output.AD))
      
      
      
      #write.table(r.gy.output.AD, paste(species,phenames[l], "r.gy.output.AD.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(the.coefficients.output.AD.int, paste(species,phenames[l], "the.coefficients.output.AD.int.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(the.coefficients.output.AD.slope, paste(species,phenames[l], "the.coefficients.output.AD.slope.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      #write.table(myCI.output.AD, paste(species,phenames[l], "myCI.output.AD.txt", sep="_"), sep="\t", quote=F, row.names = F, col.names = T)
      
      #r.gy.output.AD.Traits <- rbind(r.gy.output.AD.Traits, r.gy.output.AD)
      #the.coefficients.output.AD.int.Traits <- rbind(the.coefficients.output.AD.int.Traits, the.coefficients.output.AD.int)
      #the.coefficients.output.AD.slope.Traits <- rbind(the.coefficients.output.AD.slope.Traits, the.coefficients.output.AD.slope)
      #myCI.output.AD.Traits <- rbind(myCI.output.AD.Traits, myCI.output.AD)
    
      Compiled.Results <- rbind(r.gy.output.ADE, r.gy.output.AD, the.coefficients.output.ADE.int, the.coefficients.output.AD.int, the.coefficients.output.ADE.slope, the.coefficients.output.AD.slope, myCI.output.ADE, myCI.output.AD)
      
  } # End of trait
  
  #print("--------- Writing Final Results In Working Directory ------------")
  
  #save.image(paste(species, "GS.GS.Results.Div.5FCV.TRSvec.ADE.AND.AD.rda", sep='.'))
  
  #write.table(r.gy.output.ADE.Traits, paste(species, "r.gy.output.ADE.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(the.coefficients.output.ADE.int.Traits, paste(species, "the.coefficients.output.ADE.int.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(the.coefficients.output.ADE.slope.Traits, paste(species, "the.coefficients.output.ADE.slope.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(myCI.output.ADE.Traits, paste(species, "myCI.output.ADE.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  
  #write.table(r.gy.output.AD.Traits, paste(species, "r.gy.output.AD.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(the.coefficients.output.AD.int.Traits, paste(species, "the.coefficients.output.AD.int.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(the.coefficients.output.AD.slope.Traits, paste(species, "the.coefficients.output.AD.slope.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  #write.table(myCI.output.AD.Traits, paste(species, "myCI.output.AD.Traits", date,"txt", sep='.'), sep="\t", quote=F, row.names=F, col.names=T)
  
  return(Compiled.Results)
  
  } #
#} # End of library
# Started running 07022019 09:17 AM


