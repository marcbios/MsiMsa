
Within.Diversity.Panel.GS.Pipeline <- function(phenames=NULL, polyRAD.Geno=NULL, m_train.pheno=NULL, her.vec=NULL, dat=NULL, PCA.info=NULL, Opt.PCA=NULL,
                                               num.vec=NULL, s=NULL, Iter=NULL, burn=NULL, species=NULL, k=NULL, itrn=NULL, cycles=NULL, 
                                               trn.app=NULL, number.of.folds=NULL, geno.diplo.tetra=NULL, Groups.3.Pop3=NULL, date=NULL){
  
  r.gy.output.PCA.Traits <- NULL
  the.coefficients.output.PCA.int.Traits <- NULL
  the.coefficients.output.PCA.slope.Traits <- NULL
  myCI.output.PCA.Traits <- NULL
  myGSR.output.PCA.Traits <- NULL
  
  r.gy.output.ADE.Traits <- NULL
  the.coefficients.output.ADE.int.Traits <- NULL
  the.coefficients.output.ADE.slope.Traits <- NULL
  myCI.output.ADE.Traits <- NULL
  myGSR.output.ADE.Traits <- NULL
  
  
  CDmean.set.n.test.set.trt=NULL
  
  for(l in 1:length(phenames)){ #length(phenames)#1:length(phenames)
    
    CDmean.set.n.test.set.NULL=NULL
    r.gy.output.ADE.NULL=NULL
    the.coefficients.output.ADE.int.NULL=NULL
    the.coefficients.output.ADE.slope.NULL=NULL
    myCI.output.ADE.NULL=NULL
    myGSR.output.ADE.NULL=NULL
    
    r.gy.output.PCA.NULL=NULL
    the.coefficients.output.PCA.int.NULL=NULL
    the.coefficients.output.PCA.slope.NULL=NULL
    myCI.output.PCA.NULL=NULL
    myGSR.output.PCA.NULL=NULL
    
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], "!!!!!!!----------", sep = ""))
    
    Pheno1=as.matrix(m_train.pheno[,phenames[l]])
    rownames(Pheno1) <- m_train.pheno$Taxa
    colnames(Pheno1) <- phenames[l]
    
    G <- polyRAD.Geno
    herit <- her.vec[which(her.vec$trait==phenames[l]),1]
    
    G2 <- G
    
    r.gy.ADE <- matrix(NA, length(trn.app), k)
    myCI.ADE <- matrix(NA, length(trn.app), k)
    myGSR.ADE <- matrix(NA, length(trn.app), k)
    the.coefficients.ADE.int <- matrix(NA, length(trn.app), k)
    the.coefficients.ADE.slope <- matrix(NA, length(trn.app), k)
    
    G.PCA <- PCA.info
    G2.PCA <- G.PCA[,-1]
    rownames(G2.PCA) <- G.PCA$Taxa
    r.gy.PCA <- matrix(NA, length(trn.app), k)
    myCI.PCA <- matrix(NA, length(trn.app), k)
    myGSR.PCA <- matrix(NA, length(trn.app), k)
    the.coefficients.PCA.int <- matrix(NA, length(trn.app), k)
    the.coefficients.PCA.slope <- matrix(NA, length(trn.app), k)
    
    for(r in 1:cycles){
      
      print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles, " cycles !!!!!!!----------", sep = "")) 
      
      Pheno1.TRS1 <- as.data.frame(Pheno1)
      Pheno1.TRS1$id <- sample(1:k, nrow(Pheno1.TRS1), replace = TRUE)
      
    list <- 1:k
    
    CDmean.set.n.test.set <- matrix(NA,5,1)
    
    
    for (f in 1:k){
      print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles," --- CrossFolds ", f," !!!!!!!----------", sep = ""))
      
      # remove rows with id i from dataframe to create training set
      # select rows with id i to create test set
      #set.seed()
      pheno_trainingset <- subset(Pheno1.TRS1, id %in% list[-f])
      
      #colnames(pheno_trainingset)[1] <- phenames[l]
      Geno_trainingset <- G2[match(rownames(pheno_trainingset), rownames(G2)),]
      Geno_trainingset <- as.matrix(Geno_trainingset)
      
      
      pheno_testset <- subset(Pheno1.TRS1, id %in% list[f])
      #colnames(pheno_testset) <- phenames[l]
      Geno_testset <- G2[match(rownames(pheno_testset), rownames(G2)),]
      Geno_testset <- as.matrix(Geno_testset)
      
      # Apply CD Mean to select optimal training set from this fold
      
      Dat_train <- dat[match(rownames(Geno_trainingset), rownames(dat)),]
      Dat_train <- Dat_train-1
      
      
      A.A <- A.mat(Dat_train)#, impute.method="mean",return.imputed = TRUE
      matA1 <- A.mat(Dat_train) #, min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,shrink=FALSE,return.imputed=FALSE
      matA1=as.matrix(matA1)
      
      # Estimate variance components
      y<-pheno_trainingset[,phenames[l]]
      y[is.na(y)]<-mean(y,na.rm=T)
      model <- mixed.solve(y,K=A.A)
      
      varG = model$Vu #1.336796
      varE <- model$Ve #0.9985398
      
      print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- CD mean for CrossFolds ", f," !!!!!!!----", sep = ""))
      
      # CDmean function to return optimal 200 individuals
      sel.train.CDmean <- CDMean(matA1 = matA1, varG = varG, varE = varE, nindrep=num.vec, itrn = itrn)
      
      write.table(sel.train.CDmean, paste(species, "CDmean.Selected.Individuals.200", "Trait: ", phenames[l], "Five Folds minus fold: ", f, "txt", sep="."), sep="\t", quote=F, row.names=F, col.names = F)
      
      sel.train.CDmean.df <- data.frame(sel.train.CDmean)
      colnames(sel.train.CDmean.df) <- "Taxa"
      pheno_trainingset.2 <- pheno_trainingset
      pheno_trainingset.2$Taxa <- rownames(pheno_trainingset.2)
      sel.train.CDmean.df <- merge(sel.train.CDmean.df, pheno_trainingset.2, by="Taxa")
      sel.train.CDmean.df.2 <- sel.train.CDmean.df
      sel.train.CDmean.df.2$pops <- rep("CDmean.set", nrow(sel.train.CDmean.df.2))
      
      pheno_testset.ddf <- data.frame(pheno_testset)
      pheno_testset.ddf$Taxa <- rownames(pheno_testset.ddf)
      pheno_testset.ddf$pops <- rep("test.set", nrow(pheno_testset.ddf))
      pheno_testset.ddf[,-1]
      geno.pops <- rbind(sel.train.CDmean.df.2, pheno_testset.ddf)
      Groups.4.Pop3 <- Groups.3.Pop3[,-2]
      Groups.5.Pop3 <- merge(Groups.4.Pop3, geno.pops)
      pops <- Groups.5.Pop3$pops
      ploidy <- Groups.5.Pop3$ploidy
      
      geno.diplo.tetra2 <- geno.diplo.tetra[match(Groups.5.Pop3$Taxa, rownames(geno.diplo.tetra)),]

      res.diptetra <- pairwise.JostD.numeric(genmat=geno.diplo.tetra2, pops=pops, ploidy = ploidy)
      
      CDmean.set.n.test.set[f,1] <- mean(res.diptetra$CDmean.set_test.set, na.rm=T)
      
      # Select Random 200 individuals
      random.samp <- sample(nrow(pheno_trainingset), num.vec, replace = FALSE)
      random_train=pheno_trainingset[random.samp,]
      
      trainingset.taxa <- data.frame(rownames(pheno_trainingset))
      colnames(trainingset.taxa)[1] <- "Taxa"
      # run RRBLUP
      
      # Loop for different approaches Raw phenotypes, PCA, CDmean, DAPC
      for(trn in 1:length(trn.app)){
        
        if(trn.app[trn]=="All"){
          pheno_trainingset = pheno_trainingset; rw.nam <- rownames(pheno_trainingset);
          pheno_trainingset <- pheno_trainingset[,-2]; names(pheno_trainingset) <- rw.nam;
          rt.nam <- rownames(pheno_testset); pheno_testset <- pheno_testset[,-2]; names(pheno_testset) <- rt.nam      # 12% VAT
        } else {
          if(trn.app[trn]=="CDmean"){
            pheno_trainingset = pheno_trainingset[match(sel.train.CDmean.df$Taxa, rownames(pheno_trainingset)),]; 
            pheno_trainingset=as.matrix(pheno_trainingset); pheno_trainingset=pheno_trainingset[,-2];
            rt.nam <- rownames(pheno_testset);pheno_testset <- pheno_testset[,-2]; names(pheno_testset) <- rt.nam
          } else { 
            if(trn.app[trn]=="Random"){pheno_trainingset = pheno_trainingset[match(rownames(random_train), rownames(pheno_trainingset)),]; 
            pheno_trainingset=as.matrix(pheno_trainingset); pheno_trainingset=pheno_trainingset[,-2];
            rt.nam <- rownames(pheno_testset);pheno_testset <- pheno_testset[,-2]; names(pheno_testset) <- rt.nam
            } else{
              pheno_trainingset = pheno_trainingset; rw.nam <- rownames(pheno_trainingset);
              pheno_trainingset <- pheno_trainingset[,-2]; names(pheno_trainingset) <- rw.nam;
              rt.nam <- rownames(pheno_testset); pheno_testset <- pheno_testset[,-2]; names(pheno_testset) <- rt.nam}   # 0% VAT 
            }
          }
        
        #}  
          
        
        head(pheno_trainingset)
        head(pheno_testset)
      
        pheno_trainingset <- as.matrix(pheno_trainingset)
        colnames(pheno_trainingset) <- phenames[l]
        
        pheno_testset <- as.matrix(pheno_testset)
        colnames(pheno_testset) <- phenames[l]
        
        head(pheno_trainingset)
        head(pheno_testset)
        
        Yt=(pheno_trainingset[,1])
        
         
        #rw.test.nam <- rownames(pheno_testset)
        #pheno_testset <- as.matrix(pheno_testset)
        #pheno_testset <- pheno_testset[,-2] #remove fold id column
        #pheno_testset <- as.matrix(pheno_testset)
        #colnames(pheno_testset) <- phenames[l]
        
        ww <- pheno_testset # delete data for 1/5 of the population
        
        Geno_trainingset <- Geno_trainingset[match(rownames(pheno_trainingset), rownames(Geno_trainingset)),]
        
        Geno_trainingset.PCA <- G2.PCA[match(rownames(pheno_trainingset), rownames(G2.PCA)),c(1:Opt.PCA)]
        Geno_trainingset.PCA <- as.matrix(Geno_trainingset.PCA)
        
        Geno_testset.PCA <- G2.PCA[match(rownames(pheno_testset), rownames(G2.PCA)),c(1:Opt.PCA)]
        Geno_testset.PCA <- as.matrix(Geno_testset.PCA)
        
        print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles," --- CrossFolds ", f," RRBLUP GS ", trn.app[trn], sep = ""))
        rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=Geno_trainingset, K=NULL, SE=F, return.Hinv=F)
        mEff.Ro<-rrMod.Ro$u
        e.Ro= as.matrix(mEff.Ro)
        
        predYv.Ro = Geno_testset%*%e.Ro
        
        predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
        
        Y_valid=pheno_testset

        r.gy.ADE[trn,f] <-  cor(predYr.Ro,Y_valid,use="complete")
        
        the.fitted.model.ADE <- lm(Y_valid ~ predYr.Ro)
        the.coefficients.ADE.int[trn,f] <- the.fitted.model.ADE$coefficients[1]
        the.coefficients.ADE.slope[trn,f] <- the.fitted.model.ADE$coefficients[2]
        
        these.observed.and.predicted.phenotypic.values.Ade <- data.frame(rownames(ww), Y_valid, predYr.Ro)
        rownames(these.observed.and.predicted.phenotypic.values.Ade) <- NULL
        colnames(these.observed.and.predicted.phenotypic.values.Ade) <- c("Taxa", "Observed.Value", "Predicted.Value")
        x.p.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,3)]
        y.o.Ade=these.observed.and.predicted.phenotypic.values.Ade[,c(1,2)]
        write.table(these.observed.and.predicted.phenotypic.values.Ade, paste("these.observed.and.predicted.phenotypic.values.Ade", phenames[l], trn.app[trn],"txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
        myCI.ADE[trn,f] <- round(CI(x.p.Ade,y.o.Ade,s=s,top=T),2)
        
        these.observed.and.predicted.phenotypic.values.Ade <- these.observed.and.predicted.phenotypic.values.Ade[order(these.observed.and.predicted.phenotypic.values.Ade$Predicted.Value, decreasing = T),]
        #sel.taxa.Ade <- these.observed.and.predicted.phenotypic.values.Ade$Taxa[c(1:(round(length(these.observed.and.predicted.phenotypic.values.Ade$Taxa)*0.1,0)))]
        #sel.taxa.Ade <- as.character(sel.taxa.Ade)
        #sel.taxa.Ade.numvec <- c(sel.taxa.Ade.numvec, sel.taxa.Ade)
        corr.sqr.Ade <- round(summary(the.fitted.model.ADE)$adj.r.squared, 2)
        myGSR.ADE[trn,f] <- round((corr.sqr.Ade/herit), 2)
        
        print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles," --- CrossFolds ", f," PCA Prediction ", trn.app[trn], sep = ""))
        
        Geno_trainingset.df.PCA <- data.frame(Geno_trainingset.PCA)
        Geno_trainingset.df.PCA$Taxa <- rownames(Geno_trainingset.PCA)
        
        Geno_testset.df.PCA <- data.frame(Geno_testset.PCA)
        Geno_testset.df.PCA$Taxa <- rownames(Geno_testset.PCA)
        
        
        pheno_trainingset.df.PCA <- data.frame(pheno_trainingset)
        pheno_trainingset.df.PCA$Taxa <- rownames(pheno_trainingset)
        pheno_trainingset.df.PCA <- pheno_trainingset.df.PCA[,c(2,1)]
        
        pheno_testset.df.PCA <- data.frame(pheno_testset)
        pheno_testset.df.PCA$Taxa <- rownames(pheno_testset)
        pheno_testset.df.PCA <- pheno_testset.df.PCA[,c(2,1)]
        #X <- Geno.comb
        
        data.for.lm.train <- merge(pheno_trainingset.df.PCA, Geno_trainingset.df.PCA)
        data.for.lm.pred <- Geno_testset.df.PCA
        data.for.lm.pred <- data.for.lm.pred[,c(ncol(data.for.lm.pred),1:(ncol(data.for.lm.pred)-1))]
        
        equation.for.lm <- paste(colnames(data.for.lm.train)[2], "~", colnames(data.for.lm.train)[3],sep = "")
        
        for(index in 4:ncol(data.for.lm.train)) equation.for.lm <- paste(equation.for.lm,colnames(data.for.lm.train)[index],sep = "+")
        
        lm.model.fitted.in.train <- lm(equation.for.lm, data = data.for.lm.train)
        the.predicted.MAS.values <- predict(lm.model.fitted.in.train, newdata = data.for.lm.pred)
        
        these.observed.and.predicted.phenotypic.values <- data.frame(pheno_testset.df.PCA, the.predicted.MAS.values)
        colnames(these.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
        
        #Measure correclation between OBS and Pred in validation set (V.S.)
        r.gy.PCA[trn,f] <- cor(these.observed.and.predicted.phenotypic.values[,3], these.observed.and.predicted.phenotypic.values[,2],use="complete")
        x.p.MAS=these.observed.and.predicted.phenotypic.values[,c(1,3)]
        y.o.MAS=these.observed.and.predicted.phenotypic.values[,c(1,2)]
        myCI.PCA[trn,f] <- round(CI(x.p.MAS,y.o.MAS,s=s,top=T),2)
        the.fitted.regression.model <- lm(these.observed.and.predicted.phenotypic.values[,2] ~ these.observed.and.predicted.phenotypic.values[,3])
        the.coefficients.PCA.int[trn,f] <- the.fitted.regression.model$coefficients[1]
        the.coefficients.PCA.slope[trn,f] <- the.fitted.regression.model$coefficients[2]
        
        
      }# Raw phenotype, end of CDmean, end of PCA, end of DAPC
      
      pheno_trainingset=NULL
      pheno_testset=NULL
      Geno_trainingset=NULL
      Geno_testset=NULL
      

      pheno_trainingset.PCA=NULL
      pheno_testset.PCA=NULL
      Geno_trainingset.PCA=NULL
      Geno_testset.PCA=NULL
      
    } # End of i in k cross valid
    
    CDmean.set.n.test.set <- data.frame(CDmean.set.n.test.set)
    CDmean.set.n.test.set.NULL <- rbind(CDmean.set.n.test.set.NULL, CDmean.set.n.test.set)
    
    r.gy.output.ADE <- as.data.frame(r.gy.ADE)
    r.gy.output.ADE$Approach <- as.vector(trn.app)
    r.gy.output.ADE$Cycle <- rep(r, nrow(r.gy.output.ADE))
    r.gy.output.ADE.NULL <- rbind(r.gy.output.ADE.NULL, r.gy.output.ADE)
    
    the.coefficients.output.ADE.int <- as.data.frame(the.coefficients.ADE.int)
    the.coefficients.output.ADE.int$Approach <- trn.app
    the.coefficients.output.ADE.int$Cycle <- rep(r, nrow(the.coefficients.output.ADE.int))
    the.coefficients.output.ADE.int.NULL <- rbind(the.coefficients.output.ADE.int.NULL, the.coefficients.output.ADE.int)

    the.coefficients.output.ADE.slope <- as.data.frame(the.coefficients.ADE.slope)
    the.coefficients.output.ADE.slope$Approach <- trn.app
    the.coefficients.output.ADE.slope$Cycle <- rep(r, nrow(the.coefficients.output.ADE.slope))
    the.coefficients.output.ADE.slope.NULL <- rbind(the.coefficients.output.ADE.slope.NULL, the.coefficients.output.ADE.slope)
    
    myCI.output.ADE<- as.data.frame(myCI.ADE)
    myCI.output.ADE$Approach <- trn.app
    myCI.output.ADE$Cycle <- rep(r, nrow(myCI.output.ADE))
    myCI.output.ADE.NULL <- rbind(myCI.output.ADE.NULL, myCI.output.ADE)
    
    myGSR.output.ADE <- as.data.frame(myGSR.ADE)
    myGSR.output.ADE$Approach <- trn.app
    myGSR.output.ADE$Cycle <- rep(r, nrow(myGSR.output.ADE))
    myGSR.output.ADE.NULL <- rbind(myGSR.output.ADE.NULL, myGSR.output.ADE)
    
    
    r.gy.output.PCA <- as.data.frame(r.gy.PCA)
    r.gy.output.PCA$Approach <- as.vector(trn.app)
    r.gy.output.PCA$Cycle <- rep(r, nrow(r.gy.output.PCA))
    r.gy.output.PCA.NULL <- rbind(r.gy.output.PCA.NULL, r.gy.output.PCA)
    
    the.coefficients.output.PCA.int <- as.data.frame(the.coefficients.PCA.int)
    the.coefficients.output.PCA.int$Approach <- trn.app
    the.coefficients.output.PCA.int$Cycle <- rep(r, nrow(the.coefficients.output.PCA.int))
    the.coefficients.output.PCA.int.NULL <- rbind(the.coefficients.output.PCA.int.NULL, the.coefficients.output.PCA.int)
    
    the.coefficients.output.PCA.slope <- as.data.frame(the.coefficients.PCA.slope)
    the.coefficients.output.PCA.slope$Approach <- trn.app
    the.coefficients.output.PCA.slope$Cycle <- rep(r, nrow(the.coefficients.output.PCA.slope))
    the.coefficients.output.PCA.slope.NULL <- rbind(the.coefficients.output.PCA.slope.NULL, the.coefficients.output.PCA.slope)
    
    myCI.output.PCA<- as.data.frame(myCI.PCA)
    myCI.output.PCA$Approach <- trn.app
    myCI.output.PCA$Cycle <- rep(r, nrow(myCI.output.PCA))
    myCI.output.PCA.NULL <- rbind(myCI.output.PCA.NULL, myCI.output.PCA)
    
    myGSR.output.PCA <- as.data.frame(myGSR.PCA)
    myGSR.output.PCA$Approach <- trn.app
    myGSR.output.PCA$Cycle <- rep(r, nrow(myGSR.output.PCA))
    myGSR.output.PCA.NULL <- rbind(myGSR.output.PCA.NULL, myGSR.output.PCA)


    
    } # End of cycle
    
    CDmean.set.n.test.set.NULL$trait <- rep(phenames[l], nrow(CDmean.set.n.test.set.NULL))
    CDmean.set.n.test.set.trt <- rbind(CDmean.set.n.test.set.trt, CDmean.set.n.test.set.NULL)
    CDmean.set.n.test.set=NULL
    
    r.gy.output.ADE.NULL$Trait <- rep(phenames[l], nrow(r.gy.output.ADE.NULL))
    the.coefficients.output.ADE.int.NULL$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.int.NULL))
    the.coefficients.output.ADE.slope.NULL$Trait <- rep(phenames[l], nrow(the.coefficients.output.ADE.slope.NULL))
    myCI.output.ADE.NULL$Trait <- rep(phenames[l], nrow(myCI.output.ADE.NULL))
    myGSR.output.ADE.NULL$Trait <- rep(phenames[l], nrow(myGSR.output.ADE.NULL))
    
    
    write.table(r.gy.output.ADE.NULL, paste(species,phenames[l], "r.gy.output.RRBLUP.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(the.coefficients.output.ADE.int.NULL, paste(species,phenames[l], "the.coefficients.output.RRBLUP.int.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(the.coefficients.output.ADE.slope.NULL, paste(species,phenames[l], "the.coefficients.output.RRBLUP.slope.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(myCI.output.ADE.NULL, paste(species,phenames[l], "myCI.output.RRBLUP.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(myGSR.output.ADE.NULL, paste(species,phenames[l], "myGSR.output.RRBLUP.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    
    r.gy.output.ADE.Traits <- rbind(r.gy.output.ADE.Traits, r.gy.output.ADE.NULL)
    the.coefficients.output.ADE.int.Traits <- rbind(the.coefficients.output.ADE.int.Traits, the.coefficients.output.ADE.int.NULL)
    the.coefficients.output.ADE.slope.Traits <- rbind(the.coefficients.output.ADE.slope.Traits, the.coefficients.output.ADE.slope.NULL)
    myCI.output.ADE.Traits <- rbind(myCI.output.ADE.Traits, myCI.output.ADE.NULL)
    myGSR.output.ADE.Traits <- rbind(myGSR.output.ADE.Traits, myGSR.output.ADE.NULL)
    
    
    r.gy.output.PCA.NULL$Trait <- rep(phenames[l], nrow(r.gy.output.PCA.NULL))
    the.coefficients.output.PCA.int.NULL$Trait <- rep(phenames[l], nrow(the.coefficients.output.PCA.int.NULL))
    the.coefficients.output.PCA.slope.NULL$Trait <- rep(phenames[l], nrow(the.coefficients.output.PCA.slope.NULL))
    myCI.output.PCA.NULL$Trait <- rep(phenames[l], nrow(myCI.output.PCA.NULL))
    myGSR.output.PCA.NULL$Trait <- rep(phenames[l], nrow(myGSR.output.PCA.NULL))
    
    
    write.table(r.gy.output.PCA.NULL, paste(species,phenames[l], "r.gy.output.PCA.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(the.coefficients.output.PCA.int.NULL, paste(species,phenames[l], "the.coefficients.output.PCA.int.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(the.coefficients.output.PCA.slope.NULL, paste(species,phenames[l], "the.coefficients.output.PCA.slope.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(myCI.output.PCA.NULL, paste(species,phenames[l], "myCI.output.PCA.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    write.table(myGSR.output.PCA.NULL, paste(species,phenames[l], "myGSR.output.PCA.csv", sep="_"), sep=",", quote=F, row.names = F, col.names = T)
    
    r.gy.output.PCA.Traits <- rbind(r.gy.output.PCA.Traits, r.gy.output.PCA.NULL)
    the.coefficients.output.PCA.int.Traits <- rbind(the.coefficients.output.PCA.int.Traits, the.coefficients.output.PCA.int.NULL)
    the.coefficients.output.PCA.slope.Traits <- rbind(the.coefficients.output.PCA.slope.Traits, the.coefficients.output.PCA.slope.NULL)
    myCI.output.PCA.Traits <- rbind(myCI.output.PCA.Traits, myCI.output.PCA.NULL)
    myGSR.output.PCA.Traits <- rbind(myGSR.output.PCA.Traits, myGSR.output.PCA.NULL)
  } # End of trait
  
  print("--------- Writing Final REsults In Working Directory ------------")
  
  save.image(paste(species, "GS.GS.Results.Div.5FCV.TRSvec.rda", sep='.'))
  
  write.table(CDmean.set.n.test.set.trt, paste(species, "CDmean.set.n.test.set.trt", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  
  write.table(r.gy.output.ADE.Traits, paste(species, "r.gy.output.RRBLUP.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.int.Traits, paste(species, "the.coefficients.output.RRBLUP.int.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.slope.Traits, paste(species, "the.coefficients.output.RRBLUP.slope.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(myCI.output.ADE.Traits, paste(species, "myCI.output.RRBLUP.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(myGSR.output.ADE.Traits, paste(species, "myGSR.output.RRBLUP.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  

  write.table(r.gy.output.ADE.Traits, paste(species, "r.gy.output.PCA.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.int.Traits, paste(species, "the.coefficients.output.PCA.int.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(the.coefficients.output.ADE.slope.Traits, paste(species, "the.coefficients.output.PCA.slope.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(myCI.output.ADE.Traits, paste(species, "myCI.output.PCA.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(myGSR.output.ADE.Traits, paste(species, "myGSR.output.PCA.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  
  
} # End of library
  # Started running 05302019 09:17 AM


