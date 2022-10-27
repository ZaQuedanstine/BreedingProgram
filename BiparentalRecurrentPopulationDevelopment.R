
########Biparental/Pedigree Recurrent Population Breeding#########
#Simulation with ALPHA sim
library(AlphaSimR)
library(parallel)
library(ggplot2)
library("adegenet")
library("pegas")
library("hierfstat")
library(SNPRelate)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")
library(statgenGWAS)
inslibrary(sommer)



sim.dir="~/Google Drive/My Drive/Brian lab folder/HPB/Simulation"
setwd(sim.dir)
source(paste(sim.dir,"/SimulationFunctions.R",sep=""))
####simulation parameters ####
narrowhert=0.5
Ne=24 #starting effect population size
n.sites=1000
n.QTL=100 #per chr
n.QTL.trait2=1
n.gen=20 #number of generations/cycles of program after initial population development 
nInd=24 #number of founder individuals # we will want to vary this and compare detection rates
intial.crosses=500 #how many crosses to make initial population
m.SI=0.1 #male selection intensity for random crossing
f.SI=0.05 #female selection intensity for random crossing
Pairwise.SI=0.05 #selection intensity for pairwise crossing scheme
random.size=1000 #sample size for GWAS
gwas.threshold=0.05 #threshold for significance for calling something a hit
distance=100#distance on same chr for what we call a hit
runs=100#total iteration
enviroment.shift.generation=10 #what generation does trait 2 come under selection
run.gwas=T #weather to run GWAS or not
gwas.gen = 10 # generation to run GWAS at
inbreeding="DH" #for biparental selection do you want to self (self) or make doubled haploids (DH) #can be none if no inbreeding desired
AF=T#do you want to get MAF of each marker for each generation?
PopStr=T #get population structure descriptors each generation
maf.metric="Summary"

##run##
for(r in 1:runs){
  print(paste("Starting Run",r,sep=" "))
  start.time=Sys.time()
  #founder/intial population development
  founderPop = runMacs(nInd=nInd,nChr=10,segSites=n.sites,
                       inbred=TRUE,manualCommand = paste(
                         "1000000000 -t", #Physical length 1e8 base pairs
                         2.5/1E8*(4*Ne), #Mutation rate adjusted for Ne
                         "-r",1/1E8*(4*Ne), #Recombination rate adjusted for Ne
                         "-eN",10/(4*Ne),100/Ne), #Modeling Ne=100 at 10 generations ago
                       manualGenLen=1) #Genetic length 1 Morgan
  SP = SimParam$new(founderPop)
  #sexes
  #SP$setSexes("yes_rand")
  SP$addTraitA(nQtlPerChr=n.QTL)
  SP$addTraitA(nQtlPerChr=n.QTL.trait2)
  SP$addSnpChip(0.1*n.sites,minSnpFreq = 0.05) #add SNPs for 10% of seg sites #also want to vary this
  SP$setVarE(h2=c(narrowhert,0.8))
  SP$nThreads=detectCores()
  
  #create all pairwise crossing scheme 
  for(w in 1:25){
    if(w==1){
      crossPlan = matrix(c(rep(1,24),2:25), nrow=24, ncol=2)
    }else{
      plan.temp=matrix(c(rep(w,24),c(1:25)[-w]), nrow=24, ncol=2)
      crossPlan=rbind(crossPlan,plan.temp)  
    }
  }
  
  #loop through generations
  print(paste("Running",n.gen,"Generations"))
  for(g in 0:n.gen){
    if(g==0){
      #print(paste("Initial Breeding Population Development Generation"))
      pop = newPop(founderPop)
      #initial population
      #random crosses from founders
      pop=randCross(pop=pop, nCrosses=intial.crosses,nProgeny = 10,simParam = SP)
      
      genMean_t1=c(meanG(pop)[1])
      phenMean_t1 = c(meanP(pop)[1])
      genMean_t2=c(meanG(pop)[2])
      phenMean_t2 = c(meanP(pop)[2])
      #genetic variance/diversity
      Gvar_t1=c(varA(pop,simParam = SP)[1,1])
      Pvar_t1=c(varP(pop)[1,1])
      Gvar_t2=c(varA(pop,simParam = SP)[2,2])
      Pvar_t2=c(varP(pop)[2,2])
      
      #get initial pop str descriptions
      if(PopStr){
        SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
        pheno=pheno(pop)
        pheno=data.frame(genotype=row.names(SNP),pheno)
        SNP.Map=getSnpMap(snpChip=1, simParam=SP)
        SNP.Map=SNP.Map[,-4]
        colnames(SNP.Map)=c("Locus", "LG" ,"Position")
        
        BG.LD=data.frame(gen=g,backgroundLD(map=SNP.Map,geno=SNP,rep=10,n.SNP=10))
        rm(SNP,pheno,SNP.Map)
      }
      ####initial MAF
      #get minor allele freq 
      if(AF){
        SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
        maf=data.frame(g=g,MAF(SNP,return=maf.metric))
        rm(SNP)
      }
    
      pop=selectInd(pop,nInd=pop@nInd*Pairwise.SI,trait=1, use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
      
      #inbreeding/selfing
      if(inbreeding=="self"){
        #print(paste("Selfing"))
        for(i in 1:5){
          pop=self(pop, nProgeny = 1, parents = NULL, keepParents = FALSE, simParam = SP)
        }
      } else if(inbreeding=="DH"){
        #print(paste("Making DH"))
        pop=makeDH(pop, nDH = 1, useFemale = TRUE, keepParents = FALSE, simParam = SP)
      } ######what does nDH=1 mean?
      
      #selection 25 top DHs
      pop=selectInd(pop,nInd=25,trait=1, use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
      gc()
    } #end of if g==0 
    
    pop = makeCross(pop,crossPlan=crossPlan, nProgeny=10,simParam=SP)
    
    #record population pheno and genetic means
    #phenotypic/Genotypic mean
    genMean_t1=c(genMean_t1,meanG(pop)[1])
    phenMean_t1 = c(phenMean_t1,meanP(pop)[1])
    genMean_t2=c(genMean_t2,meanG(pop)[2])
    phenMean_t2 = c(phenMean_t2,meanP(pop)[2])
    #genetic variance/diversity
    Gvar_t1=c(Gvar_t1,varA(pop,simParam = SP)[1,1])
    Pvar_t1=c(Pvar_t1,varP(pop)[1,1])
    Gvar_t2=c(Gvar_t2,varA(pop,simParam = SP)[2,2])
    Pvar_t2=c(Pvar_t2,varP(pop)[2,2])
    #get pop str descriptions
    if(PopStr){
      SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
      pheno=pheno(pop)
      pheno=data.frame(genotype=row.names(SNP),pheno)
      SNP.Map=getSnpMap(snpChip=1, simParam=SP)
      SNP.Map=SNP.Map[,-4]
      colnames(SNP.Map)=c("Locus", "LG" ,"Position")
      
      BG.LD=rbind(BG.LD,data.frame(gen=g+1,backgroundLD(map=SNP.Map,geno=SNP,rep=10,n.SNP=10)))
      rm(SNP,pheno,SNP.Map)
    }
    #get minor allele freq 
    if(AF){
      SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
      maf=rbind(maf,data.frame(g=g+1,MAF(SNP,return=maf.metric)))
      rm(SNP)
    }
    #remove geno for memory 
    ##Run GWAS
    if(run.gwas && g==gwas.gen){
      QTLmap=getQtlMap(trait=2,simParam = SP)
      #pull QTL states
      QTLgeno=pullQtlGeno(pop, trait = 2, chr = NULL, asRaw = FALSE, simParam =SP)
      #get MAF
      ns = nrow(QTLgeno)
      ss = apply(QTLgeno, 2, sum)
      maf.matrix = rbind((.5*ss/ns), (1-(0.5*ss/ns)))
      maf.temp = apply(maf.matrix, 2, min)
      length(which(maf.temp==0))/c(n.QTL.trait2*10)
      if(r==1){
        QTL.Loss.rate=length(which(maf.temp==0))/c(n.QTL.trait2*10)
      } else{
        QTL.Loss.rate=c(QTL.Loss.rate,length(which(maf.temp==0))/c(n.QTL.trait2*10))
      }
      
      SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
      pheno=pheno(pop)
      pheno=data.frame(genotype=row.names(SNP),pheno)
      n.ind=nrow(pheno)
      SNP.Map=getSnpMap(snpChip=1, simParam=SP)
      SNP.Map=SNP.Map[,-4]
      colnames(SNP.Map)=c("SNP.names" ,"chr" ,"pos")
      SNP.Map$allele1="A"
      SNP.Map$allele2="T"
      
      SNP.Map$chr=as.integer(SNP.Map$chr)
      SNP.Map$SNP.names=paste("X",SNP.Map$SNP.names,sep="")
      
      print(paste("Running Gwas at generation ",g,sep=""))
      gwas=sim.gwas(SNP=SNP,n.ind=pop@nInd,random.size=1000,p=gwas.threshold,distance)
      
      if(r==1){
        hit.rate=c(gwas/c(n.QTL.trait2*10))
      } else{
        hit.rate=c(hit.rate,gwas/c(n.QTL.trait2*10))
      }
    }
    
    #select top individuals
    if(g<enviroment.shift.generation){
      pop=selectInd(pop,nInd=pop@nInd*Pairwise.SI,trait=1,use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
    } else { 
      if(g==enviroment.shift.generation){
        print("Environment Shift to trait 2")
      }
      #selection on index (Trait 1, Trait 2)*b; b=(0.5,0.9)
      pop=selectInd(pop,nInd=pop@nInd*Pairwise.SI,trait=selIndex(Y =matrix(c(1,2),nrow=1), b=c(0.5,0.9), scale = FALSE),use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
    }
    #inbreeding/Selfing/DHs
    if(inbreeding=="self"){
      for(i in 1:5){
        #print(paste("Self",i))
        pop=self(pop, nProgeny = 1, parents = NULL, keepParents = FALSE, simParam = SP)
      }
    } else if(inbreeding=="DH"){
      #print(paste("Making DH"))
      pop=makeDH(pop, nDH = 1, useFemale = TRUE, keepParents = FALSE, simParam = SP)
    }
    #selection 25 top DHs
    pop=selectInd(pop,nInd=25,trait=1, use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
    
    gc()
  }
  
  #save information as recurrent population items (rec)
  if(r==1){
    biph2_t1=data.frame(run=r,gen=0:c(n.gen+1),h2_t1=Gvar_t1/(Pvar_t1+Gvar_t1))
    bipGvar_t1=data.frame(run=r,gen=0:c(n.gen+1),Gvar_t1)
    bipPvar_t1=data.frame(run=r,gen=0:c(n.gen+1),Pvar_t1)
    bipgenMean_t1=data.frame(run=r,gen=0:c(n.gen+1),genMean_t1)
    bipphenMean_t1=data.frame(run=r,gen=0:c(n.gen+1),phenMean_t1)
    biph2_t2=data.frame(h2_t2=Gvar_t2/(Pvar_t2+Gvar_t2))
    bipGvar_t2=data.frame(Gvar_t2)
    bipPvar_t2=data.frame(Pvar_t2)
    bipgenMean_t2=data.frame(genMean_t2)
    bipphenMean_t2=data.frame(phenMean_t2)
    if(AF){
      bipMAF=data.frame(r=r,maf)
    }
    if(PopStr){
      biPopStr=data.frame(r=r,BG.LD)
    }
  } else {
    biph2_t1=rbind(biph2_t1,data.frame(run=r,gen=0:c(n.gen+1),h2_t1=Gvar_t1/(Pvar_t1+Gvar_t1)))
    bipGvar_t1=rbind(bipGvar_t1,data.frame(run=r,gen=0:c(n.gen+1),Gvar_t1))
    bipPvar_t1=rbind(bipPvar_t1,data.frame(run=r,gen=0:c(n.gen+1),Pvar_t1))
    bipgenMean_t1=rbind(bipgenMean_t1,data.frame(run=r,gen=0:c(n.gen+1),genMean_t1))
    bipphenMean_t1=rbind(bipphenMean_t1,data.frame(run=r,gen=0:c(n.gen+1),phenMean_t1))
    biph2_t2=rbind(biph2_t2,data.frame(h2_t2=Gvar_t2/(Pvar_t2+Gvar_t2)))
    bipGvar_t2=rbind(bipGvar_t2,data.frame(Gvar_t2))
    bipPvar_t2=rbind(bipPvar_t2,data.frame(Pvar_t2))
    bipgenMean_t2=rbind(bipgenMean_t2,data.frame(genMean_t2))
    bipphenMean_t2=rbind(bipphenMean_t2,data.frame(phenMean_t2))
    if(AF){
      bipMAF=rbind(bipMAF,data.frame(r=r,maf))
    }
    if(PopStr){
      biPopStr=rbind(biPopStr,data.frame(r=r,BG.LD))
    }
  }
  
  end.time=Sys.time()
  print(end.time-start.time)
  
  if(r==runs){
    breeding_results=list(
      h2=cbind(biph2_t1,biph2_t2),
      Gvar=cbind(bipGvar_t1,bipGvar_t2),
      Pvar=cbind(bipPvar_t1,bipPvar_t2),
      genMean=cbind(bipgenMean_t1,bipgenMean_t2),
      phenMean=cbind(bipphenMean_t1,bipphenMean_t2)
    )
    if(run.gwas){
      b.hit.rate=hit.rate
    }
    b.QTL.Loss.rate=QTL.Loss.rate
    if(run.gwas){
      save(breeding_results,b.QTL.Loss.rate,b.hit.rate,file="BiparentalSelectionProgram.Rdata")
    } else {
      save(breeding_results,file="BiparentalSelectionProgram.Rdata")
    } 
    if(AF){
      write.table(bipMAF,file="BiparentalSelectionProgram.MAF.Results.txt",row.names = F,col.names = T,quote=F)
    }
    if(PopStr){
      write.table(biPopStr,file="BiparentalSelectionProgram.PopStr.Results.txt",row.names = F,col.names = T,quote=F)
    }
    rm(QTLgeno,QTLmap,plan.temp,crossPlan,founderPop,pop,b.hit.rate,hit.rate,bipMAF,biPopStr,BG.LD,
       biph2_t1,biph2_t2,bipGvar_t1,bipGvar_t2,bipPvar_t1,bipPvar_t2,bipgenMean_t1,bipgenMean_t2,bipphenMean_t1,bipphenMean_t2)
  }
} ######## for(r in 1:runs){ ######
