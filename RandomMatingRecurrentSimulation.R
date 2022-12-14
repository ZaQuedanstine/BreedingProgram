#Simulation with ALPHA sim
library(AlphaSimR)
library(parallel)
library(ggplot2)
library("adegenet")
library("pegas")
library("hierfstat")
library(SNPRelate)
library(statgenGWAS)
library(sommer)
library(foreach)
library(doMC)
c1 = detectCores()
registerDoMC(c1)
#  library(doParallel)
# cl = makeCluster(detectCores())
# registerDoParallel(cl)



# compiles SimulationFunctions, found in current working directory
source(paste(getwd(),"/SimulationFunctions.R",sep=""))
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
runs=5#total iteration
enviroment.shift.generation=10 #what generation does trait 2 come under selection
run.gwas=F #whether to run GWAS or not
gwas.gen = 10 # generation to run GWAS at
AF=T#do you want to get MAF of each marker for each generation?
PopStr=T #get population structure descriptors each generation
maf.metric="Summary"

start.time=Sys.time()



data <- foreach(r = 1:runs, .combine="rbind") %do%{
  print(paste("Starting Run",r,sep=" "))
#founder population
  founderPop = runMacs(nInd=nInd,nChr=10,segSites=n.sites,
                     inbred=TRUE,manualCommand = paste(
                       "1000000000 -t", #Physical length 1e8 base pairs
                       2.5/1E8*(4*Ne), #Mutation rate adjusted for Ne
                       "-r",1/1E8*(4*Ne), #Recombination rate adjusted for Ne
                       "-eN",10/(4*Ne),100/Ne), #Modeling Ne=100 at 10 generations ago
                     manualGenLen=1) #Genetic length 1 Morgan  
  SP <<- SimParam$new(founderPop)
  #sexes
  SP$setSexes("yes_rand")
  SP$addTraitA(nQtlPerChr=n.QTL)
  SP$addTraitA(nQtlPerChr=n.QTL.trait2)
  SP$addSnpChip(0.1*n.sites,minSnpFreq = 0.05) #add SNPs for 10% of seg sites #also want to vary this
  SP$setVarE(h2=c(narrowhert,0.8))
  SP$nThreads=detectCores()
### reccurent selection breeding program ###
#loop through generations
print(paste("Running",n.gen,"Generations"))
hit.rate_temp = 0
QTL.Loss.rate_temp = 0
for(g in 0:n.gen){
  #print(paste("Generation",g))
  if(g==0){
    pop = newPop(founderPop)
    pop = randCross(pop=pop, nCrosses=intial.crosses,nProgeny = 10,simParam = SP)
    #phenotypic/Genotypic mean
    genMean_t1=c(meanG(pop)[1])
    phenMean_t1 = c(meanP(pop)[1])
    genMean_t2=c(meanG(pop)[2])
    phenMean_t2 = c(meanP(pop)[2])
    #genetic variance/diversity
    Gvar_t1=c(varA(pop,simParam = SP)[1,1])
    Pvar_t1=c(varP(pop)[1,1])
    Gvar_t2=c(varA(pop,simParam = SP)[2,2])
    Pvar_t2=c(varP(pop)[2,2])
    #get intial pop str descriptions
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
    #remove geno for memory 
    gc()
  } #end of  if(g==0)
  
#phenotypically select 5% of females
female=pop[isFemale(pop)]
#print(paste("females ",female@nInd))
male=pop[isMale(pop)]
#print(paste("Males ",male@nInd))
if(g<enviroment.shift.generation){
female=selectInd(female,nInd=2000*f.SI,trait=1,use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
male=selectInd(male,nInd=2000*m.SI,trait=1,use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
} else {
  if(g==enviroment.shift.generation){
  print("Environment Shift to trait 2")
  }
  female=selectInd(female,nInd=2000*f.SI,trait=selIndex(Y =matrix(c(1,2),nrow=1), b=c(0.5,0.9), scale = FALSE),use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
  male=selectInd(male,nInd=2000*m.SI,trait=selIndex(Y =matrix(c(1,2),nrow=1), b=c(0.5,0.9), scale = FALSE),use = "pheno",simParam=SP, selectTop = TRUE, returnPop = TRUE)
}
#merge populations
pop=mergePops(list(female,male))

rm(female,male)
#random cross 5 males to each female to form initial population
m=which(isMale(pop)==T)
f=which(isFemale(pop)==T)
#create crossing plan
for(w in 1:length(f)){
  if(w==1){
    crossPlan = matrix(c(rep(1,5),sample(m,5,replace=F)), nrow=5, ncol=2)
  }else{
    plan.temp=matrix(c(rep(w,5),sample(m,5,replace=F)), nrow=5, ncol=2)
    crossPlan=rbind(crossPlan,plan.temp)  
  }
}
#make crosses and select 10 seed from each cross
pop = makeCross(pop,crossPlan=crossPlan, nProgeny=10,simParam=SP)
#report:
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
    ####Obtain the mafs of all SNPs###
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
 
  #population structure
    #LDdecay
  ##Run GWAS
  if(run.gwas && g==gwas.gen){
    QTLmap=getQtlMap(trait=2,simParam = SP)
    #pull QTL states and get proportion lost at time of environmental stress
    QTLgeno=pullQtlGeno(pop, trait = 2, chr = NULL, asRaw = FALSE, simParam =SP)
    #get MAF
    ns = nrow(QTLgeno)
    ss = apply(QTLgeno, 2, sum)
    maf.matrix = rbind((.5*ss/ns), (1-(0.5*ss/ns)))
    maf.temp = apply(maf.matrix, 2, min)
    length(which(maf.temp==0))/c(n.QTL.trait2*10)

    QTL.Loss.rate_temp=length(which(maf.temp==0))/c(n.QTL.trait2*10)

    SNP=pullSnpGeno(pop, snpChip = 1, simParam = SP)
    pheno=pheno(pop)
    pheno=data.frame(genotype=row.names(SNP),pheno)
    n.ind=nrow(pheno)
    #pheno=pheno[random.temp,]
    
    SNP.Map=getSnpMap(snpChip=1, simParam=SP)
    SNP.Map=SNP.Map[,-4]
    colnames(SNP.Map)=c("SNP.names" ,"chr" ,"pos")
    SNP.Map$allele1="A"
    SNP.Map$allele2="T"
    
    SNP.Map$chr=as.integer(SNP.Map$chr)
    SNP.Map$SNP.names=paste("X",SNP.Map$SNP.names,sep="")
    
    print(paste("Running Gwas at generation ",g,sep=""))
    gwas=sim.gwas(SNP=SNP,n.ind=pop@nInd,random.size=1000,p=gwas.threshold,distance=distance)

    hit.rate_temp=c(gwas/c(n.QTL.trait2*10))
  } 
#remove geno for memory 
    gc()
}# end of loop for(g in 0:n.gen){

#save information as recurrent population items (rec)
if(AF){
  ranMAF_temp= data.frame(r=r,maf)
}
else{
  ranMAF_temp = 0
}
if(PopStr){
  ranPopStr_temp=data.frame(r=r,BG.LD)
}
else{
  ranPopStr_temp = 0
}

list(rech2_t1=data.frame(run=r,gen=0:c(n.gen+1),h2_t1=Gvar_t1/(Pvar_t1+Gvar_t1)),
    recGvar_t1=data.frame(run=r,gen=0:c(n.gen+1),Gvar_t1),
    recPvar_t1=data.frame(run=r,gen=0:c(n.gen+1),Pvar_t1),
    recgenMean_t1=data.frame(run=r,gen=0:c(n.gen+1),genMean_t1),
    recphenMean_t1=data.frame(run=r,gen=0:c(n.gen+1),phenMean_t1),
    rech2_t2=data.frame(h2_t2=Gvar_t2/(Pvar_t2+Gvar_t2)),
    recGvar_t2=data.frame(Gvar_t2),
    recPvar_t2=data.frame(Pvar_t2),
    recgenMean_t2=data.frame(genMean_t2),
    recphenMean_t2=data.frame(phenMean_t2),
    ranMAF = ranMAF_temp,
    ranPopStr = ranPopStr_temp,
    QTL.Loss.rate = QTL.Loss.rate_temp,
    hit.rate = hit.rate_temp)

} ######## for(r in 1:runs){ ######
end.time=Sys.time()
print(end.time-start.time)

rech2_t1 = data[1, "rech2_t1"]
rech2_t2 = data[1, "rech2_t2"]
recGvar_t1 = data[1, "recGvar_t1"]
recGvar_t2 = data[1, "recGvar_t2"]
recPvar_t1 = data[1, "recPvar_t1"]
recPvar_t2 = data[1, "recPvar_t2"]
recgenMean_t1 = data[1, "recgenMean_t1"]
recgenMean_t2 = data[1, "recgenMean_t2"]
recphenMean_t1 = data[1, "recphenMean_t1"]
recphenMean_t2 = data[1, "recphenMean_t2"]
hit.rate = data[1, "hit.rate"]
QTL.Loss.rate = data[1, "QTL.Loss.rate"]
ranMAF = data[1, "ranMAF"]
ranPopStr = data[1, "ranPopStr"]
for(i in 2:runs)
{
  rech2_t1 = rbind(rech2_t1, data[i, "rech2_t1"])
  rech2_t2 = rbind(rech2_t2, data[i, "rech2_t2"])
  recGvar_t1 = rbind(recGvar_t1, data[i, "recGvar_t1"])
  recGvar_t2 = rbind(recGvar_t2, data[i, "recGvar_t2"])
  recPvar_t1 = rbind(recPvar_t1, data[i, "recPvar_t1"])
  recPvar_t2 = rbind(recPvar_t2, data[i, "recPvar_t2"])
  recgenMean_t1 = rbind(recgenMean_t1, data[i, "recgenMean_t1"])
  recgenMean_t2 = rbind(recgenMean_t2, data[i, "recgenMean_t2"])
  recphenMean_t1 = rbind(recphenMean_t1, data[i, "recphenMean_t1"])
  recphenMean_t2 = rbind(recphenMean_t2, data[i, "recphenMean_t2"])
  hit.rate = rbind(hit.rate, data[i, "hit.rate"])
  QTL.Loss.rate = rbind(QTL.Loss.rate, data[i, "QTL.Loss.rate"])
  ranMAF = rbind(ranMAF, data[i, "ranMAF"])
  ranPopStr = rbind(ranPopStr, data[i, "ranPopStr"])
}


    breeding_results=list(
    h2=cbind(rech2_t1, rech2_t2),
    Gvar=cbind(recGvar_t1, recGvar_t2),
    Pvar=cbind(recPvar_t1, recPvar_t2),
    genMean=cbind(recgenMean_t1, recgenMean_t2),
    phenMean=cbind(recphenMean_t1, recphenMean_t2)
  )
   if(run.gwas){
   r.hit.rate=hit.rate
  }
 r.QTL.Loss.rate=QTL.Loss.rate
 if(run.gwas){
  save(breeding_results,r.QTL.Loss.rate,r.hit.rate,file="RecurrentSelectionProgram.Rdata")
 } else {
   save(breeding_results,file="RecurrentSelectionProgram.Rdata")
 } 
 if(AF){
   write.table(ranMAF,file="RecurrentSelectionProgram.MAF.Results.txt",row.names = F,col.names = T,quote=F)
 }
 if(PopStr){
   write.table(ranPopStr,file="RecurrentSelectionProgram.PopStr.Results.txt",row.names = F,col.names = T,quote=F)
 }
#  rm(QTLgeno,QTLmap,plan.temp,crossPlan,founderPop,pop,BG.LD, data)
