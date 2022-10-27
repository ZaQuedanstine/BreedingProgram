#fuctions for simulation

#Run GWAS
#pop=pop
#SNP=SNP

#n.ind=5000
#random.size=1000

sim.gwas=function(SNP,n.ind,random.size,p,distance){

random.temp=sample(1:n.ind,size = random.size)
pheno.ran=pheno[random.temp,]
SNP.ran=SNP[random.temp,]


gDataDrops = createGData(geno = SNP.ran, map = SNP.Map, pheno = pheno.ran)
start.time=Sys.time()
GWASDrops = runSingleTraitGwas(gData = gDataDrops,
                               traits = c("Trait2"))
end.time=Sys.time()
end.time-start.time

results=as.data.frame(GWASDrops$GWAResult)
#remove markers with AF of 0 or 1
results=results[-c(which(results$pheno.ran.allFreq==0)),]
results=results[-c(which(results$pheno.ran.allFreq==1)),]
results$FDR=p.adjust(results$pheno.ran.pValue,method = "BH")
#marker passing FDR 0.05 threshold
results=results[which(results$FDR<=p),]

if(nrow(results)>0){
  for(i in 1:nrow(results)){
    #difference between hit and QTL
    temp.results=data.frame(
      chr=results$pheno.ran.chr[i],
      hit=i,
      diff=abs(QTLmap$site[which(QTLmap$chr==results$pheno.ran.chr[i])] - results$pheno.ran.pos[i]),
      if(abs(QTLmap$site[which(QTLmap$chr==results$pheno.ran.chr[i])] - results$pheno.ran.pos[i])<=distance){
        hit=TRUE} else {
          hit=FALSE
        }
    )      
    colnames(temp.results)[4]="success"

    if(i==1){     
      hit.results=temp.results
    } else{
      hit.results=rbind(hit.results,temp.results)
    }
    rm(temp.results)
    colnames(hit.results)[4]="success"
  } # for(i in 1:nrow(gwas)){
  #record T rate
  tmp=hit.results[which(hit.results$success==TRUE),]
  return(length(unique(tmp$chr)))
} else{ #if(nrow(gwas)>0){
  return(0)
} #else{ #if(nrow(gwas)>0){

#return(results, length(unique(tmp$chr)))
       
}

#add in option to use kinship and rrblup
###its slower since it gets kinship each time but may be prefered later
# SNP=as.data.frame(pullSnpGeno(pop, snpChip = 1, simParam = SP))
# pheno=pheno(pop)
# pheno=data.frame(genotype=row.names(SNP),pheno)
# 
# SNP=SNP[random.temp,]
# SNP.Map=getSnpMap(snpChip=1, simParam=SP)
# SNP.Map$id=paste("X",SNP.Map$id,sep="")
# SNP=data.frame(SNP.Map[,-4],t(SNP))
# View(SNP[1:100,1:100])
# 
# pheno$genotype=paste("X",pheno$genotype,sep="")
# pheno=pheno[random.temp,]
# 
# start.time=Sys.time()
# library(rrBLUP)
# GWAS.rrblup=GWAS(pheno[,c(1,3)], SNP, fixed=NULL, K=NULL, n.PC=0,
#                  min.MAF=0.05, n.core=1, P3D=TRUE, plot=F)
# end.time=Sys.time()
# end.time-start.time
# 
# results=as.data.frame(GWAS.rrblup)
# #remove markers with AF of 0 or 1
# results=results[-c(which(results$pheno.allFreq==0)),]
# results=results[-c(which(results$pheno.allFreq==1)),]
# results$FDR=p.adjust(results$pheno.pValue,method = "BH")

########compute background LD from HBP GBS data########
#take 25 random SNPs from each chromosome


backgroundLD=function(map,geno,rep,n.SNP){
for(j in 1:rep){ #run random loop j times
  if(j==1){
    p50=c()
    p90=c()
    p95=c() 
  }
  for(i in 1:10){
    #print(i)
    if(i==1){
      temp.geno=map[which(map$LG==1),]
      randomSNPs=temp.geno[sample(1:nrow(temp.geno),n.SNP,replace=F),]
    }else {
      temp.geno=map[which(map$LG==i),]
      randomSNPs.temp=temp.geno[sample(1:nrow(temp.geno),n.SNP,replace=F),]
      randomSNPs=rbind(randomSNPs,randomSNPs.temp)
    }
  }
  #compute pairwise LD
  geno_random=geno[,randomSNPs$Locus]
  
  LDdecay=LD.decay(geno_random,randomSNPs,silent=T,unlinked=F)
  
  LD=LDdecay$all.LG
  if(nrow(LD)>1){
  LD=LD[-c(which(LD$d==0)),]
  #the 95th percentile value is the estimate of background LD
  if(length(which(is.na(LD)==T))>0){
    LD=na.omit(LD)
  }
   p50=c(p50,quantile(LD$r2,probs = 0.50))
  p90=c(p90,quantile(LD$r2,probs = 0.90))
  p95=c(p95,quantile(LD$r2,probs = 0.95)) 
  } else{
    p50=c(p50,NA)
    p90=c(p90,NA)
    p95=c(p95,NA) 
  }
  
}
return(data.frame(percentile=c(50,90,95),mean=c(mean(p50),mean(p90),mean(p95))))
} #end of function


#getting minor allele freq
MAF=function(geno,return){
####Obtain the mafs of all SNPs###
#Total number of lines
ns = nrow(geno)
#Sum of the allele scores for each SNP
ss = apply(geno, 2, sum)
#Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
maf.matrix = rbind((.5*ss/ns), (1-(0.5*ss/ns)))
#Copy the minor allele frequencies for all SNPs
maf = apply(maf.matrix, 2, min)
maf=as.data.frame(maf)
row.names(maf)=colnames(geno)
if(return=="All"){
  return(maf)
} else {
  return(data.frame(
    mean=mean(maf[,1]),
     median= median(maf[,1]),
      sd=sd(maf[,1]))  
  )
}
}



