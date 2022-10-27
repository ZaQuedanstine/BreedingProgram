#plotting simulation results
library(ggplot2)
library(hrbrthemes)

##########plot GWAS hit rate############
####plotting breeding/variance parameters############
load("BiparentalSelectionProgram.Rdata")
BiparentalBreedingResults=breeding_results
rm(breeding_results)
b.hit.rate=data.frame(Biparental=b.hit.rate)
b.QTL.Loss.rate=data.frame(Biparental=b.QTL.Loss.rate)

load("RecurrentSelectionProgram.Rdata")
RecurrentBreedingResults=breeding_results
rm(breeding_results)

r.hit.rate=data.frame(Biparental=r.hit.rate)
r.QTL.Loss.rate=data.frame(Biparental=r.QTL.Loss.rate)

QTL.Loss.rate=cbind(b.QTL.Loss.rate,r.QTL.Loss.rate)
colnames(QTL.Loss.rate)=c("Biparental","Recurent")

hit.rate=cbind(b.hit.rate,r.hit.rate)
colnames(hit.rate)=c("Biparental","Recurent")
rm(r.hit.rate,b.hit.rate,r.QTL.Loss.rate,b.QTL.Loss.rate)

###### need to edit so y axis is absolute values
hit.rate.plot = ggplot(hit.rate, aes(x=x) ) +
  # Top
  geom_density( aes(x = Recurent, y = ..density..), fill="#00bfc4" ) +
  geom_label( aes(x=.5, y=2.5, label="Random Mating"), color="#00bfc4") +
  # Bottom
  geom_density( aes(x = Biparental, y = -..density..), fill= "#E69F00") +
  geom_label( aes(x=.5, y=-2.5, label="Pairwise"), color="#E69F00") +
  theme_ipsum() +
  xlab("Proportion of QTL Detected")+
  ylab("Density")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

pdf(file="GWAShitrate.pdf",useDingbats = F)
print(hit.rate.plot)
dev.off()

jpeg(file="GWAShitrate.jpeg")
print(hit.rate.plot)
dev.off()



###### need to edit so y axis is absolute values
##### can we make distributions smoothed?- look normal?
QTL.Loss.rate.plot = ggplot(QTL.Loss.rate, aes(x=x) ) +
  # Top
  geom_density( aes(x = Recurent, y = ..density..), fill="#00bfc4" ) +
  geom_label( aes(x=.25, y=4, label="Random Mating"), color="#00bfc4") +
  # Bottom
  geom_density( aes(x = Biparental, y = -..density..), fill= "#E69F00") +
  geom_label( aes(x=.5, y=-2.5, label="Pairwise"), color="#E69F00") +
  theme_ipsum() +
  xlab("Proportion of QTL Fixed")+
  ylab("Density")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

pdf(file="QTLLossrate.pdf",useDingbats = F)
print(QTL.Loss.rate.plot)
dev.off()

jpeg(file="QTLLossrate.jpeg")
print(QTL.Loss.rate.plot)
dev.off()


#plot genetic variance
plot=ggplot() + 
  geom_smooth(method = "loess",data=RecurrentBreedingResults$Gvar,aes(x=gen,y=Gvar_t1),se=T,color="#00bfc4")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$Gvar,aes(x=gen,y=Gvar_t1),se=T,color="#E69F00")+
  geom_smooth(method = "loess",data=RecurrentBreedingResults$Gvar,aes(x=gen,y=Gvar_t2),se=T,color="#00bfc4",linetype = "dashed")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$Gvar,aes(x=gen,y=Gvar_t2),se=T,color="#E69F00",linetype = "dashed")+
  labs(x="Generation",y="Genetic Variance")+
  theme_classic()+theme(text = element_text(size = 20))


pdf(file="GeneticVarianceBiParentalVsRecurrant.pdf",useDingbats = F)
print(plot)
dev.off()

jpeg(file="GeneticVarianceBiParentalVsRecurrant.jpeg")
print(plot)
dev.off()

#Genetic Gain
#plot Genetic mean
plot=ggplot() + 
  geom_smooth(method = "loess",data=RecurrentBreedingResults$genMean,aes(x=gen,y=genMean_t1),se=T,color="#00bfc4")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$genMean,aes(x=gen,y=genMean_t1),se=T,color="#E69F00")+
  geom_smooth(method = "loess",data=RecurrentBreedingResults$genMean,aes(x=gen,y=genMean_t2),se=T,color="#00bfc4",linetype = "dashed")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$genMean,aes(x=gen,y=genMean_t2),se=T,color="#E69F00",linetype = "dashed")+
  labs(x="Generation",y="Genetic Mean")+
  theme_classic()+theme(text = element_text(size = 20))


pdf(file="GeneticMeanBiParentalVsRecurrant.pdf",useDingbats = F)
print(plot)
dev.off()

jpeg(file="GeneticMeanBiParentalVsRecurrant.jpeg")
print(plot)
dev.off()

#plot Phenotypic mean
plot=ggplot() + 
  geom_smooth(method = "loess",data=RecurrentBreedingResults$phenMean,aes(x=gen,y=phenMean_t1),se=T,color="#00bfc4")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$phenMean,aes(x=gen,y=phenMean_t1),se=T,color="#E69F00")+
  geom_smooth(method = "loess",data=RecurrentBreedingResults$phenMean,aes(x=gen,y=phenMean_t2),se=T,color="#00bfc4",linetype = "dashed")+ 
  geom_smooth(method = "loess",data=BiparentalBreedingResults$phenMean,aes(x=gen,y=phenMean_t2),se=T,color="#E69F00",linetype = "dashed")+
  labs(x="Generation",y="Phenotypic Mean")+
  theme_classic()+theme(text = element_text(size = 20))


pdf(file="PhenoMeanBiParentalVsRecurrant.pdf",useDingbats = F)
print(plot)
dev.off()

jpeg(file="PhenoMeanBiParentalVsRecurrant.jpeg")
print(plot)
dev.off()



bi.maf=read.table("BiparentalSelectionProgram.MAF.Results.txt",header=T)
ran.maf=read.table("RecurrentSelectionProgram.MAF.Results.txt",header=T)
#MAF mean
plot=ggplot() + 
  geom_smooth(method = "loess",data=ran.maf,aes(x=g,y=mean),se=T,color="#00bfc4")+ 
  geom_smooth(method = "loess",data=bi.maf,aes(x=g,y=mean),se=T,color="#E69F00")+
  labs(x="Generation",y="Mean MAF")+
  theme_classic()+
    theme(text = element_text(size = 20))

pdf(file="MAFBiParentalVsRecurrant.pdf",useDingbats = F)
print(plot)
dev.off()

jpeg(file="MAFBiParentalVsRecurrant.jpeg")
print(plot)
dev.off()

bi.popstr=read.table("BiparentalSelectionProgram.PopStr.Results.txt",header=T)
bi.popstr=bi.popstr[which(bi.popstr$percentile=="95"),]
  ran.popstr=read.table("RecurrentSelectionProgram.PopStr.Results.txt",header=T)
  ran.popstr=ran.popstr[which(ran.popstr$percentile=="95"),]
  
#Plot LD Decay
plot=ggplot() + 
  geom_smooth(method = "loess",data=ran.popstr,aes(x=gen,y=mean),se=T,color="#00bfc4")+ 
  geom_smooth(method = "loess",data=bi.popstr,aes(x=gen,y=mean),se=T,color="#E69F00")+
  labs(x="Generation",y="Mean LD (95th)")+
  theme_classic()+
  theme(text = element_text(size = 20))

pdf(file="PopSt950pParentalVsRecurrant.pdf",useDingbats = F)
print(plot)
dev.off()

jpeg(file="PopStr95pParentalVsRecurrant.jpeg")
print(plot)
dev.off()



#for each chromosome get a matrix of Minor Allele around QTL
trait=2
QTL=getQtlMap(trait = trait,simParam=SP)

apply(biMAF,2,mean)

apply(recMAF,2,mean)
chr=1
pos=1
window=10
temp=biMAF[paste(chr,"_",1:1000,sep=""),]

temp=as.data.frame(cbind(1:1000,temp))
temp=temp[c(c(which(rownames(temp)==QTL$id[which(QTL$chr==chr)][pos])-window):
              c(which(rownames(temp)==QTL$id[which(QTL$chr==chr)][pos])+window)),
]

temp2=recMAF[paste(chr,"_",1:1000,sep=""),]
temp2=as.data.frame(cbind(1:1000,temp2))
temp2=temp2[c(c(which(rownames(temp2)==QTL$id[which(QTL$chr==chr)][pos])-window):
                c(which(rownames(temp2)==QTL$id[which(QTL$chr==chr)][pos])+window)),
]


if(trait==1){
  plot=ggplot() + 
    geom_smooth(method = "loess",data=temp2,aes(x=V1,y=V2),se=FALSE,color="#00bfc4")+
    geom_vline(xintercept=QTL$site[which(QTL$chr==chr)[pos]],size = 2, linetype="dashed", color = "#f8766d")+
    geom_smooth(method = "loess",data=temp,aes(x=V1,y=V2),se=FALSE,color="#E69F00")+
    labs(x="Position",y="MAF")+
    theme_classic()+theme(text = element_text(size = 20))
  
  pdf(file="MAFaroundQTLBiParentalVsRecurrant_SelectedTrait.pdf",useDingbats = F)
  print(plot)
  dev.off()
  
  jpeg(file="MAFaroundQTLBiParentalVsRecurrant_SelectedTrait.jpeg")
  print(plot)
  dev.off()
}
if(trait==2){
  plot=ggplot() + 
    geom_smooth(method = "loess",data=temp2,aes(x=V1,y=V2),se=FALSE,color="#00bfc4",linetype="dashed")+
    geom_vline(xintercept=QTL$site[which(QTL$chr==chr)[pos]],size = 2, linetype="dashed", color = "#f8766d")+
    geom_smooth(method = "loess",data=temp,aes(x=V1,y=V2),se=FALSE,color="#E69F00",linetype="dashed")+
    labs(x="Position",y="MAF")+
    theme_classic()+theme(text = element_text(size = 20))
  
  pdf(file="MAFaroundQTLBiParentalVsRecurrant_nonSelectedTrait.pdf",useDingbats = F)
  print(plot)
  dev.off()
  
  jpeg(file="MAFaroundQTLBiParentalVsRecurrant_nonSelectedTrait.jpeg")
  print(plot)
  dev.off()
}









