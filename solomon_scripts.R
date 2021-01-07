load("nebrodensis.21.12.2020.RData")

library(dplyr)
library(outliers)
library(strataG)
library(pegas)

#Preparación del data set
padresDB<-nebrodf3[nebrodf3@other$Plantula_madre==0 & nebrodf3@other$Vivero_campo==0]
padresDB_loci<-genind2loci(padresDB)
summary(padresDB_loci)
#tabla de frecuencia
# af<-rraf(padresDB,res = "data.frame")
# af$frequency<-round(af$frequency,3)
#makefreq(as.genpop(cleaned_mothers2))
#write.table(af[,-3], "/home/fbalao/solomon/motherAF.txt")
#write.loci(padresDB_loci, loci.sep = "\t", quote = FALSE, col.names = T, file="PADRES_nebro.txt")


padresDB_solomon<-as.data.frame(padresDB_loci)
padresDB_solomon2<-alleleSplit(padresDB_solomon[,-1], sep = "/")
rownames(padresDB_solomon2)<-rownames(padresDB_solomon)


plantulasNDB<-nebrodf3[nebrodf3@other$Plantula_madre==1 & nebrodf3@other$Vivero_campo==0]
plantulasNDB_loci<-genind2loci(plantulasNDB)
summary(plantulasNDB_loci)
#tabla de frecuencia
# af<-rraf(plantulasNDB,res = "data.frame")
# af$frequency<-round(af$frequency,3)
#makefreq(as.genpop(cleaned_mothers2))
#write.table(af[,-3], "/home/fbalao/solomon/motherAF.txt")
#write.loci(plantulasNDB_loci, loci.sep = "\t", quote = FALSE, col.names = T, file="plantulasN_nebro.txt")


plantulasNDB_solomon<-as.data.frame(plantulasNDB_loci)
plantulasNDB_solomon2<-alleleSplit(plantulasNDB_solomon[,-1], sep = "/")
rownames(plantulasNDB_solomon2)<-rownames(plantulasNDB_solomon)


plantulasNDB_solomon2[is.na(plantulasNDB_solomon2)] <- 0
padresDB_solomon2[is.na(padresDB_solomon2)] <- 0
changeNucleotides2numbers<-function(juve){
juve[juve=="A"]<-'100'
juve[juve=="T"]<-'101'
juve[juve=="C"]<-'102'
juve[juve=="G"]<-'103'
juve[is.na(juve)] <-'0'
juve
}
padresDB_solomon3<-changeNucleotides2numbers(padresDB_solomon2)
plantulasNDB_solomon3<-changeNucleotides2numbers(plantulasNDB_solomon2)


write.table(padresDB_solomon3,"./padresDB_solomon.txt",
            sep="\t",quote = T)
write.table(plantulasNDB_solomon3,"./plantulasNDB_solomon.txt",
            sep="\t",quote=T)

library(SOLOMON)
solomon()

# SOLOMON mother known
df2<-df[rep(seq_len(nrow(df)), each = 10), ]
write.table(df2, "/home/fbalao/solomon/MomsKnown.txt", row.names = F)

#SOLOMON NA add
juve<-read.table("/home/fbalao/solomon/Juveniles_sim.txt", header=T)
juve2 <- apply (juve[,-1], 2, function(x) {x[sample( c(1:300), floor(300/10))] <- NA; x} )
juve3<-cbind(juve[,1], juve2)
f<-function(x){sum(is.na(x))}
mean(apply(juve3,1, f)/204)
write.table(juve3, "/home/fbalao/solomon/Juveniles_0.1NA.txt", row.names = F, sep="\t", quote = F, na="0")


juve[juve=='100']<-'A'
juve[juve=='101']<-'T'
juve


f.jazurro <- function(mydf) {
  odd <- seq(1, ncol(mydf), 2);
  lapply(odd, function(x) paste(mydf[,x], mydf[,x+1], sep = "/")) %>% 
    do.call(cbind,.)
}

juveL<-f.jazurro(juve[,-1])
juveL<-cbind(juve$IDs, juveL)
loci<-paste0("loci",01:102)
loci<-c("Genotype", loci)
colnames(juveL)<-loci
juveL<-as.data.frame(juveL)
juveL$key<-rep("Off",300)
juveL<-juveL[,c(1,104,2:103)]


parents<-read.table("/home/fbalao/solomon/todo_sim.txt", header=T)

parents[parents=='100']<-'A'
parents[parents=='101']<-'T'

parentsL<-f.jazurro(parents[,-1])
parentsL<-cbind(parents$IDs, parentsL)
colnames(parentsL)<-loci
parentsL<-as.data.frame(parentsL)
parentsL$key<-c(rep("Pa", 30), rep("Fa",30))
parentsL<-parentsL[,c(1,104,2:103)]

todo<-rbind(parentsL, juveL)
library(outliers)
source("apparent.R")
apparentOUT <- apparent(todo, MaxIdent=0.10, alpha=0.10, nloci=80, self=F, plot=TRUE, Dyad=FALSE)
        