library(pedantics)

nebroped<-read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/colony/nebrodi_results3/nebrodiParesdePadresR.txt", header=T)
colnames(nebroped)<-c("id","father","mother", "prob")

nebroped118<-read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/colony/nebrodi118plantulas/nebrodi118_results/nebrodi118ParesdePadres.txt", 
                        na.strings="#", sep=",",header=T)
colnames(nebroped118)<-c("id","father","mother", "prob")
nebroped118<-nebroped118[nebroped118$prob>0.4,]

comPar<-PedCompare(nebroped,nebroped118)

comPar$Counts["GX",, ]

nebroped2<-GeneticsPed:::extend(nebroped)


nebroped2<-nebroped2[,c(1,3,2)]
colnames(nebroped2)<-c("id","dam","sire")

nebroped3<-prepPed(nebroped2)
pedigreeStats(nebroped2, lowMem=T)


unkfather<-sum(is.na(nebroped$father))/length(nebroped$id)
unkmother<-sum(is.na(nebroped$mother))/length(nebroped$id)
