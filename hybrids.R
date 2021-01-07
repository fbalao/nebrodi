
library(vcfR)
library(hierfstat)
library(poppr)
library(ggplot2)
library(factoextra)

#import vcf into geneind object

vcf_abies <- read.vcfR("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/hybrids/R0.3populations_TODOsnps.recode.vcf", verbose = FALSE )
markers<-read.table("Abies_nebrodi/openarray_data/hybrids/markers.recode.txt", header=T)

abies_genind <- vcfR2genind(vcf_abies, return.alleles=T)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
abies_genind@pop<-as.factor(pop)
abiesnames<-read.table("namesabiespecies.txt", sep=",")
indNames(abies_genind)<-abiesnames[match(indNames(abies_genind), abiesnames$V1), 2]
locNames(abies_genind)<-markers[match(locNames(abies_genind), markers$ID), "loci"]

abies_genindc<-abies_genind[6:35]

x <-indpca(abies_genindc) 
plot(x, cex =1.2, col=rainbow(30))

#### extraer las plantulas sospechosas de ser hibridas

indNames(plantulasNDB)<-padresnames[match(indNames(plantulasNDB), padresnames$Codigo_muestra), "ID_LIFE"]
#indNames(padresDB)<-padresnames[match(indNames(padresDB), padresnames$Codigo_muestra), "ID_LIFE"]

plantulas_sospechosas<-c("16.2.P","16.3.P","18.10.P","18.13.P","20.1.P",
                         "21.3.P","21.4.P","21.5.P","22.18.P")

plantulas_sosDB<-plantulasNDB[plantulas_sospechosas]


xloci<-intersect(locNames(padresDB), locNames(abies_genindc))



abies_genindc_sub<-abies_genindc[loc = xloci]
padresDB_sub<-padresDB[loc= xloci]
padresDB_sub@pop<-as.factor(rep("adults",30))

plantulas_sosDB_sun<-plantulas_sosDB[loc= xloci]
plantulas_sosDB_sun@pop<-as.factor(rep("seedlings",9))

poolabies<-repool(abies_genindc_sub, padresDB_sub,plantulas_sosDB_sun)


y <-indpca(poolabies, ind.labels=indNames(poolabies)) 
plot(y, cex =0.8, col=funky(6)[1:4][poolabies@pop])

abline(v=0, lty=2, col="grey75", lwd=0.5)
abline(h=0, lty=2, col="grey75", lwd=0.5)

legend(-0.6,-1,legend = c("A. alba", "A. cephalonica", 
                          "A. nebrodensis adults", "A. nebrodensis seedlings"), 
       pch=16, col=funky(6)[c(2,1,3,4)], cex=0.9, 
       box.lwd = 0, border = "white")

nebrotab<-tab(poolabies, freq=T, NA.method="mean")
pca.neb <- dudi.pca(nebrotab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.neb$li, fac=pop(poolabies),
        col=rainbow(30), clabel = 0.9)

rownames(pca.neb$li)<-padresDB@other$Madre_procedencia

colx<-rainbow(4)[pop(poolabies)]
labels = c("A. nebrodensis adults", "A. nebrodensis seedlings", "A. alba", "A. cephalonica")

fviz_pca_ind(pca.neb,
             col.ind = rainbow(4)[pop(poolabies)],
             show.legend.text = FALSE,
             labelsize=4,
             addEllipses=F,
             label="all",# Color by the quality of representation
           #  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title="",
           legend.title = "Species",
           legend.labels=c("1","2","a","b"),             repel = TRUE) +
 scale_color_manual(name = "Species", labels = c("A. nebrodensis adults", "A. nebrodensis seedlings", "A. alba", "A. cephalonica"),
                   values= c("red","blue","orange","forestgreen"))




kneb<-find.clusters(poolabies)
neb_dapc<-dapc(poolabies,pop = kneb$grp, n.pca = 4,n.da = 2)   
compoplot(neb_dapc,border=T,
          cex.names=.61,    
          col=funky(6),show.lab=T, legend=F)
   temp <- optim.a.score(neb_dapc)   




###############3 subset loci hybrids

abies_genind <- vcfR2genind(vcf_abies, return.alleles=T)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
abies_genind@pop<-as.factor(pop)
indNames(abies_genind)<-substr(indNames(abies_genind),21,32)

introgressmarkers<-read.table("Abies_nebrodi/openarray_data/hybrids/introgressmarkers.txt", header=F)
abiesintrogress_genind<-abies_genind[loc=introgressmarkers$V1]

locNames(abiesintrogress_genind)<-markers[match(locNames(abiesintrogress_genind), markers$ID), "loci"]

xloci<-intersect(locNames(padresDB), locNames(abiesintrogress_genind))


abiesintrogress_genind_sub<-abiesintrogress_genind[loc = xloci]
padresDB_sub<-padresDB[loc= xloci]
padresDB_sub@pop<-as.factor(rep("adults",30))

plantulas_sosDB_sun<-plantulas_sosDB[loc= xloci]
plantulas_sosDB_sun@pop<-as.factor(rep("seedlings",9))

poolabies2<-repool(abiesintrogress_genind_sub, padresDB_sub,plantulas_sosDB_sun)

y <-indpca(poolabies2, ind.labels=indNames(poolabies2)) 
plot(y, cex =0.8, col=rainbow(5)[poolabies2@pop])
   legend(1.1,-1.1,legend = c("A. alba", "A. cephalonica", "A. nebrodensis old",
                        "A. nebrodensis adults", "A. nebrodensis seedlings"), 
       pch=16, col=rainbow(5)[c(3,2,1,4,5)], cex=0.8, 
       box.lwd = 0, border = "white")


   kneb<-find.clusters(poolabies2)
   neb_dapc<-dapc(poolabies2,pop = kneb$grp, n.pca = 4,n.da = 3)   
   compoplot(neb_dapc,border=T,
             cex.names=.7
             ,
             col=funky(6),show.lab=T, legend=F)
   temp <- optim.a.score(neb_dapc)   
   
   compoplot(neb_dapc,border=T,
             cex.names=.61,    
             col=funky(6),show.lab=T, legend=F)
   
genind2structure(poolabies, "poolabies.str")


fs<-read.table("poolabies_fasts.2.meanQ", sep=",")

barplot(t(as.matrix(fs[,2:3])),beside=F, names.arg=indNames(poolabies2), las=2, cex.names=0.6, col=funky(6))



fs2<-read.table("poolabies2_fasts.2.meanQ", sep=",")

barplot(t(as.matrix(fs2[,2:3])),beside=F, names.arg=indNames(poolabies2), las=2, cex.names=0.6, col=funky(6))
   





layout(matrix(c(1,2),ncol=1))
par(oma=c(3.5,0.5,0,1))
par(mar=c(4,4,0.5,0.5) + 0.1)
y <-indpca(poolabies, ind.labels=indNames(poolabies)) 
plot(y, cex =0.8, col=funky(6)[1:4][poolabies@pop])

abline(v=0, lty=2, col="grey75", lwd=0.5)
abline(h=0, lty=2, col="grey75", lwd=0.5)

legend(-0.40,-0.30,legend = c("A. alba", "A. cephalonica", 
                          "A. nebrodensis adults", "A. nebrodensis seedlings"), 
        pch=16, col=funky(6)[c(2,1,3,4)], cex=0.7, 
       box.lwd = 0, border = "white")  

par(mar=c(2,4,2,0.5) + 0.1)

compoplot(neb_dapc,border=T,
          cex.names=.61,    
          col=funky(6),show.lab=T, legend=F)

