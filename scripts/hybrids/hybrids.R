# libraries ----
library(vcfR)
library(hierfstat)
library(poppr)
library(ggplot2)
library(factoextra)
library(here)
library(RColorBrewer)
library(GISTools)

#import vcf into geneind object ----
load(here("envs","nebrodensis.21.12.2020.RData"))

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






compoplot(neb_dapc,border=T,
          cex.names=.61,    
          col=funky(6),show.lab=T, legend=F)


# PCA especies, adultos y plantulas sospechosas ----
# subset loci hibridacion

abies_genind <- vcfR2genind(vcf_abies, return.alleles=T)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
abies_genind@pop<-as.factor(pop)
indNames(abies_genind)<-substr(indNames(abies_genind),21,32)
abies_genind<-abies_genind[-c(1:5)] # quitar nebrodensis_old

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

y2 <-indpca(poolabies2, ind.labels=indNames(poolabies2)) 
plot(y2, cex =0.8, col=rainbow(5)[poolabies2@pop])
legend(10,-5,legend = c("A. alba", "A. cephalonica", "A. nebrodensis old",
                           "A. nebrodensis adults", "A. nebrodensis seedlings"), 
       pch=16, col=rainbow(5)[c(3,2,1,4,5)], cex=0.8, 
       box.lwd = 0, border = "white")



col1=brewer.pal(5,"Spectral")[c(3,2,4,1)]

col1[4]<-"#e58f90"

plot(1:4, cex=4, col=col1,pch=16)

nebrotab2<-tab(poolabies2, freq=T, NA.method="mean")
pca.neb2 <- dudi.pca(nebrotab2, center = T, scale = T, scannf = F, nf = 2)

summary(pca.neb2)

x<-pca.neb2$li$Axis1<-pca.neb2$li$Axis1*(-1)
y<-pca.neb2$li$Axis2<-pca.neb2$li$Axis2*(-1)


circ<-rep("grey65",69)
circ[c(39,59,65)]<-"black"

col3<-col1



bgcol<-col1[pop(poolabies2)]
#bgcol[c(39,59,65)]<-c(col2[1],col2[1], col2[2])

#Figure 3
plot(x,-y, bg=bgcol, pch=21, cex=1.7, ylim=c(-28,28), col=circ,
     xlab="PC1 (28.6%)", ylab="PC2 (13.4%)", cex.axis=1.3, cex.lab=1.4)

legend(10,-15,legend = c(expression(italic("A. alba")),
                      expression(italic("A. cephalonica")),
                      expression(paste(italic("A. nebrodensis")," adults")),
                      expression(paste(italic("A. nebrodensis"), " seedlings"))),
       pch=21, cex=1.2, col="grey35", pt.bg=col1, pt.cex=1.5,
       box.lwd = 0, border = "white", bg="transparent")

legend(11,-15,legend = c("Silver fir", "Greek fir", "Sicilian fir adults", "Sicilian fir seedlings"), 
       pch=21, cex=1.2, col="grey35", pt.bg=col1, pt.cex=1.5,
       box.lwd = 0, border = "white", bg="transparent")

#text(x,y,labels = rownames(pca.neb2$li))

s.class(pca.neb2$li, fac=pop(poolabies2),
        col=col1, clabel = 0, cstar=0, cellipse = 0,cpoint = 2, contour=T) 

rownames(pca.neb$li)<-padresDB@other$Madre_procedencia

colx<-rainbow(4)[pop(poolabies2)]
labels = c("A. nebrodensis adults", "A. nebrodensis seedlings", "A. alba", "A. cephalonica")

fviz_pca_ind(pca.neb2,
             col.ind = "black", 
             fill.ind=rainbow(5)[pop(poolabies2)],
             show.legend.text = FALSE,
             labelsize=0,
             addEllipses=F,
             geom = c("point", "text"),
             alpha.ind = 0.5,
             label="all",
             title="",
             legend.title = "Species",
             legend.labels=c("1","2","a","b","z"),repel = TRUE)

+
  scale_color_manual(name = "Species", labels = c("A. nebrodensis adults", "A. nebrodensis seedlings", "A. alba", "A. cephalonica"),
                     values= c("red","blue","orange","forestgreen"))


