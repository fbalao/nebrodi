# libraries ----
library(vcfR)
library(hierfstat)
library(ggplot2)
library(factoextra)
library(here)
library(RColorBrewer)


load(here("envs","nebrodensis.21.12.2020.RData"))



# PCA especies, adultos y plantulas sospechosas ----
# subset loci hibridacion
# Se necesita PadresDB, plantulas_sosDB


vcf_abies <- read.vcfR("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/hybrids/R0.3populations_TODOsnps.recode.vcf", verbose = FALSE )
markers<-read.table("Abies_nebrodi/openarray_data/hybrids/markers.recode.txt", header=T)

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

padresnames<-read.table("Abies_nebrodi/openarray_data/diversity/muestras2.csv", header=T)
indNames(plantulasNDB)<-padresnames[match(indNames(plantulasNDB), padresnames$Codigo_muestra), "ID_LIFE"]
plantulas_sospechosas<-c("16.2.P","16.3.P","18.10.P","18.13.P","20.1.P",
                         "21.3.P","21.4.P","21.5.P","22.18.P")
plantulas_sosDB<-plantulasNDB[plantulas_sospechosas]


plantulas_sosDB_sun<-plantulas_sosDB[loc= xloci]
plantulas_sosDB_sun@pop<-as.factor(rep("seedlings",9))

poolabies2<-repool(abiesintrogress_genind_sub, padresDB_sub,plantulas_sosDB_sun)


col1=brewer.pal(5,"Spectral")[c(3,2,4,1)]
col1[4]<-"#e58f90"

nebrotab2<-tab(poolabies2, freq=T, NA.method="mean")
pca.neb2 <- dudi.pca(nebrotab2, center = T, scale = T, scannf = F, nf = 2)

summary(pca.neb2)

x<-pca.neb2$li$Axis1<-pca.neb2$li$Axis1*(-1)
y<-pca.neb2$li$Axis2<-pca.neb2$li$Axis2*(-1)


circ<-rep("grey65",69)
circ[c(39,59,65)]<-"black"
col3<-col1
bgcol<-col1[pop(poolabies2)]

#Figure 3
plot(x,-y, bg=bgcol, pch=21, cex=1.7, col=circ,
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


