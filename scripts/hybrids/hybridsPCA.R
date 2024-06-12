# libraries ----
library(vcfR)
library(hierfstat)
library(ggplot2)
library(factoextra)
library(here)
library(RColorBrewer)
library(adegenet)

# we need nebrodf3_118 from AGP
load(here("envs","AGP_results8.11.2023.RData"))



# PCA especies, adultos y plantulas sospechosas ----
# subset loci hibridacion
# Se necesita PadresDB, plantulas_sosDB


vcf_abies <- read.vcfR(here("data","openarray_design","data","R0.3populations_TODOsnps.recode.vcf"), verbose = FALSE )
markers<-read.table(here("data","openarray_design","data","markers.recode.txt"), header=T)

abies_genind <- vcfR2genind(vcf_abies, return.alleles=T)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
abies_genind@pop<-as.factor(pop)
indNames(abies_genind)<-substr(indNames(abies_genind),21,32)
abies_genind<-abies_genind[-c(1:5)] # quitar nebrodensis_old

introgressmarkers<-read.table(here("data","openarray_design","data","introgressmarkers.txt"), header=F)
abiesintrogress_genind<-abies_genind[loc=introgressmarkers$V1]



locNames(abiesintrogress_genind)<-markers[match(locNames(abiesintrogress_genind), markers$ID), "loci"]

nebrodf3_118introgress_genind<-nebrodf3_118[loc=locNames(abiesintrogress_genind)]

# # Comparar los nombres de alelos
# all_names_1 <-abiesintrogress_genind@all.names
# all_names_2 <- nebrodf3_118introgress_genind@all.names
# comparison <- mapply(function(x, y) identical(x, y), all_names_1, all_names_2)
# 
# all_names_1["a00064811_27192"]
# all_names_2["a00064811_27192"]
# 
# test<-repool(abiesintrogress_genind,nebrodf3_118introgress_genind)
# 
# nebrodf3_118introgress_genind[loc="a00064811_27192"]@tab[1,]
# abiesintrogress_genind[loc="a00064811_27192"]@tab[1,]
# test[loc="a00064811_27192"]@tab[c("01_2019_2121","_cep_39_2_sb"),]




xloci<-intersect(locNames(nebrodf3_118introgress_genind), locNames(abiesintrogress_genind))
#21 loci



padresDB_sub<-nebrodf3_118introgress_genind[nebrodf3_118introgress_genind@other$Vivero_campo==0 & nebrodf3_118introgress_genind@other$Plantula_madre==0]
padresnames<-read.table(here("data","openarray_design","data","muestras2.csv"), header=T)
padresDB_sub@pop<-as.factor(rep("adults",30))
indNames(padresDB_sub)<-padresnames[match(indNames(padresDB_sub), padresnames$Codigo_muestra), "ID_LIFE"]


plantulasNDB_sub<-nebrodf3_118introgress_genind[nebrodf3_118introgress_genind@other$Vivero_campo==0 & nebrodf3_118introgress_genind@other$Plantula_madre==1]
indNames(plantulasNDB_sub)<-padresnames[match(indNames(plantulasNDB_sub), padresnames$Codigo_muestra), "ID_LIFE"]
plantulas_sospechosas<-c("16.2.P","16.3.P","18.10.P","18.13.P","20.1.P",
                         "21.3.P","21.4.P","21.5.P","22.18.P")
plantulas_sosDB<-plantulasNDB_sub[plantulas_sospechosas]
plantulas_sosDB@pop<-as.factor(rep("seedlings",9))

poolabies2<-repool(abiesintrogress_genind, padresDB_sub,plantulas_sosDB)


col1=brewer.pal(5,"Spectral")[c(3,2,4,1)]
col1[4]<-"#e58f90"

col1<-c("#FDAE61","#FFFFBF", "#ABDDA4", "#e58f90")

nebrotab2<-tab(poolabies2, freq=T, NA.method="mean")
pca.neb2 <- dudi.pca(nebrotab2, center = T, scale = T, scannf = F, nf = 2)

summary(pca.neb2)

x<-pca.neb2$li$Axis1<-pca.neb2$li$Axis1*(-1)
y<-pca.neb2$li$Axis2<-pca.neb2$li$Axis2*(-1)

eigenv<-get_eigenvalue(pca.neb2)

circ<-rep("grey65",69)
circ[c(39,59,65)]<-"black"
col3<-col1
bgcol<-col1[pop(poolabies2)]

#Figure 3

tiff(filename=here("results","Figures","Figure3.tiff"), width=800,
     height = 680)

plot(x,y, bg=bgcol, pch=21, cex=1.7, col=circ, 
     xlab="PC1 (29.6%)", ylab="PC2 (14.1%)", cex.axis=1.3, cex.lab=1.4)

legend(2,-3.5,legend = c(expression(italic("A. alba")),
                         expression(italic("A. cephalonica")),
                         expression(paste(italic("A. nebrodensis")," adults")),
                         expression(paste(italic("A. nebrodensis"), " seedlings"))),
       pch=21, cex=1.2, col="grey35", pt.bg=col1, pt.cex=1.5,
       box.lwd = 0, border = "white", bg="transparent")
dev.off()

# legend(11,-15,legend = c("Silver fir", "Greek fir", "Sicilian fir adults", "Sicilian fir seedlings"), 
#        pch=21, cex=1.2, col="grey35", pt.bg=col1, pt.cex=1.5,
#        box.lwd = 0, border = "white", bg="transparent")


