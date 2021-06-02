
library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(here)


# En este script sigo los pasos de Fran para cargar y preparar los datos del vivero. La entrega de noviembre incluia las madres y casi 
# 400 plantulas del vivero. La segunda entrega (abril 2021) incluia el resto de muestras del vivero. En este segundo entregable hay un
# total de 7 SNPs que varian respecto al primer entregable. Asi, lo que hay que hacer es cargar las muestras del primer lote y crear
# un objeto genind y las muestras del segundo lote y crear por separado otro objeto genind. Luego se fusionaran ambos objetos genid

setwd("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/nursery")  # [1] "/Users/delValle/Dropbox/R/r nebrodi/jun21"
library(here) 
list.files()



# Placas 7-23 (entregable de abril en 2 archivos excel separados):
data<-here("abetos_p7p8p9p10p11p12P13P14_20210413_Results.csv")
x<-read.csv(data, header=T, skip = 17)
head(x)
data2<-here("ABETOS_P15P16P17P19P20P21P22P23_20210413_Results.csv")
y<-read.csv(data2, header=T, skip = 17)
head(y)

xy <- rbind(x,y)
dim(xy)  # [1] 184320      4
head(xy)
xy <-unite(xy, Genotype, Allele1, Allele2, sep="")

xy2<-xy %>%   group_by(SampleID, SNP_ID) 

xy3<-xy2 %>% dplyr::na_if("NOAMPNOAMP") 

xy4<-xy3 %>% dplyr::na_if("PRAPRA") 
GENOTYPE_CALLS_P7_23 <-xy4 %>% dplyr::na_if("UNDUND") 

dim(GENOTYPE_CALLS_P7_23)   # [1] 184320     3   
# Este numero cuadra con lo que sale al dividir el total de filas (184320) entre los 120 SNPs = 1536 muestras 

GENOTYPE_CALLS_P7_23 <- as.data.frame(GENOTYPE_CALLS_P7_23)
GENOTYPE_CALLS_P7_23<-as.data.frame(GENOTYPE_CALLS_P7_23[,c(2,1,3)]) 
# SNP_ID, Sample.ID y Genotype
GENOTYPE_CALLS_P7_23<-GENOTYPE_CALLS_P7_23[order(GENOTYPE_CALLS_P7_23$SNP_ID),] 
# Lo que obtenemos es una lista con las 1536 muestras ordenadas por orden alfabetico para cada uno de los 120 SNPs
head(GENOTYPE_CALLS_P7_23)



aa <- 1536 #change this number by the number of samples in your dataset  
cc <- 120  #change this number by the number of SNPs in your dataset
new_GENOTYPE_CALLS_P7_23 <- matrix(NA, aa, cc) #create empty matrix (en este caso de 1536x120)
a <- 1536  #change this number by the number of samples in your dataset
b <- 1 #initialization
for(j in 1:cc){
  b <- b
  a <- a  
  print(b)
  print(a)
  new_GENOTYPE_CALLS_P7_23[1:aa,j] <- GENOTYPE_CALLS_P7_23[b:a,3]
  b <- a+1
  a <- a+aa
}


colnames(new_GENOTYPE_CALLS_P7_23) <- unique(GENOTYPE_CALLS_P7_23[,1]) #assing SNP ID
rownames(new_GENOTYPE_CALLS_P7_23) <- unique(GENOTYPE_CALLS_P7_23[1:aa,2]) #assing Samples ID
new_GENOTYPE_CALLS_P7_23 <- as.data.frame(new_GENOTYPE_CALLS_P7_23) #output file
dim(new_GENOTYPE_CALLS_P7_23)  # 1536 120
head(new_GENOTYPE_CALLS_P7_23)

#              a00001151_105453 a00001795_101025 a00002911_27399 a00002963_19883 a00003814_72402 a00005964_11488 a00007680_38992 a00008388_38342
# 11_2008_0385               AG               AA              TC              CT              TT              AA              GC              GG
# 11_2008_0387               AA               AA              TC              CC              TT              AA              GC              GG
# 11_2008_0388               AG               AA              TC              CC              TG              AA              GC              AA
# 11_2008_0389               AG               AA              CC              CT              TT              AA              GC              GA
# 11_2008_0390               AA               AA              CC              CC              TT              AA              CC              GG
# 11_2008_0391             <NA>             <NA>            <NA>            <NA>            <NA>            <NA>            <NA>            <NA>






# Plascas 0 y 3-6 (entregable de finales de 2020):
# Nota: aqui se utilizaron algunos SNPs que posteriormente se descartaron porque fallaban en todas la muestras. Esto nos obliga
# a crear un objeto genind independiente para ambos lotes de muestras.

data3<-here("abetos_p0-p6_2020100_Results_vivaio.csv")
z<-read.csv(data3, header=T, sep=";")
head(z)
z<-unite(z, Genotype, Allele1, Allele2, sep="")

z2<-z %>%   group_by(SampleID, SNP_ID) 

z3<-z2 %>% dplyr::na_if("NOAMPNOAMP") 

z4<-z3 %>% dplyr::na_if("PRAPRA") 
GENOTYPE_CALLS_P3_6 <-z4 %>% dplyr::na_if("UNDUND") 

dim(GENOTYPE_CALLS_P3_6)   # [1] 48958     4  
# Aqu?? tengo replica en la placa 0, que es la que usamos a modo probatina la primera vez, que contenia 12 muestras
unique(GENOTYPE_CALLS_P3_6[c("SampleID")])  # 396 valores unicos en la columna SampleID (= 96 muestras x 8 placas)
# Este numero NO cuadra con lo que sale al dividir el total de filas entre los 120 SNPs: 48958/120 SNPs = 407.98 
# Deberia salir 408, que con 396 valores unicos + 12 del duplicado = 408
# Esto ocurre porque faltaban un par de SNPs de la muestra 08_2013_0002

new<-data.frame(SampleID=c("08_2013_0002","08_2013_0002"), SNP_ID=c("a00087283_35345","a00161543_9085"),Genotype=c(NA,NA)  )
GENOTYPE_CALLS_P3_6 <- as.data.frame(GENOTYPE_CALLS_P3_6)
GENOTYPE_CALLS_P3_6<-rbind(GENOTYPE_CALLS_P3_6,new)
dim(GENOTYPE_CALLS_P3_6)  # [1] 48960     3
# 48960 / 120 = 408 (= 12 repetidos + 396 valores unicos)



# Ahora es cuando aparece el problema de trabajar con duplicados. Lo que tengo que hacer es aplicar la solucion de Fran:
labels<-unite(GENOTYPE_CALLS_P3_6, label, SNP_ID, SampleID, sep="-" )$label
# Crea el objeto "labels" con una combinacion del nombre del SNP y sample, unidos por un guion
rs<-make.names(labels,unique=T)
# SNPs y samples estan separados por un punto. 
# Crea un elemento con etiquetas no repetidos que corresponden a nombres de fila unicos; de esta forma, se anade 
# .1, .2, etc., a la combinacion repetida de SNP+sample cuando son replicas (ej.  a00273547_1562.01_2019_2123.1)
Sample.ID2<-str_split(rs, "[.]", simplify = T, n = 2)[,2]
# Separa el vector rs creado anteriormente en una matriz de la que se conserva unicamente lo que queda
# a la derecha del punto (e.g. a00001151_105453.01_2019_2121 --> 01_2019_2121). Asi, obtenemos una matriz de 
# longitud 86400 filas con los nombres de cada sample para cada una de las filas


GENOTYPE_CALLS_P3_6$Sample.ID2<-Sample.ID2  
GENOTYPE_CALLS_P3_6<-as.data.frame(GENOTYPE_CALLS_P3_6[,c(2,4,3)]) 
# SNP_ID, Sample.ID2 y Genotype
GENOTYPE_CALLS_P3_6<-GENOTYPE_CALLS_P3_6[order(GENOTYPE_CALLS_P3_6$SNP_ID),] 
# Lo que obtenemos es una lista con las 408 muestras ordenadas por orden alfabetico para cada uno de los 120 SNPs
head(GENOTYPE_CALLS_P3_6)



aa <- 408 #change this number by the number of samples in your dataset  
cc <- 120  #change this number by the number of SNPs in your dataset
new_GENOTYPE_CALLS_P3_6 <- matrix(NA, aa, cc) #create empty matrix (en este caso de 768x120)
a <- 408  #change this number by the number of samples in your dataset
b <- 1 #initialization
for(j in 1:cc){
  b <- b
  a <- a  
  print(b)
  print(a)
  new_GENOTYPE_CALLS_P3_6[1:aa,j] <- GENOTYPE_CALLS_P3_6[b:a,3]
  b <- a+1
  a <- a+aa
}


colnames(new_GENOTYPE_CALLS_P3_6) <- unique(GENOTYPE_CALLS_P3_6[,1]) #assing SNP ID
rownames(new_GENOTYPE_CALLS_P3_6) <- unique(GENOTYPE_CALLS_P3_6[1:aa,2]) #assing Samples ID
new_GENOTYPE_CALLS_P3_6 <- as.data.frame(new_GENOTYPE_CALLS_P3_6) #output file
dim(new_GENOTYPE_CALLS_P3_6)  # 408 120


head(new_GENOTYPE_CALLS_P3_6)
#              a00001151_105453 a00001795_101025 a00002911_27399 a00002963_19883 a00003814_72402 a00005964_11488 a00007680_38992 a00008388_38342
# 13_2008_1210               AG               AA              TC              CT              TT              AA              GC              GG
# 13_2008_1211               GG               AA              TC              CT              TG              AA              GG            <NA>
# 13_2008_1212               AG               AA              CC              CC              TT              AA              GC              GG
# 13_2008_1213               AA               AA              CC              CC              TG              AA            <NA>            <NA>
# 13_2008_1214               AA               AA              TT              CC              TG              AA              CC              GG
# 13_2008_1215               AA               AA              CC              CC              TT              AA              CC              GA   



# Ya tengo las matrices preparadas. Voy a crear los 2 objetos genind:

library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)

nebrodf_P7_23 <-df2genind(new_GENOTYPE_CALLS_P7_23, ncode = 1,
                   NA.char = "NA",
                   ploidy = 2,
                   type = "codom")
# Con la matriz preparada (new_GENOTYPE_CALLS_P7_23), que es un data frame, crea un objeto genind con los siguientes parametros:
# 1) ncode = 1; numero de caracteres necesarios para codificar un genotipo en un locus
# 2) NA.char = NA; se les ha llamado NA, son los missing data 
# 3) ploidy = 2; el nivel de ploid??a
# 4) type = "codom"; son alelos codominantes


# Warning messages:
#   1: In df2genind(new_GENOTYPE_CALLS_P7_23, ncode = 1, NA.char = "NA",  :
#   Markers with no scored alleles have been removed
#   2: In df2genind(new_GENOTYPE_CALLS_P7_23, ncode = 1, NA.char = "NA",  :
#   Individuals with no scored loci have been removed

nebrodf_P7_23
# /// GENIND OBJECT /////////
#   
#   // 1,381 individuals; 119 loci; 230 alleles; size: 1.4 Mb   --> Teniamos 1536 individuos (se han perdido 155) y 120 loci (se ha perdido 1 locus)
summary(nebrodf_P7_23)
# Percentage of missing data: 8.91 %

failedind<-new_GENOTYPE_CALLS_P7_23[(rownames(new_GENOTYPE_CALLS_P7_23)%in%rownames(nebrodf_P7_23@tab))==FALSE,]
## 155 ind

# loci que han fallado completamente
failedloc<-colnames(new_GENOTYPE_CALLS_P7_23)[colnames(new_GENOTYPE_CALLS_P7_23)%in%names(nebrodf_P7_23@all.names)==FALSE]

#                a3
# 1 a00181436_13305
rm(GENOTYPE_CALLS_P7_23_delete_me)




nebrodf_P3_6 <-df2genind(new_GENOTYPE_CALLS_P3_6, ncode = 1,
                          NA.char = "NA",
                          ploidy = 2,
                          type = "codom")



nebrodf_P3_6
# /// GENIND OBJECT /////////
#   
#   // 408 individuals; 112 loci; 217 alleles; size: 443.8 Kb
summary(nebrodf_P3_6)
# Percentage of missing data: 11.95 %


locioverlapping<-locNames(nebrodf_P7_23)[locNames(nebrodf_P7_23)%in%locNames(nebrodf_P3_6)==TRUE]

nursery<-repool(nebrodf_P7_23[loc=locioverlapping], nebrodf_P3_6[loc=locioverlapping], list=FALSE) 

summary(nursery)


# CONSENSO DUPLICADOS

samples<-str_split(rownames(nursery@tab), "[.]", simplify = T, n = 2)[,1]

DF<-as.data.frame(nursery@tab)
DF$ID=samples



DF2<-DF %>%
  group_by(ID) %>%
  summarise_each(funs(median(., na.rm = T)))

DF3<-as.data.frame(DF2)
DF3[DF3==-Inf]<-NA
rownames(DF3)<-DF3$ID
DF4<-DF3[,-1]
DF4<-DF4[sort(rownames(DF4)),]

#dataset sin duplicados
nursery2<-as.genind(DF4)

summary(nursery2)



#FORMATEO PARA COLONY

library(dplyr)
library(outliers)
library(strataG)
library(pegas)

#Preparación del data set



nursery_loci<-genind2loci(nursery2)
summary(nursery_loci)



nursery_loci_colony<-as.data.frame(nursery_loci)
nursery_loci_colony2<-alleleSplit(nursery_loci_colony[,-1], sep = "/")
rownames(nursery_loci_colony2)<-rownames(nursery_loci_colony)


nursery_loci_colony2[is.na(nursery_loci_colony2)] <- 0

changeNucleotides2numbers<-function(juve){
  juve[juve=="A"]<-'100'
  juve[juve=="T"]<-'101'
  juve[juve=="C"]<-'102'
  juve[juve=="G"]<-'103'
  juve[is.na(juve)] <-'0'
  juve
}

nursery_loci_colony3<-changeNucleotides2numbers(nursery_loci_colony2)

write.table(nursery_loci_colony3,"./nursery_colony",
            sep="\t",quote=T)


# comprobar orden loci
# asignar grupos de plantas a madres
