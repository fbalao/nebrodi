# Tareas
# 1. Limpiar madres y seleccionar loci que funcionen
# 2. Limpiar resto de datos
# 3. Comprobar los loci
#   3.1 Homocigotos -> set hibridación?
#   3.2 No funcionan -> elegir nuevos?
# 4. Comprobar individuos que han fallado
# 5. Comprobar el poder de asignación SOLOMON
#   

# libraries ----
library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(here)
library(adegenet)
library(poppr)
library(plotrix)
library(dartR)
# Carga, trasnformación y renombrado de réplicas ----

data<-here("Abies_nebrodi", "openarray_data","abetos_p0-p6_2020100_Results.csv")
data2<-here("Abies_nebrodi", "openarray_data","abetos_p1p2triplicados_20201105_Results.csv")
x<-read.csv(data, header=T, skip = 17)
xr<-read.csv(data2, header=T, skip = 17)
colnames(xr)<-colnames(x)

x<-rbind(x,xr)
  
x<-unite(x, Genotype,Allele1,Allele2, sep="")

x2<-x %>%   group_by(SampleID, SNP_ID) 
x2[x2 =="NOAMPNOAMP"] <- NA
x2[x2 =="PRAPRA"] <- NA
x2[x2 =="UNDUND"] <- NA
GENOTYPE_CALLS<-x2


#Arreglo provisional
dim(GENOTYPE_CALLS)
new<-data.frame(Sample.ID2=c("08_2013_0002_HXL05","08_2013_0002_HXL05"),SNP_ID=c("a00087283_35345","a00161543_9085"),Genotype=c(NA,NA)  )
GENOTYPE_CALLS<-rbind(GENOTYPE_CALLS,new)
dim(GENOTYPE_CALLS)


labels<-unite(GENOTYPE_CALLS, label, SNP_ID, SampleID,, sep="-" )$label
rs<-make.names(labels,unique=T)
Sample.ID2<-str_split(rs, "[.]", simplify = T, n = 2)[,2]

GENOTYPE_CALLS$Sample.ID2<-Sample.ID2
GENOTYPE_CALLS2<-as.data.frame(GENOTYPE_CALLS[,c(3,5,4)])
GENOTYPE_CALLS2<-GENOTYPE_CALLS2[order(GENOTYPE_CALLS2$SNP_ID),]

aa <- 1248 #change this number by the number of samples in your dataset
cc <- 120  #change this number by the number of SNPs in your dataset
new_GENOTYPE_CALLS <- matrix(NA, aa, cc) #create empty matrix
a <- 1248  #change this number by the number of samples in your dataset
b <- 1 #initialization
for(j in 1:cc){
  b <- b
  a <- a  
  print(b)
  print(a)
  new_GENOTYPE_CALLS[1:aa,j] <- GENOTYPE_CALLS2[b:a,3]
  b <- a+1
  a <- a+aa
}
colnames(new_GENOTYPE_CALLS) <- unique(GENOTYPE_CALLS2[,1]) #assing SNP ID
rownames(new_GENOTYPE_CALLS) <- unique( GENOTYPE_CALLS2[1:aa,2]) #assing Samples ID
new_GENOTYPE_CALLS <- as.data.frame(new_GENOTYPE_CALLS)#output file


# Importación a adegenet ----
library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)

nebrodf<-df2genind(new_GENOTYPE_CALLS,  ncode = 1,
                   NA.char = "NA",
                   ploidy = 2,
                   type = "codom")

# loci que han fallado completamente
failedloc<-colnames(new_GENOTYPE_CALLS[colnames(new_GENOTYPE_CALLS)%in%locNames(nebrodf)==FALSE])


#individuos que han fallado completamente
failedind<-rownames(new_GENOTYPE_CALLS[(rownames(new_GENOTYPE_CALLS)%in%rownames(nebrodf@tab))==FALSE,])

# loci homocigotos
names(nebrodf$loc.n.all)[nebrodf$loc.n.all==1]

#Asignar las madres ----

samples<-str_split(rownames(nebrodf@tab), "[.]", simplify = T, n = 2)[,1]

# CONSENSO

DF<-as.data.frame(nebrodf@tab)
DF$ID=samples

DF[DF$ID=="01_2019_2257", "a00161543_9085.G" ]<-0 #Detectado una replica no concordante con genotipo 1/1 en vez de 2/0. Cambiado
# HAY MÁS DISCORDANCIAS, REPASAR

# DF2<-DF %>%
#   group_by(ID) %>%
#   summarise_each(funs(max(., na.rm = T)))
# 

DF2<-DF %>%
  group_by(ID) %>%
  summarise_each(list(~median(., na.rm = T)))

DF3<-as.data.frame(DF2)
DF3[DF3==-Inf]<-NA
rownames(DF3)<-DF3$ID
DF4<-DF3[,-1]


tx<-function(x){sum(is.na(x))}
mean(apply(DF4,1,tx)/219)
mean(apply(nebrodf@tab,1,tx)/219)
# Cambio nombre

DF4<-DF4[sort(rownames(DF4)),]
nebrodf2<-as.genind(DF4)


#Filtrado  # se puede filtrar tras arreglar las madres

#nebrodf<-missingno(nebrodf, type = "geno", cutoff = 0.50, quiet = FALSE, freq = FALSE)
# se eliminan 10 muestras más
ty<-function(x){sum(is.na(x))}

y<-matrix(nrow = 21, ncol=2)
for (z in 1:21){
  i<-seq(0,1,0.05)[z]
  db<-missingno(nebrodf2, type = "loci", cutoff = i, quiet = F, freq = FALSE)
  y[z,1]<-nlevels(db@loc.fac)
  y[z,2]<-mean(apply(db@tab,1,ty)/(y[z,1]*2))
}

cuttoff<-cbind(seq(0,1,0.05),y)
cuttNA<-cbind(seq(0,1,0.05),y[,2])
plot( cuttoff,type="b", ylab="N loci", xlab="cut off call rate", ylim=c(0,120))
plot( cuttNA,type="b", ylab="NA %", xlab="cut off call rate")
abline(v=0.25, lty=2)


y<-matrix(nrow = 21, ncol=2)
for (z in 1:21){
  i<-seq(0,1,0.05)[z]
  db<-missingno(nebrodf2, type = "geno", cutoff = i, quiet = F, freq = FALSE)
  y[z,1]<-dim(db@tab)[1]
  y[z,2]<-mean(apply(db@tab,1,ty)/(dim(db@tab)[2]))
}

cuttoff<-cbind(seq(0,1,0.05),y[,1])
cuttNA<-cbind(seq(0,1,0.05),y[,2])
plot( cuttoff,type="b", ylab="N individuals", xlab="cut off call rate")
plot(cuttNA, type="b", col=2)
abline(v=0.25, lty=2)

nebrodf2<-missingno(nebrodf2, type = "geno", cutoff = 0.30, quiet = FALSE, freq = FALSE)

nebrodf2<-missingno(nebrodf2, type = "loci", cutoff = 0.50, quiet = FALSE, freq = FALSE)

codes<-here("Abies_nebrodi", "openarray_data","muestras.csv")
codes<-read.table(codes, header=T)
codes$Codigo_muestra%in%rownames(nebrodf2@tab)

db<-codes[match(rownames(nebrodf2@tab), codes$Codigo_muestra),]

nebrodf2@other$Plantula_madre<-as.factor(db$Plantula_madre)
nebrodf2@other$Vivero_campo<-as.factor(db$Vivero_campo)
nebrodf2@other$Madre_procedencia<-as.factor(db$Madre_procedencia)



nebrodf2@pop<-as.factor(rep("Anebrodi",600))


# Estima del fallo de genotipado ----

jd<-gi2gl(nebrodf2)
nakd<-NA.posi(jd)

hist(table(unlist(nakd))/598)

mean(table(unlist(nakd))/598)
std.error(as.vector(table(unlist(nakd))/598))

# 0.05760301 +- 0.009546861 basado en 598 muestras
# la mediana es 0.01170569

error<-DF %>%
  group_by(ID) %>%
  summarise_each(list(~sum(is.na(.))))


col_odd <- seq_len(ncol(DF)) %% 2
DF_odd <- DF[ , col_odd == 1] 
error<-DF_odd %>%
  summarise_each(list(~sum(is.na(.)/1227)))
