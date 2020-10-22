# Tareas
# 1. Limpiar madres y seleccionar loci que funcionen
# 2. Limpiar resto de datos
# 3. Comprobar los loci
#   3.1 Homocigotos -> set hibridación?
#   3.2 No funcionan -> elegir nuevos?
# 4. Comprobar individuos que han fallado
# 5. Comprobar el poder de asignación SOLOMON
#   


library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(here)

# Carga, trasnformación y renombrado de réplicas
data<-here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","abetos_p0-p6_2020100_Results.csv")
x<-read.csv(data, header=T, skip = 17)
x<-unite(x, Genotype,Allele1,Allele2, sep="")
x<-unite(x, Sample.ID2,SampleID,PlateBarcode, sep="_")
x2<-x %>%   group_by(Sample.ID2) 

x3<-x2 %>% dplyr::na_if("NOAMPNOAMP") 

x4<-x3 %>% dplyr::na_if("PRAPRA") 
GENOTYPE_CALLS<-x4 %>% dplyr::na_if("UNDUND") 


#Arreglo provisional
dim(GENOTYPE_CALLS)
new<-data.frame(Sample.ID2=c("08_2013_0002_HXL05","08_2013_0002_HXL05"),SNP_ID=c("a00087283_35345","a00161543_9085"),Genotype=c(NA,NA)  )
GENOTYPE_CALLS<-rbind(GENOTYPE_CALLS,new)
dim(GENOTYPE_CALLS)

labels<-unite(GENOTYPE_CALLS, label, SNP_ID, Sample.ID2, sep="-" )$label
rs<-make.names(labels,unique=T)
Sample.ID2<-str_split(rs, "[.]", simplify = T, n = 2)[,2]

GENOTYPE_CALLS$Sample.ID2<-Sample.ID2
GENOTYPE_CALLS2<-as.data.frame(GENOTYPE_CALLS[,c(2,1,3)])
GENOTYPE_CALLS2<-GENOTYPE_CALLS2[order(GENOTYPE_CALLS2$SNP_ID),]

aa <- 720 #change this number by the number of samples in your dataset
cc <- 120  #change this number by the number of SNPs in your dataset
new_GENOTYPE_CALLS <- matrix(NA, aa, cc) #create empty matrix
a <- 720  #change this number by the number of samples in your dataset
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


# Importación a adegenet
library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)

nebrodf<-df2genind(new_GENOTYPE_CALLS,  ncode = 1,
                   NA.char = "NA",
                   ploidy = 2,
                   type = "codom")

# loci que han fallado completamente
failedloc<-colnames(new_GENOTYPE_CALLS)[colnames(new_GENOTYPE_CALLS)%in%names(nebrodf@all.names)==FALSE]


#individuos que han fallado completamente
failedind<-new_GENOTYPE_CALLS[(rownames(new_GENOTYPE_CALLS)%in%rownames(nebrodf@tab))==FALSE,]

# loci homocigotos
names(nebrodf$loc.n.all)[nebrodf$loc.n.all==1]

#Asignar las madres

mothers<-read.table(here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","parental_list.txt"),header=T)

#nombres cudadriplicados
parentals_names<-c(mothers$parentals,paste0(rep(mothers$parentals,4),".", 1:4))

type<-factor(rep("sibling", 708), levels=c("sibling","parental"))
table(rownames(nebrodf@tab)%in%parentals_names) # hay 48 padres
type[rownames(nebrodf@tab)%in%parentals_names]="parental"

nebrodf@pop<-type


#Filtrado

#nebrodf<-missingno(nebrodf, type = "geno", cutoff = 0.50, quiet = FALSE, freq = FALSE)
# se eliminan 10 muestras más

nebrodf<-missingno(nebrodf, type = "loci", cutoff = 0.25, quiet = FALSE, freq = FALSE)

 y<-numeric(21)
 for (z in 1:21){
   i<-seq(0,1,0.05)[z]
   y[z]<-nlevels((missingno(nebrodf, type = "loci", cutoff = i, quiet = F, freq = FALSE))@loc.fac)
 }
 
 cuttoff<-cbind(seq(0,1,0.05),y)
 plot( cuttoff,type="b", ylab="N loci", xlab="cut off call rate")
 abline(v=0.25, lty=2)
 


#separar los padres e hijos
nebrodf2<-seppop(nebrodf)



tx<-function(x){sum(is.na(x))}
apply(nebrodf@tab,1,tx)

mothernames<-substr(rownames(nebrodf@tab),1,12)
nebrodf@other$mother<-mothernames




DF<-as.data.frame(nebrodf@tab)
DF$ID=nebrodf@other$mother

DF2<-DF %>%
  group_by(ID) %>%
  summarise_each(funs(max(., na.rm = T)))

DF3<-as.data.frame(DF2)
DF3[DF3==-Inf]<-NA
rownames(DF3)<-DF3$ID
DF4<-DF3[,-1]
apply(DF4,1,tx)
# Cambio nombre

DF4<-DF4[sort(rownames(DF4)),]

m<-read.table(here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","codigosmadre.csv"), sep=" ", header=T)

rownames(DF4)<-m$mother_code
nebrodf<-as.genind(DF4)
nebrodf@pop<-factor(rep("P",30))

x<-missingno(nebrodf, type = "loci", cutoff = 0.25, quiet = FALSE, freq = FALSE)
summary(x)

### missing
missing_data <- info_table(nebrodf, type = "missing")
sum(missing_data["Total", 1:104] > 0)
barplot(missing_data["Total", 1:104], xlab = "Locus", ylab = "% Fallo")
abline(h=0.2, lty=2)

# cambio codigos madre


#diversidad

div<-summary(nebrodf)

basicstat <- basic.stats(nebrodf, diploid = TRUE, digits = 2) 

#Prueba del Equilibrio de Hardy-Weinberg

hw<-hw.test(nebrodf, B =1000)
out<-hw[,3]<0.0001
plot(div$Hobs, div$Hexp, pch=16, col=as.numeric(out)+1)
abline(0,1)


boot.ppfis(nebrodf)


x <-indpca(nebrodf, ind.labels = row.names(nebrodf$tab)) 
plot(x, cex =0.4)

nebrotab<-tab(nebrodf, freq=T, NA.method="mean")
pca.neb <- dudi.pca(nebrotab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.neb$li, fac=as.factor(row.names(nebrotab)),
        col=rainbow(30), clabel = 0.7)

distgen<- dist(nebrodf, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)    
plot(bionjs(distgen), type="u", cex=0.7)

nebrodf2 <-genind2loci(nebrodf)


distgenDIFF <- dist.gene(nebrodf2, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
plot(bionj(distgenDIFF), type="p", cex=0.7)


distgenDISS <- diss.dist(nebrodf, percent = FALSE, mat = FALSE)
plot(bionj(distgenDISS), type="p", cex=0.7)









