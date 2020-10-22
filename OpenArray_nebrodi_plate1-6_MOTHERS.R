library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(here)


# 1. Cargar los datos
# 2. Transformar datos a adegenet
# 3. Limpiar los datos
#   3.1 Eliminar marcadores fallidos
#   3.2 Eliminar individuos fallidos
# 4. Padres
#   4.1 Consenso de réplicas
#   4.2 Diversidad


data<-here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","abetos_p0-p6_2020100_Results.csv")
x<-read.csv(data, header=T, skip = 17)
x<-unite(x, Genotype,Allele1,Allele2, sep="")

x2<-x %>%   group_by(SampleID, SNP_ID) 

x3<-x2 %>% dplyr::na_if("NOAMPNOAMP") 

x4<-x3 %>% dplyr::na_if("PRAPRA") 
GENOTYPE_CALLS<-x4 %>% dplyr::na_if("UNDUND") 


#Arreglo provisional
dim(GENOTYPE_CALLS)
new<-data.frame(SampleID=c("08_2013_0002","08_2013_0002"), PlateBarcode=c("HXL05","HXL05"),SNP_ID=c("a00087283_35345","a00161543_9085"),Genotype=c(NA,NA)  )
GENOTYPE_CALLS<-rbind(GENOTYPE_CALLS,new)
dim(GENOTYPE_CALLS)

labels<-unite(GENOTYPE_CALLS, label, SNP_ID, SampleID,, sep="-" )$label
rs<-make.names(labels,unique=T)
Sample.ID2<-str_split(rs, "[.]", simplify = T, n = 2)[,2]

GENOTYPE_CALLS$Sample.ID2<-Sample.ID2
GENOTYPE_CALLS2<-as.data.frame(GENOTYPE_CALLS[,c(3,5,4)])
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

# library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)

nebrodf<-df2genind(new_GENOTYPE_CALLS,  ncode = 1,
  NA.char = "NA",
  ploidy = 2,
  type = "codom")

#Asignar las madres

mothers<-read.table(here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","parental_list.txt"),header=T)

#nombres cudadriplicados
parentals_names<-c(mothers$parentals,paste0(rep(mothers$parentals,4),".", 1:4))

type<-factor(rep("sibling", 708), levels=c("sibling","parental"))
 table(rownames(nebrodf@tab)%in%parentals_names) # hay 48 padres
 type[rownames(nebrodf@tab)%in%parentals_names]="parental"
 
nebrodf@pop<-type


# #Filtrado
# nebrodf<-missingno(nebrodf, type = "geno", cutoff = 0.50, quiet = FALSE, freq = FALSE)
# 
# nebrodf<-missingno(nebrodf, type = "loci", cutoff = 0.45, quiet = FALSE, freq = FALSE)
# 


#separar los padres e hijos
nebrodf2<-seppop(nebrodf)

parentaldf<-nebrodf2$parental
#parentaldf<-missingno(parentaldf, type = "loci", cutoff = 0.45, quiet = FALSE, freq = FALSE)
parentaldf

tx<-function(x){sum(is.na(x))}
apply(parentaldf@tab,1,tx)

mothernames<-substr(rownames(parentaldf@tab),1,12)
parentaldf@other$mother<-mothernames




DF<-as.data.frame(parentaldf@tab)
DF$ID=parentaldf@other$mother

DF2<-DF %>%
  group_by(ID) %>%
  summarise_each(funs(max(., na.rm = T)))

DF3<-as.data.frame(DF2)
DF3[DF3==-Inf]<-NA ###
rownames(DF3)<-DF3$ID
DF4<-DF3[,-1]
apply(DF4,1,tx)
# Cambio nombre
           
DF4<-DF4[sort(rownames(DF4)),]

m<-read.table(here("Datos","ARTICULO","Abies_nebrodi", "openarray_data","codigosmadre.csv"), sep=" ", header=T)
            
rownames(DF4)<-m$mother_code


cleaned_mothers<-as.genind(DF4)
cleaned_mothers@pop<-factor(rep("P",30))

cleaned_mothers2<-informloci(cleaned_mothers, cutoff = 0, MAF = 0.01, quiet = FALSE)

loci<-genind2loci(cleaned_mothers2)
summary(loci)# check a00255117_6348
write.loci(loci, loci.sep = "\t", quote = FALSE, col.names = T, file="moGeno_nebro.txt")



x<-missingno(cleaned_mothers, type = "loci", cutoff = 0.25, quiet = FALSE, freq = FALSE)
summary(x)

### missing
missing_data <- info_table(cleaned_mothers, type = "missing")
sum(missing_data["Total", 1:104] > 0)
barplot(missing_data["Total", 1:104], xlab = "Locus", ylab = "% Fallo")
abline(h=0.2, lty=2)

# cambio codigos madre


#diversidad

div<-summary(cleaned_mothers)

basicstat <- basic.stats(cleaned_mothers, diploid = TRUE, digits = 2) 

#Prueba del Equilibrio de Hardy-Weinberg

hw<-hw.test(cleaned_mothers, B =1000)
out<-hw[,3]<0.01
plot(div$Hobs, div$Hexp, pch=16, col=as.numeric(out)+1)
abline(0,1)


boot.ppfis(cleaned_mothers)


x <-indpca(cleaned_mothers, ind.labels = row.names(cleaned_mothers$tab)) 
plot(x, cex =0.7, col=rainbow(30))

 nebrotab<-tab(cleaned_mothers, freq=T, NA.method="mean")
 pca.neb <- dudi.pca(nebrotab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
 s.class(pca.neb$li, fac=as.factor(row.names(nebrotab)),
         col=rainbow(30), clabel = 0.7)
 
distgen<- dist(cleaned_mothers, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)    
plot(bionjs(distgen), type="p", cex=0.7)

cleaned_mothers2 <-genind2loci(cleaned_mothers)


distgenDIFF <- dist.gene(cleaned_mothers2, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
plot(bionj(distgenDIFF), type="p", cex=0.7)


distgenDISS <- diss.dist(cleaned_mothers, percent = FALSE, mat = FALSE)
plot(bionj(distgenDISS), type="p", cex=0.7)









