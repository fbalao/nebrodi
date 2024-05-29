##### Import A. nebrodensis nursery seedlings genotypes in Openarray format to genind

#Libraries and functions ----
library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(here)
library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)
library(dplyr)
library(outliers)
library(strataG)
library(pegas)

# changeNucleotides2numbers is a function to change nucleotides to numbers
changeNucleotides2numbers <- function(juve){
  juve[juve=="A"] <- '100'
  juve[juve=="T"] <- '101'
  juve[juve=="C"] <- '102'
  juve[juve=="G"] <- '103'
  juve[is.na(juve)] <- '0'
  juve
}

# Load OpenArray genotypes, transform and rename replicates  ----
# Plates 7-23 (nursery seedlings genotypes)
data<-here("Abies_nebrodi","openarray_data","nursery","abetos_p7p8p9p10p11p12P13P14_20210413_Results.csv")
x<-read.csv(data, header=T, skip = 17)
head(x)
data2<-here("Abies_nebrodi","openarray_data","nursery","ABETOS_P15P16P17P19P20P21P22P23_20210413_Results.csv")
y<-read.csv(data2, header=T, skip = 17)
head(y)

xy <- rbind(x,y)
dim(xy)  # Dimensions [1] 184320 4
head(xy)
xy <-unite(xy, Genotype, Allele1, Allele2, sep="")

xy2<-xy %>%   group_by(SampleID, SNP_ID) 

xy2[xy2 =="NOAMPNOAMP"] <- NA
xy2[xy2 =="PRAPRA"] <- NA
xy2[xy2 =="UNDUND"] <- NA

GENOTYPE_CALLS_P7_23<-xy2 

dim(GENOTYPE_CALLS_P7_23)   # Dimensions [1] 184320 3
# This number matches what comes out when dividing the total rows (184320) by the 120 SNPs = 1536 samples

GENOTYPE_CALLS_P7_23 <- as.data.frame(GENOTYPE_CALLS_P7_23)
GENOTYPE_CALLS_P7_23<-as.data.frame(GENOTYPE_CALLS_P7_23[,c(2,1,3)]) 
# SNP_ID, Sample.ID, and Genotype
GENOTYPE_CALLS_P7_23<-GENOTYPE_CALLS_P7_23[order(GENOTYPE_CALLS_P7_23$SNP_ID),] 
# What we get is a list with the 1536 samples sorted alphabetically for each of the 120 SNPs
head(GENOTYPE_CALLS_P7_23)


# Plates 0 and 3-6 (deliverable from late 2020):
# Note: some SNPs were used here that were later discarded because they failed in all samples. This requires
# creating a separate genind object for both sets of samples.

data3<-here("Abies_nebrodi","openarray_data","nursery","abetos_p0-p6_2020100_Results_vivaio.csv")
z<-read.csv(data3, header=T, sep=";")
head(z)
z<-unite(z, Genotype, Allele1, Allele2, sep="")

z2<-z %>%   group_by(SampleID, SNP_ID) 

z2[z2 =="NOAMPNOAMP"] <- NA
z2[z2 =="PRAPRA"] <- NA
z2[z2 =="UNDUND"] <- NA 
GENOTYPE_CALLS_P3_6 <-z2

dim(GENOTYPE_CALLS_P3_6)   # Dimensions [1] 48958 4
# Here I have a replicate on plate 0, which we used as a trial the first time, containing 12 samples
unique(GENOTYPE_CALLS_P3_6[c("SampleID")])  # 396 unique values in the SampleID column (= 96 samples x 8 plates)
# This number does NOT match what comes out when dividing the total rows by the 120 SNPs: 48958/120 SNPs = 407.98
# Should be 408, which with 396 unique values + 12 duplicates = 408
# This occurs because a couple of SNPs from sample 08_2013_0002 were missing

new<-data.frame(SampleID=c("08_2013_0002","08_2013_0002"),
                SNP_ID=c("a00087283_35345","a00161543_9085"),Genotype=c(NA,NA)) # Fix the empty genotypes
GENOTYPE_CALLS_P3_6 <- as.data.frame(GENOTYPE_CALLS_P3_6)
GENOTYPE_CALLS_P3_6<-rbind(GENOTYPE_CALLS_P3_6,new)
dim(GENOTYPE_CALLS_P3_6)  # Dimensions [1] 48960 3
# 48960 / 120 = 408 (= 12 duplicates + 396 unique values)

# Now the issue of working with duplicates arises. What I need to do is apply Fran's solution:
labels<-unite(GENOTYPE_CALLS_P3_6, label, SNP_ID, SampleID, sep="-" )$label
# Creates the "labels" object with a combination of the SNP name and sample, joined by a hyphen
rs<-make.names(labels,unique=T)
# SNPs and samples are separated by a period.
# Creates a unique element with unique row names; this way, .1, .2, etc., are added to the repeated SNP+sample combination when they are duplicates (e.g. a00273547_1562.01_2019_2123.1)
Sample.ID2<-str_split(rs, "[.]", simplify = T, n = 2)[,2]
# Splits the previously created rs vector into a matrix, keeping only what is to the right of the period (e.g., a00001151_105453.01_2019_2121 --> 01_2019_2121). Thus, we obtain a matrix of
# 86400 rows with the names of each sample for each row

GENOTYPE_CALLS_P3_6$Sample.ID2<-Sample.ID2  
GENOTYPE_CALLS_P3_6<-as.data.frame(GENOTYPE_CALLS_P3_6[,c(2,4,3)]) 
# SNP_ID, Sample.ID2, and Genotype
GENOTYPE_CALLS_P3_6<-GENOTYPE_CALLS_P3_6[order(GENOTYPE_CALLS_P3_6$SNP_ID),] 
# What we get is a list with the 408 samples sorted alphabetically for each of the 120 SNPs
head(GENOTYPE_CALLS_P3_6)

# Transform the genind objects ----

nebrodf_P7_23 <- df2genind(new_GENOTYPE_CALLS_P7_23, ncode = 1,
                           NA.char = "NA",
                           ploidy = 2,
                           type = "codom")

nebrodf_P7_23
# /// GENIND OBJECT /////////
#   
#   // 1,381 individuals; 119 loci; 230 alleles; size: 1.4 Mb   --> Originally 1536 individuals (155 lost) and 120 loci (1 locus lost)
summary(nebrodf_P7_23)
# Percentage of missing data: 8.91 %

failedind <- new_GENOTYPE_CALLS_P7_23[(rownames(new_GENOTYPE_CALLS_P7_23)%in%rownames(nebrodf_P7_23@tab))==FALSE,]
## 155 individuals

# loci that have completely failed
failedloc <- colnames(new_GENOTYPE_CALLS_P7_23)[colnames(new_GENOTYPE_CALLS_P7_23)%in%names(nebrodf_P7_23@all.names)==FALSE]

#                a3
# 1 a00181436_13305
rm(GENOTYPE_CALLS_P7_23_delete_me)

nebrodf_P3_6 <- df2genind(new_GENOTYPE_CALLS_P3_6, ncode = 1,
                          NA.char = "NA",
                          ploidy = 2,
                          type = "codom")

nebrodf_P3_6
# /// GENIND OBJECT /////////
#   
#   // 408 individuals; 112 loci; 217 alleles; size: 443.8 Kb
summary(nebrodf_P3_6)
# Percentage of missing data: 11.95 %


locioverlapping <- locNames(nebrodf_P7_23)[locNames(nebrodf_P7_23)%in%locNames(nebrodf_P3_6)==TRUE]

nursery <- repool(nebrodf_P7_23[loc=locioverlapping], nebrodf_P3_6[loc=locioverlapping], list=FALSE) 

summary(nursery)


# Consensus for duplicates ----

samples <- str_split(rownames(nursery@tab), "[.]", simplify = T, n = 2)[,1]

DF <- as.data.frame(nursery@tab)
DF$ID = samples

DF2 <- DF %>%
  group_by(ID) %>%
  summarise_each(funs(median(., na.rm = T)))

DF3 <- as.data.frame(DF2)
DF3[DF3==-Inf] <- NA
rownames(DF3) <- DF3$ID
DF4 <- DF3[,-1]
DF4 <- DF4[sort(rownames(DF4)),]

# Dataset without duplicates
nursery2 <- as.genind(DF4)

summary(nursery2)

save.image("nursery.RData")

# Formatting genotypes to Colony2  ----
# Preparing the dataset

nursery_loci <- genind2loci(nursery2)
summary(nursery_loci)

nursery_loci_colony <- as.data.frame(nursery_loci)
nursery_loci_colony2 <- alleleSplit(nursery_loci_colony[,-1], sep = "/")
rownames(nursery_loci_colony2) <- rownames(nursery_loci_colony)

nursery_loci_colony2[is.na(nursery_loci_colony2)] <- 0



nursery_loci_colony3 <- changeNucleotides2numbers(nursery_loci_colony2)

write.table(nursery_loci_colony3, file= here("data","formats","nursery_colony"),
            sep="\t", quote=T)

# Check the order of loci
# Assign plant groups to mothers
