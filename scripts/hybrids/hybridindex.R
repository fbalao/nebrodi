# Libraries ----
library(triangulaR)
library(gghybrid)
library(coda)
library(here)
library(vcfR)
library(radiator)
source(here("scripts","hybrids","plot_h2.R"))

load(here("envs","nebrodensis.21.12.2020.RData"))

# Pool Abies species and A. nebrodensis hybrids
vcf_abies <- read.vcfR("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/hybrids/R0.3populations_TODOsnps.recode.vcf", verbose = FALSE )
markers<-read.table("Abies_nebrodi/openarray_data/hybrids/markers.recode.txt", header=T)

abies_genind <- vcfR2genind(vcf_abies, return.alleles=T)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
abies_genind@pop<-as.factor(pop)
abiesnames<-read.table("namesabiespecies.txt", sep=",")
indNames(abies_genind)<-abiesnames[match(indNames(abies_genind), abiesnames$V1), 2]
locNames(abies_genind)<-markers[match(locNames(abies_genind), markers$ID), "loci"]

abies_genindc<-abies_genind[6:35]


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

hybrid_alba<-poolabies[16:69]

hybrid_ceph<-poolabies[c(1:15,31:69)]

# TriangulaR ----

# Export to vcf - a intermediate step to triangulaR

genomic_converter(
  data =hybrid_alba, output = "vcf", filename = "hybrid_alba")

genomic_converter(
  data =hybrid_ceph, output = "vcf", filename = "hybrid_ceph")


# Import VCF

vcf_hybrid_alba <- read.vcfR(here("data", "hybrids","hybrid_alba.vcf"), verbose = FALSE )
vcf_hybrid_alba@fix[,3]<-as.character(1:95)

vcf_hybrid_ceph <- read.vcfR(here("data", "hybrids","hybrid_ceph.vcf"), verbose = FALSE )
vcf_hybrid_ceph@fix[,3]<-as.character(1:95)

#popmap

popabies<-c(rep("A. nebrodensis", 39 ),rep("A. alba", 15))
popabies[grep(".P", colnames(vcf_hybrid_alba@gt)[-1])]<-"A. nebrodensis seedlings"
popabies[33]<-"A. nebrodensis hybrid"
vcf_hybrid_alba_popmap<-data.frame(id=colnames(vcf_hybrid_alba@gt)[-1], pop=popabies)


popceph<-c(rep("A. nebrodensis", 39 ),rep("A. cephalonica", 15))
popceph[grep(".P", colnames(vcf_hybrid_ceph@gt)[-1])]<-"seedlings"
popceph[33]<-"A. nebrodensis hybrid"

vcf_hybrid_ceph_popmap<-data.frame(id=colnames(vcf_hybrid_ceph@gt)[-1], pop=popceph)



# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
alba.vcfR.diff <- alleleFreqDiff(vcfR = vcf_hybrid_alba,
                                    pm = vcf_hybrid_alba_popmap,
                                    p1 = "A. nebrodensis",
                                    p2 = "A. alba",
                                    difference = 0.72)

ceph.vcfR.diff <- alleleFreqDiff(vcfR = vcf_hybrid_ceph,
                                 pm = vcf_hybrid_ceph_popmap,
                                 p1 = "A. nebrodensis",
                                 p2 = "A. cephalonica",
                                 difference = 0.2)

# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
hi.het.alb <- hybridIndex(vcfR = alba.vcfR.diff, pm = vcf_hybrid_alba_popmap, p1 = "A. nebrodensis", p2 = "A. alba")

hi.het.ceph <- hybridIndex(vcfR = ceph.vcfR.diff, pm = vcf_hybrid_ceph_popmap, p1 = "A. nebrodensis hybrids", p2 = "A. alba")


# Generate colors (or leave blank to use default)
#cols <- c("#762a83", "#1b7837", "#af8dc3", "#7fbf7b", "#bababa", "#878787")
cols <- c("#762a83", "#1b7837", "#878787", "#bababa")
# View triangle plot
triangulaR::triangle.plot(hi.het.alb, colors = cols, cex=3.5, jitter=0.025)


triangulaR::triangle.plot(hi.het.ceph, colors = cols, cex=3.5, jitter=0.025)


# gghybrid bayesian hybrid index----

# alba - neb
dat_seedlings <- read.data(file=here("data","hybrids","seedlings_alba_hybrids.str"),nprecol=2,NUMLOCI=18,NUMINDS=54,MISSINGVAL=NA)


prepdata<- data.prep(
  data=dat_seedlings$data,               #part of the read.data output object#
  loci=dat_seedlings$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_seedlings$alleles,         #part of the read.data output object#
  S0="neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="al",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
  precols=dat_seedlings$precols,         #part of the read.data output object#
  return.genotype.table=TRUE,  #This returns an table, with loci in columns, of diploid (or other ploidy) genotypes, coded as number of copies of the 
  #allele with higher frequency in S1 than S0#
  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values 
  #calculated within this function#
)


hindlabelalba<-esth(
  data.prep.object=prepdata$data.prep,
  read.data.precols=dat_seedlings$precols,
  include.Source=TRUE,	         #Leave at default TRUE if you want hybrid indices for the parental reference individuals#
  nitt=5000,                                                  #Testing suggests nitt=3000,burnin=1000 are plenty for accurate posterior estimation#
  burnin=2000
)


colorshyb<-brewer.pal(5,"Spectral")

layout(matrix(c(1,2),ncol=1))

# New Figure 3 
par(mar=c(4,4.1,1,1))
bc <- plot_h2(data=hindlabelalba$hi,
             test.subject=hindlabelalba$test.subject,
             mean.h.by="POPID",			                    #Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
             #individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
             col.group="POPID",
             group.sep=NULL,
             fill.source=F,
             basic.lines=T,
             labelx="",
             labely=substitute(paste("Hybrid index", italic(" A. alba"))),
             source.col=colorshyb[c(4,3)],
           #  source.limits=c("NA","NA"),
             #custom.abline=abline(h=0.5,col="black",lwd=1, lty=2),
             cex=1.5,pch=21,
             cex.lab=1.1,cex.axis=1.1,ylim=c(0,1))


# ceph - neb
dat_ceph <- read.data(file=here("data","hybrids","seedlings_ceph_hybrids.str"),nprecol=2,NUMLOCI=18,NUMINDS=54,MISSINGVAL=NA)


prepdata_ceph<- data.prep(
  data=dat_ceph$data,               #part of the read.data output object#
  loci=dat_ceph$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_ceph$alleles,         #part of the read.data output object#
  S0="neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="ceph",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
  precols=dat_ceph$precols,         #part of the read.data output object#
  return.genotype.table=TRUE,  #This returns an table, with loci in columns, of diploid (or other ploidy) genotypes, coded as number of copies of the 
  #allele with higher frequency in S1 than S0#
  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values 
  #calculated within this function#
)


hindlabelceph<-esth(
  data.prep.object=prepdata_ceph$data.prep,
  read.data.precols=dat_ceph$precols,
  include.Source=TRUE,	         #Leave at default TRUE if you want hybrid indices for the parental reference individuals#
  nitt=5000,                                                  #Testing suggests nitt=3000,burnin=1000 are plenty for accurate posterior estimation#
  burnin=2000
)


hceph <- plot_h2(data=hindlabelceph$hi,
              test.subject=hindlabelceph$test.subject,
              mean.h.by="POPID",	#Calculate the mean hybrid index for each value of the "POPID" column#
              sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
              #individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
              col.group="POPID",
              group.sep=NULL,
              fill.source=F,
              labely=substitute(paste("Hybrid index", italic(" A. cephalonica"))),
              basic.lines=T,
              source.col=colorshyb[c(4,2)],
              source.limits=c("NA","NA"),
              #custom.abline=abline(h=0.5,col="black",lwd=1, lty=2),
              cex=1.5,pch=21,
              cex.lab=1.1,cex.axis=1.1,ylim=c(0,1))
