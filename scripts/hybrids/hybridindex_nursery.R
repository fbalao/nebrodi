# Hybrid index seedlings from nursery

# Libraries ----
library(gghybrid)
library(coda)
library(here)
library(vcfR)
library(RColorBrewer)
library(radiator)
library(adegenet)
source(here("scripts","Diversity","genind2structrure.R"))
source(here("scripts","hybrids","plot_h2.R"))


#  we need nursery genotypes from Analisis_OpenArray_nebrodi_nursery_seedlings.R
# we need abies exotic species genotypes (abiesintrogress_genind) and A nebrodi parentals
# (padresDB_sub) from hybridsPCA.R
# adegenet nursery2 to structure format. select hybrids and snps ----
load(here("envs","nursery.RData"))


nursery3<-nursery2[loc=locNames(abiesintrogress_genind)]

nursery3@pop<-factor(rep("seedling", times=1777), levels=c("seedling",
                        "hybridneb", "ospecies","unknown"))


nursery3@pop<-factor(rep("seedling", times=1777))

selected_seedlings<-read.table(here("data","hybrids","nurseryhybrids.txt"), header=T)

nursery4<-nursery3[selected_seedlings$ID]

nursery4@pop<-as.factor(selected_seedlings$category)
nursery4<-nursery4[indNames(nursery4)!="11_2008_0426"] # This individual has no information for the hybrid loci

poolabies_nursery<-repool(abiesintrogress_genind, padresDB_sub,nursery4)


genind2structure(poolabies_nursery, file=here("data","hybrids","nurseryhybrids.str"), pops=T)


# # no hybrid seedlings
# selected_seedlings$ID <- as.character(selected_seedlings$ID)
# nursery_nhyb<-nursery3[!indNames(nursery3) %in% selected_seedlings$ID]
# genind2structure(nursery_nhyb, file=here("data","hybrids","nursery_nohybrids.str"), pops=T)




# gghybrid bayesian hybrid index----

# alba - neb
dat_nurseryalba <- read.data(file=here("data","hybrids","nursery_alba_hybrids.str"),nprecol=2,NUMLOCI=21,NUMINDS=923,MISSINGVAL=NA)


prepdata<- data.prep(
  data=dat_nurseryalba$data,               #part of the read.data output object#
  loci=dat_nurseryalba$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_nurseryalba$alleles,         #part of the read.data output object#
  S0="1_neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="3_alb",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
  precols=dat_nurseryalba$precols,         #part of the read.data output object#
  return.genotype.table=TRUE,  #This returns an table, with loci in columns, of diploid (or other ploidy) genotypes, coded as number of copies of the 
  #allele with higher frequency in S1 than S0#
  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values 
  #calculated within this function#
)


hindlabelalba<-esth(
  data.prep.object=prepdata$data.prep,
  read.data.precols=dat_nurseryalba$precols,
  include.Source=TRUE,	         #Leave at default TRUE if you want hybrid indices for the parental reference individuals#
  nitt=5000,                                                  #Testing suggests nitt=3000,burnin=1000 are plenty for accurate posterior estimation#
  burnin=2000
)

tapply(hindlabelalba$hi$h_posterior_mode, hindlabelalba$hi$Source, median)

hybridsalba<-hindlabelalba$hi[hindlabelalba$hi$h_posterior_mode>0.25]

hybridsalba[hybridsalba$POPID=="2_seedlings_nursery",]

colorshyb<-brewer.pal(5,"Spectral")

layout(matrix(c(1,2),ncol=1))

# New Figure S4  
par(mar=c(4,4.1,1,1))
bc <- plot_h2(data=hindlabelalba$hi,
              test.subject=hindlabelalba$test.subject,
              mean.h.by="POPID",			                    #Calculate the mean hybrid index for each value of the "POPID" column#
              sort.by=c("POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
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
              cex=1,pch=21,
              cex.lab=1.1,cex.axis=1.1,ylim=c(0,1))


# ceph - neb
dat_nurseryceph <- read.data(file=here("data","hybrids","nursery_ceph_hybrids.str"),nprecol=2,NUMLOCI=21,NUMINDS=923,MISSINGVAL=NA)


prepdata_ceph<- data.prep(
  data=dat_nurseryceph$data,               #part of the read.data output object#
  loci=dat_nurseryceph$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_nurseryceph$alleles,         #part of the read.data output object#
  S0="1_neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="3_cep",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
  precols=dat_nurseryceph$precols,         #part of the read.data output object#
  return.genotype.table=TRUE,  #This returns an table, with loci in columns, of diploid (or other ploidy) genotypes, coded as number of copies of the 
  #allele with higher frequency in S1 than S0#
  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values 
  #calculated within this function#
)


hindlabelceph<-esth(
  data.prep.object=prepdata_ceph$data.prep,
  read.data.precols=dat_nurseryceph$precols,
  include.Source=TRUE,	         #Leave at default TRUE if you want hybrid indices for the parental reference individuals#
  nitt=5000,                                                  #Testing suggests nitt=3000,burnin=1000 are plenty for accurate posterior estimation#
  burnin=2000
)



hceph <- plot_h2(data=hindlabelceph$hi,
                 test.subject=hindlabelceph$test.subject,
                 mean.h.by="POPID",	#Calculate the mean hybrid index for each value of the "POPID" column#
                 sort.by=c("POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
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




# summary

hindex_seedlingceph<-hindlabelceph$hi[hindlabelceph$hi$POPID=="2_seedlings_nursery",]

names_hybrids_neb_ceph<-hindex_seedlingceph[hindex_seedlingceph$h_posterior_mode>0.25,2]$INDLABEL

hybrids_neb_ceph<-hindex_seedlingceph[hindex_seedlingceph$h_posterior_mode!=0,]


selected_seedlings[which(selected_seedlings$ID%in%names_hybrids_neb_ceph),]
