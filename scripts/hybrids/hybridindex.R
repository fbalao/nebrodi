# Libraries ----
library(gghybrid)
library(coda)
library(here)
source(here("scripts","hybrids","plot_h2.R"))
source(here("scripts","Diversity","genind2structure.R"))

#  transform genind to structure format ----
# poolabies2 is coming from hybridsPCA.R
genind2structure(poolabies2, file=here("data","hybrids","seedlings_hybrids.str"),
                 pops=T)
#modify the structure file: change pop names, remove cephalonica/alba and missiing data code

# gghybrid bayesian hybrid index----

# alba - neb
dat_seedlings <- read.data(file=here("data","hybrids","seedlings_alba_hybrids_new2.str"),nprecol=2,NUMLOCI=21,NUMINDS=54,MISSINGVAL=NA)


prepdata<- data.prep(
  data=dat_seedlings$data,               #part of the read.data output object#
  loci=dat_seedlings$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_seedlings$alleles,         #part of the read.data output object#
  S0="1_neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="3_alb",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
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



# New Figure 3 
layout(matrix(c(1,1,2,3),ncol=2))
par(mar=c(4,4.1,1,1), pty="m")

# PCA from hybridsPCA.R
# col1<-c("#FDAE61","#FFFFBF", "#ABDDA4", "#e58f90")
# plot(x,y, bg=bgcol, pch=21, cex=2, col=circ,
#      xlab="PC1 (29.6%)", ylab="PC2 (14.1%)", cex.axis=1.3, cex.lab=1.4)
# 
# legend(1,-3.5,legend = c(expression(italic("A. alba")),
#                          expression(italic("A. cephalonica")),
#                          expression(paste(italic("A. nebrodensis")," adults")),
#                          expression(paste(italic("A. nebrodensis"), " seedlings"))),
#        pch=21, cex=1.4, col="grey35", pt.bg=col1, pt.cex=1.6,
#        box.lwd = 0, border = "white", bg="transparent")

# hybrid index plots 
par(mar=c(4,7,1,1), pty="m")

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
             source.col=colorshyb[c(4,2)],
           #  source.limits=c("NA","NA"),
             #custom.abline=abline(h=0.5,col="black",lwd=1, lty=2),
             cex=1.5,pch=21,
             cex.lab=1.1,cex.axis=1.1,ylim=c(0,1))


# ceph - neb
dat_ceph <- read.data(file=here("data","hybrids","seedlings_ceph_hybrids_new.str"),nprecol=2,NUMLOCI=21,NUMINDS=54,MISSINGVAL=NA)


prepdata_ceph<- data.prep(
  data=dat_ceph$data,               #part of the read.data output object#
  loci=dat_ceph$loci,               #part of the read.data output object#
  sourceAbsent = FALSE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
  alleles=dat_ceph$alleles,         #part of the read.data output object#
  S0="1_neb",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
  S1="3_cep",                        #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
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
              sort.by=c("POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
              #individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
              col.group="POPID",
              group.sep=NULL,
              fill.source=F,
              labely=substitute(paste("Hybrid index", italic(" A. cephalonica"))),
              basic.lines=T,
              source.col=colorshyb[c(4,3)],
              source.limits=c("NA","NA"),
              #custom.abline=abline(h=0.5,col="black",lwd=1, lty=2),
              cex=1.5,pch=21,
              cex.lab=1.1,cex.axis=1.1,ylim=c(0,1))
