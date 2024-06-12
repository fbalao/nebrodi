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

# # Para alba
# 6:     S0          31M           1_neb             0.87             0.44
# 17:     S0          21M           1_neb             0.54             0.20
# 18:   TEST       20.1.P 2_neb_seedlings             0.33             0.10
# 19:     S0           1M           1_neb             0.23             0.05
# 20:     S0           9M           1_neb             0.14             0.03
# 21:   TEST      18.10.P 2_neb_seedlings             0.13             0.02
# 22:   TEST      22.18.P 2_neb_seedlings             0.13             0.02



# # Para ceph
# 1:     S0          31M           1_neb             1.00             0.69 
# TEST       20.1.P 2_neb_seedlings             0.78             0.37             0.95 0.70699341 0.023669543
# 18:     S0          21M           1_neb             0.49             0.17             0.82 0.49218511 0.028899046
# 19:     S0          26M           1_neb             0.47             0.16             0.81 0.47750202 0.029629240
# # 20:   TEST      22.18.P 2_neb_seedlings             0.42             0.13             0.78 0.44198390 0.029916819
# 1:   TEST      18.10.P 2_neb_seedlings             0.21             0.04             0.67 0.29837468 0.028032072
# 22:   TEST      18.13.P 2_neb_seedlings             0.14             0.02             0.68 0.27739999 0.032138874
# 23:   TEST       21.5.P 2_neb_seedlings             0.13             0.02             0.68 0.27053319 0.031651412