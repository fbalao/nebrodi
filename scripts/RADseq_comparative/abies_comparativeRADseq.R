# Libraries ----
# RAD Pop diversity poplevel2 4ind - nebrodi 25.14.3
library(RColorBrewer)
library(gridExtra)
library(lemon)
library(ggplot2)
library(here)
library(tidyverse)
library(vcfR)
library(snpR)
library(dartR)
library(here)



##### All diversity statistics

datarad<-read.table("data/diversity/populations.sumstats_summary.tsv", header=T)

datarad$sp<-as.factor(c("nebrodensis","cephalonica","cephalonica","cephalonica","alba","alba","alba"))


dataradtest<-rbind(datarad, datarad[1,], datarad[1,])
#levels(dataradtest$sp)<-c("A. alba", "A. cephalonica", "A. nebrodensis")
levels(dataradtest$sp)<-c("Silver fir", "Greek fir", "Sicilian fir")


legend<- ggplot(dataradtest, aes(x=sp, y=Pi, fill=sp)) +
  geom_violin(trim=FALSE, show.legend = T) +
  scale_fill_brewer(palette="Spectral") +
theme(legend.title=element_blank(), legend.text = element_text(size=14,face = "italic")) 


legendp <- g_legend(legend + theme(legend.position='bottom'))
  

dpi <-ggplot(datarad, aes(x=sp, y=Pi, fill=sp)) + 
  geom_violin(trim=FALSE, show.legend = F) +
  geom_point(data = filter(datarad, sp == "nebrodensis"), size=6, 
             color=brewer.pal(n = 4, name = "Spectral")[3], show.legend = F) +
  ylim(0.185,0.245) +
  scale_x_discrete(labels=expression("A. alba", "A. cephalonica", "A. nebrodensis")) +
    labs(x="", y = expression(Pi))+
  geom_boxplot(width=0.1, fill="white")  +
  scale_fill_brewer(palette="Spectral", direction=1) + theme_classic() + 
  theme(axis.text.x=element_blank(),axis.text=element_text(size=12),
   axis.title=element_text(size=14,face="bold"))
dpi

dpriv <-ggplot(datarad, aes(x=sp, y=Private, fill=sp)) + 
  geom_violin(trim=FALSE, show.legend = T) +
  geom_point(data = filter(datarad, sp == "nebrodensis"), size=6, 
             color=brewer.pal(n = 4, name = "Spectral")[3], show.legend = F) +
  #ylim(0.185,0.245) +
  scale_x_discrete(labels=expression("A. alba", "A. cephalonica", "A. nebrodensis")) +
  labs(x="", y = "Private alleles")+
  geom_boxplot(width=0.1, fill="white")  +
  scale_fill_brewer(palette="Spectral") + theme_classic() + 
  theme(axis.text.x=element_blank(),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

dpriv


dfis <-ggplot(datarad, aes(x=sp, y=Fis, fill=sp)) + 
  geom_violin(trim=FALSE, show.legend = F) +
  geom_point(data = filter(datarad, sp == "nebrodensis"), size=6, 
             color=brewer.pal(n = 4, name = "Spectral")[3], show.legend = F) +
    #ylim(0.185,0.245) +
  scale_x_discrete(labels=expression("A. alba", "A. cephalonica", "A. nebrodensis")) +
  labs(x="", y = "Fis")+
  geom_boxplot(width=0.1, fill="white")  +
  scale_fill_brewer(palette="Spectral") + theme_classic() + 
  theme(axis.text.x=element_blank(),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dfis

dhet <-ggplot(datarad, aes(x=sp, y=Exp_Het, fill=sp)) + 
  geom_violin(trim=FALSE, show.legend = F) +
  geom_point(data = filter(datarad, sp == "nebrodensis"), size=6, 
             color=brewer.pal(n = 4, name = "Spectral")[3], show.legend = F) +
   #ylim(0.185,0.245) +
  scale_x_discrete(labels=expression("A. alba", "A. cephalonica", "A. nebrodensis")) +
  labs(x="", y = "Hs")+
  geom_boxplot(width=0.1, fill="white")  +
  scale_fill_brewer(palette="Spectral") + theme_classic() + 
  theme(axis.text.x=element_blank(),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dhet


grid.arrange(dpi+theme(legend.position='hidden'), dhet+theme(legend.position='hidden'),
             dpriv+theme(legend.position='hidden'),dfis+theme(legend.position='hidden'),
             bottom=legendp$grobs[[1]], layout_matrix=matrix(c(1,3,2,4), ncol=2))


### analysis

# kruskal-wallis test

kruskal.test(datarad$Pi, datarad$sp)
kruskal.test(datarad$Private, datarad$sp)
kruskal.test(datarad$Exp_Het, datarad$sp)
kruskal.test(datarad$Fis, datarad$sp)

# kruskal-wallis effect size

library(rstatix)

datarad %>% kruskal_effsize(Pi ~ sp)
datarad %>% kruskal_effsize(Private ~ sp)
datarad %>% kruskal_effsize(Exp_Het ~ sp)
datarad %>% kruskal_effsize(Fis ~ sp)

library(rcompanion)

epsilonSquared(x = datarad$Pi,
               g = datarad$sp, ci=T, histogram=T, 
               type="perc", reportIncomplete=T)
epsilonSquared(x = datarad$Private,
               g = datarad$sp, ci=T, histogram=T, 
               type="perc", reportIncomplete=T, B=100000)
epsilonSquared(x = datarad$Exp_Het,
               g = datarad$sp)
epsilonSquared(x = datarad$Fis,
               g = datarad$sp)

# Genomic differentiation ----
fst_genomewide<-read.table("Abies_nebrodi/genomic_analysis/Aneb_populations_all/populations.phistats_nebrodensis-alba1.tsv", header=T)
fst_genomewide$LocusID<-1:length(fst_genomewide$LocusID)

fst_genomewide<-fst_genomewide[1:500,]
p <- ggplot(data = fst_genomewide, aes(x = LocusID, y = -log10(SmoothedPhi_stPvalue))) +
  geom_bar(stat = "identity") + labs(x = "Genome", y = "Fst")

p
# Ne estimation from vcf for 5 ind nebrodensis ----

vcf <- read.vcfR(here("Abies_nebrodi","genomic_analysis","Ne_5ind_RAD", "nebrodi_snps_final.vcf"))
chrom <- create.chromR(name='RAD_data', vcf=vcf)
nebrodisnp<-import.snpR.data(vcf)
sample.meta(nebrodisnp)
sample.meta(nebrodisnp)$pop<-c(rep("nebrodi",5), rep("cephalonica",15), rep("alba",15))

nebrodisnp.dat <- calc_maf(nebrodisnp, "pop")
nebrodisnp.dat <- calc_pi(nebrodisnp.dat, "pop")
nebrodisnp.dat <- calc_ho(nebrodisnp.dat, "pop")
nebrodisnp.dat <- calc_fis(nebrodisnp.dat, "pop")
nebrodisnp.dat <- calc_he(nebrodisnp.dat, "pop")
nebrodisnp.dat <- calc_hwe(nebrodisnp.dat, "pop")
nebrodisnp.dat <- calc_ne(nebrodisnp.dat, "pop",
                          NeEstimator_path="Ne2-1.exe")

nebrodi.div<-get.snpR.stats(nebrodisnp.dat, "pop",stats =c("ho","he", "fis","pi",
                                                      "prop_poly","private", "hwe"))

nebro.gl<-gl.read.vcf(here("Abies_nebrodi","genomic_analysis","Ne_5ind_RAD", "nebrodi_snps_pruned.vcf"), verbose = NULL)
nebro.gl$pop<-as.factor(rep("nebrodi",5))
gl.report.heterozygosity(nebro.gl)
nebro2.gl <- gl.edit.recode.ind(nebro.gl)
ne_neb<-gl.LDNe(nebro2.gl,neest.path ="/home/fbalao/Programas/Ne_estimator/", critical = c(0,0.05,0.1))


