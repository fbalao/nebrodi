##### AGF

# Load R libraries ----
library(adegenet)
library(poppr)
library(dartR)
library(hierfstat)
library(RColorBrewer)
library(plotrix)
library(strataG)
library(ggplot2)
library(forcats)
library(gridExtra)
library(ggpubr)
library(snpR)
library(here)
#source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
source(here("scripts","AGF","hybridize3.R"))



# Data format ----
## nebrodf3_118 and nebrodf3
#load(here("envs","nebrodi_29.11.2020.RData"))
#load(here("envs","AGP_results8.11.2023.RData"))

loci100<-locNames(nebrodf3)
nebrodf4<-nebrodf3_118
plantulasNDB<-nebrodf4[nebrodf4@other$Plantula_madre==1]



# Synthetic AGP population ----

nebrodf4<-nebrodf4[1:148]
nebrodf4<-nebrodf4[loc = loci100]

db4<-codes[match(rownames(nebrodf4@tab), codes$Codigo_muestra),]

nebrodf4@other$Plantula_madre<-as.factor(db4$Plantula_madre)
nebrodf4@other$Vivero_campo<-as.factor(db4$Vivero_campo)
nebrodf4@other$Madre_procedencia<-as.factor(db4$Madre_procedencia)


### loop

codespadres<-codes[codes$Plantula_madre==0 & codes$Vivero_campo==0,c(1,5)]


listacruces<-read.table(here("data","agf","12crosses_nohybs.txt"))
colnames(listacruces)<-c("M1","M2")




syntheticpop<-list()
for (i in 1: dim(listacruces)[1]){
  print(paste(listacruces$M1[i],listacruces$M2[i], sep="_"))
  syntheticpop[[i]]<-hybridize3(nebrodf4[listacruces$M1[i]],
                                nebrodf4[listacruces$M2[i]],
                                nebrodf4,
                                res.type="genind",quiet=T,n=10,
                                pop=paste(listacruces$M1[i],listacruces$M2[i], sep="_"))
}

syntallloci<-lapply(syntheticpop,locNames)
syntaredloci<-Reduce(intersect, syntallloci)


synt_red<-lapply(syntheticpop,function(x) x[loc=syntaredloci])


# seedlings outbreed from nursery
#### For next generation, outbreed plantulas from the nursery would be included

load(here("envs","nursery.RData"))

outbreed_nursery<-read.table(here("data","agf","outbreed_nursery.txt"), header=T)
nursery3<-nursery2[loc=syntaredloci]
nurseryagp<-nursery3[outbreed_nursery$ID]


# reduce parents' genind to loci and parentals in simulation  

pr<-padresDB[loc=syntaredloci]
pr<-pr[unique(c(listacruces[,1],listacruces[,2]))]


# reduced all db to loci in simulation

nebrodf4<-nebrodf4[loc=syntaredloci]
pop(nebrodf4)<-nebrodf4@other$Plantula_madre
pop(nebrodf4)<-rownames(nebrodf4@tab)

padrescruces<-nebrodf4[c(listacruces[,1],listacruces[,2])]
synthdb<-repool(synt_red)
pooled<-repool(synthdb, padrescruces)

# remove NA to run
pooled.tab<- missingno(pooled, type = "loci", cutoff = 0, quiet = FALSE, freq = FALSE)
pca.neb <- dudi.pca(pooled.tab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.neb$li, fac=pop(pooled.tab),
        col=rainbow(148), clabel = 0.9)


# combine synthetics and wild seedlings and nursery
pop(nebrodf4)<-nebrodf4$other$Plantula_madre

agp<-repool(synthdb,nebrodf4[pop=2],nurseryagp)
pop(agp)<-as.factor(rep("agp",nInd(agp)))


# Diversity indices with snpR  -------
# https://github.com/hemstrow/snpR

pop(padresDB)<-rep("adults",30)
pop(plantulasNDB)<-rep("seedlings",118)
pop(synthdb)<-rep("agp",120)
pop(agp)<-rep("nextg",257)

todo<-repool(padresDB,plantulasNDB,synthdb,agp)

padressnp<-import.snpR.data(padresDB)
plantulasnp<-import.snpR.data(plantulasNDB)
agpsnp<-import.snpR.data(synthdb)
nextgsnp<-import.snpR.data(agp)
todosnp<-import.snpR.data(todo)

sample.meta(padressnp)$pop<-rep("adults",30)
sample.meta(plantulasnp)$pop<-rep("seedlings",118)
sample.meta(agpsnp)$pop<-rep("agp",120)
sample.meta(nextgsnp)$pop<-rep("nextg",257)

padressnp.dat <- calc_maf(padressnp, "pop")
padressnp.dat <- calc_pi(padressnp.dat, "pop")
padressnp.dat <- calc_ho(padressnp.dat, "pop")
padressnp.dat <- calc_fis(padressnp.dat, "pop")
padressnp.dat <- calc_he(padressnp.dat, "pop")
padressnp.dat <- calc_hwe(padressnp.dat, "pop")


padres.div<-get.snpR.stats(padressnp.dat, "pop",stats =c("ho","he", "fis","pi", "hwe"))


plantulasnp.dat <- calc_maf(plantulasnp, "pop")
plantulasnp.dat <- calc_pi(plantulasnp.dat, "pop")
plantulasnp.dat <- calc_ho(plantulasnp.dat, "pop")
plantulasnp.dat <- calc_fis(plantulasnp.dat, "pop")
plantulasnp.dat <- calc_he(plantulasnp.dat, "pop")
plantulasnp.dat <- calc_hwe(plantulasnp.dat, "pop")


plantulas.div<-get.snpR.stats(plantulasnp.dat, "pop",stats =c("ho","he", "fis","pi", "hwe"))


agpsnp.dat <- calc_maf(agpsnp, "pop")
agpsnp.dat <- calc_pi(agpsnp.dat, "pop")
agpsnp.dat <- calc_ho(agpsnp.dat, "pop")
agpsnp.dat <- calc_fis(agpsnp.dat, "pop")
agpsnp.dat <- calc_he(agpsnp.dat, "pop")
agpsnp.dat <- calc_hwe(agpsnp.dat, "pop")

agp.div<-get.snpR.stats(agpsnp.dat, "pop",stats =c("ho","he", "fis","pi", "hwe"))


nextgsnp.dat <- calc_maf(nextgsnp, "pop")
nextgsnp.dat <- calc_pi(nextgsnp.dat, "pop")
nextgsnp.dat <- calc_ho(nextgsnp.dat, "pop")
nextgsnp.dat <- calc_fis(nextgsnp.dat, "pop")
nextgsnp.dat <- calc_he(nextgsnp.dat, "pop")
nextgsnp.dat <- calc_hwe(nextgsnp.dat, "pop")

nextg.div<-get.snpR.stats(nextgsnp.dat, "pop",stats =c("ho","he", "fis","pi", "hwe"))



todo.div<-rbind(padres.div$single,plantulas.div$single , agp.div$single, nextg.div$single )

p <- todo.div %>%
  mutate(subfacet = fct_relevel(subfacet, 
                                "adults", "seedlings", "agp", "nextg")) 



todosnp.dat <- calc_maf(todosnp, "pop")
todosnp.dat <- calc_pi(todosnp.dat, "pop")
todosnp.dat <- calc_ho(todosnp.dat, "pop")
todosnp.dat <- calc_fis(todosnp.dat, "pop")
todosnp.dat <- calc_he(todosnp.dat, "pop")
#todosnp.dat <- calc_ne(todosnp.dat, "pop")
todosnp.dat <- calc_hwe(todosnp.dat, "pop")
todosnp.dat <- calc_prop_poly(todosnp.dat, "pop")
todosnp.dat <- calc_private(todosnp.dat, "pop")



# Ne estimation among adults, seedlings, and synthetic pop ----

# snpR run in Windows 10

padressnp.dat <- calc_ne(padressnp.dat, "pop", 
                         NeEstimator_path="Ne2-1.exe")

plantulasnp.dat <- calc_ne(plantulasnp.dat, "pop", 
                           NeEstimator_path="Ne2-1.exe")


agpsnp.dat <- calc_ne(agpsnp.dat, "pop", 
                      NeEstimator_path="Ne2-1.exe")

nextgsnp.dat <- calc_ne(nextgsnp.dat, "pop", 
                        NeEstimator_path="Ne2-1.exe")

# result folder are moved to Results/NeEstimator
# merge results

Ne_results<- tibble(read.table("~/Datos/R/Rpackages/nebrodi/results/NeEstimator/Nesummary.txt", header=T))

Ne_results <- Ne_results %>%
  mutate(generation = fct_relevel(generation, 
                                  "adults", "seedlings", "agp", "nextg"))

# StrataG (not run in all populations) REVIEW
nebrodf4_gtypes <- genind2gtypes(padresDB)
Ne_nebrodf4 <- ldNe(nebrodf4_gtypes, maf.threshold = 0, by.strata = F, ci = 0.95, drop.missing = TRUE, num.cores = 15)


plantulas_gtypes <- genind2gtypes(nebrodf4[pop=2])
Ne_plantulas <- ldNe(plantulas_gtypes , maf.threshold = 0, by.strata = F, ci = 0.95, drop.missing = TRUE, num.cores = 15)


synthdb_gtypes <- genind2gtypes(synthdb)
Ne_synthdb <- ldNe(synthdb_gtypes, maf.threshold = 0, by.strata = F, ci = 0.95, drop.missing = TRUE, num.cores = 15)


agp2<-missingno(agp,cutoff = 0)
agp_gtypes <- genind2gtypes(agp2)
Ne_agp <- ldNe(agp_gtypes, maf.threshold = 0, by.strata = F, ci = 0.95, drop.missing = TRUE, num.cores = 15)

# Diversity comparisons ----

kruskal.test(p$pi~p$subfacet)
kruskal.test(p$ho~p$subfacet)
kruskal.test(p$fis~p$subfacet)
kruskal.test(p$fis~p$subfacet)

TukeyHSD(aov(p$ho~p$subfacet))
TukeyHSD(aov(p$fis~p$subfacet))


library(rcompanion)

epsilonSquared(x = p$ho,
               g = p$subfacet, ci=T, histogram=T, 
               type="perc", reportIncomplete=T)


epsilonSquared(x = p$fis,
               g = p$subfacet, ci=T, histogram=T, 
               type="perc", reportIncomplete=T)

epsilonSquared(x = p$pi,
               g = p$subfacet, ci=T, histogram=T, 
               type="perc", reportIncomplete=T)
library(rstatix)

p %>% kruskal_effsize(ho ~ subfacet)



# Diversity violin plot ----

my_comparisons <- list( c("adults", "seedlings"), c("seedlings", "agp"), 
                        c("adults", "agp"), c("adults", "nextg"), 
                        c("seedlings", "nextg"), c("nextg", "agp") )


dpho <-ggplot(p, aes(x=subfacet, y=ho, fill=subfacet)) + 
  geom_violin(trim=F, show.legend = FALSE) +
  scale_y_continuous(limits=c(0,1.3), breaks = c(0,0.25,0.50,0.75,1)) +
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Ho")+
  geom_boxplot(width=0.1, fill="white") 
# stat_compare_means(comparisons = my_comparisons, paired = T, label = "p.signif", position="jitter",
#                    label.y=c(1.05, 1.15 ,1.25, 1.30,1.15, 1.05))


dpho2<-dpho + scale_fill_brewer(palette="Spectral") + theme_classic() +
  theme(axis.text.x=element_blank(), axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

dpho2<-dpho2+ annotate(geom="text", x=c(1,2,3,4), y=c(0.95,0.85,1.30,0.75), label=c("a","b","c","a"),
                       color="black", size=4.5)


dpfis <-ggplot(p, aes(x=subfacet, y=fis, fill=subfacet)) + 
  geom_hline(yintercept = 0, col="grey15", linetype="dashed") +
  geom_violin(trim=FALSE, show.legend = FALSE) +
  scale_y_continuous(limits=c(-0.8,1.45), breaks = c(-0.5,0,0.5,1))+
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Fis")+
  geom_boxplot(width=0.1, fill="white") 
# stat_compare_means(comparisons = my_comparisons, paired = T, label = "p.signif", position="jitter",
#                    label.y=c(1.05, 1.15 ,1.25, 1.30,1.15, 1.05))


dpfis2<-dpfis + scale_fill_brewer(palette="Spectral") + theme_classic() +
  theme(axis.text.x=element_blank(), axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

dpfis2<-dpfis2+ annotate(geom="text", x=c(1,2,3,4), y=c(1.45,1.45,0.65,0.85), label=c("a","b","c","a"),
                         color="black", size=4.5)

dppi <-ggplot(p, aes(x=subfacet, y=pi, fill=subfacet)) + 
  geom_violin(trim=FALSE, show.legend = FALSE) +
  #ylim(-0.05,1) +
  labs(x="", y = expression(pi))+
  geom_boxplot(width=0.1, fill="white") 
# stat_compare_means(comparisons = my_comparisons, paired = T, label = "p.signif", position="jitter",
#                    label.y=c(1.05, 1.15 ,1.25, 1.30,1.15, 1.05))


dppi2<-dppi + scale_fill_brewer(palette="Spectral") + theme_classic() + 
  theme(axis.text.x=element_blank(),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



dpne<-ggplot(Ne_results, aes(generation, Ne)) +
  geom_pointrange(
    aes(ymin = lCI/1.5, ymax = uCI*1.5, color =generation), show.legend=F,
    position = position_dodge(1), linewidth = 1.3, size=0.8) +
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Ne")

dpne2<-dpne+ scale_color_brewer(palette="Spectral")+theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


grid.arrange(dpho2,dpfis2, dppi2, dpne2, ncol = 1)



# Diversity comparison

Ne_results

# Diversity bar plot ----

todo.summary <- todo.div %>%
  group_by(subfacet) %>%
  summarise(
    se.ho = std.error(ho, na.rm = TRUE),
    ho = mean(ho),
    se.fis = std.error(fis, na.rm = TRUE),
    fis = mean(fis),
    se.pi = std.error(pi, na.rm = TRUE),
    pi = mean(pi)
  )

todo.summary <- todo.summary %>%
  mutate(subfacet = fct_relevel(subfacet, 
                                "adults", "seedlings", "agp", "nextg"))

bpho<-ggplot(todo.summary, aes(subfacet, ho, fill=subfacet)) +
  geom_bar(stat="identity", color="black", show.legend = FALSE) +
  coord_cartesian(ylim = c(0.05, 0.6)) +
  geom_errorbar(aes(ymin = ho, ymax = ho+se.ho), width = 0.2) +
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Ho")

bpho2<-bpho  + scale_fill_brewer(palette="Spectral") + theme_classic()  + theme(axis.text.x=element_blank())


bpfis<-ggplot(todo.summary, aes(subfacet, fis, fill=subfacet)) +
  geom_hline(yintercept = 0, col="grey15", linetype="dashed") +
  geom_bar(stat="identity", color="black", show.legend = FALSE) +
  coord_cartesian(ylim = c(-0.20, 0.5)) +
  geom_errorbar(aes(ymin = fis, ymax = fis+se.fis), width = 0.2) +
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Fis")

bpfis2<-bpfis  + scale_fill_brewer(palette="Spectral") + theme_classic()  + theme(axis.text.x=element_blank())


bppi<-ggplot(todo.summary, aes(subfacet, pi, fill=subfacet)) +
  geom_bar(stat="identity", color="black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = pi, ymax = pi+se.pi), width = 0.2) +
  scale_x_discrete(labels=c("Adults", "Seedlings", "AGF-Pop", "Next-Gen")) +
  labs(x="", y = "Nucleotide diversity")

bppi2<-bppi  + scale_fill_brewer(palette="Spectral") + theme_classic()

grid.arrange(bpho2,bpfis2, bppi2, ncol = 1)



# Save Results to RData ----
save.image("AGP_results8.11.2023.RData")

