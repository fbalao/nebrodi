# Libraries ---------------------------------------------------------------
library(ape)
library(related)
library(corrplot)
library(dendextend)
library(geodist)
library(adegenet)
library(ecodist)
library(ade4)
library(pegas)
library(poppr)
library(hierfstat)
library(ade4)
library(factoextra)
library(dartR)
library(corrplot)
library(RColorBrewer)
library(related)
library(kableExtra)
library(igraph)
library(dplyr)
library(ggsci)
source("GENHETv3.1.R")
# --------#############

db3<-codes[match(rownames(padresDB@tab), codes$Codigo_muestra),]

padresDB@other$Plantula_madre<-as.factor(db3$Plantula_madre)
padresDB@other$Vivero_campo<-as.factor(db3$Vivero_campo)
padresDB@other$Madre_procedencia<-as.factor(db3$Madre_procedencia)

padresDB@pop<-as.factor(rep("padres",30))
#diversidad

div<-summary(padresDB)

basicstat <- basic.stats(padresDB, diploid = TRUE, digits = 2) 


mlg(padresDB)
#Prueba del Equilibrio de Hardy-Weinberg

hw<-hw.test(padresDB, B =1000)
out<-hw[,3]<0.0001
plot(div$Hobs, div$Hexp, pch=16, col=as.numeric(out)+1)
abline(0,1)


boot.ppfis(padresDB)


x <-indpca(padresDB, ind.labels=padresDB@other$Madre_procedencia) 
plot(x, cex =1.2, col=rainbow(30)[padresDB@other$Madre_procedencia])



nebrotab<-tab(padresDB, freq=T, NA.method="mean")
pca.neb <- dudi.pca(nebrotab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.neb$li, fac=as.factor(padresDB@other$Madre_procedencia),
        col=rainbow(30), clabel = 0.9)

rownames(pca.neb$li)<-padresDB@other$Madre_procedencia

fviz_pca_ind(pca.neb,
col.ind = "x",
labelsize=6,
label="all",# Color by the quality of representation
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
title="Abies nebrodensis",
repel = TRUE)

distgen<- dist(padresDB, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)    
plot.phylo(bionjs(distgen), type="u", cex=0.7,show.tip.label = F )
tiplabels(padresDB@other$Madre_procedencia)

padresDBgl<-gi2gl(padresDB)

pc <- gl.pcoa(padresDBgl, nfactors=5)
gl.pcoa.plot(pc, padresDBgl, labels="ind", xaxis=1, yaxis=2)


padresDBloci <-genind2loci(padresDB)


hist(1-poppr::bitwise.dist(padresDB, scale_missing=T), main="", xlab="Genotypic Similarity", 
     freq=F, ylim=c(0,10), col=rgb(0,0,1,1/5))
hist(1-poppr::bitwise.dist(plantulasNDB, scale_missing=T), add=T, col=rgb(1,0,0,1/5), freq=F)

legend(x=0.75, y=9, legend=c("adults", "seedlings"),
       col=c(rgb(0,0,1,1/5),rgb(1,0,0,1/5)), pch=16,box.col = "white",
          cex=1.1)

multimode::locmodes(1-poppr::bitwise.dist(plantulasNDB, scale_missing=T), mod0=2,display=TRUE)

multimode::locmodes(1-poppr::bitwise.dist(padresDB, scale_missing=T), mod0=1,display=TRUE)


# Genotyping distances ----------------------------------------------------



##### Comparando G similaridad entre adultos y plantulas
colnpg<-pal_npg("nrc", alpha = 0.5)(2)

Gsim_padres<-1-poppr::bitwise.dist(padresDB, scale_missing=T)
Gsim_seedl<-1-poppr::bitwise.dist(plantulasNDB, scale_missing=T)

hist(Gsim_padres, main="", xlab="Genotypic Similarity", 
     ylab="Frequency Density", xlim=c(0.4,1),
     freq=F, ylim=c(0,10), col=colnpg[2])
hist(Gsim_seedl, add=T, col=colnpg[1], freq=F)

legend(x=0.75, y=9, legend=c("adults", "seedlings"),
       col=c(colnpg[2],colnpg[1]), pch=16,box.col = "white",
       cex=1.1, pt.cex = 2.4)

text(0.8,9.5,labels = "KS-test's D = 0.17  p < 0.001" )

multimode::locmodes(1-poppr::bitwise.dist(plantulasNDB, scale_missing=T), mod0=2,display=TRUE)

multimode::locmodes(1-poppr::bitwise.dist(padresDB, scale_missing=T), mod0=1,display=TRUE)



ks.test(Gsim_seedl,Gsim_padres)


# Identity by descendent

padresDBgl<-dartR::gi2gl(padresDB)
ibd_padres<-dartR::gl.grm(padresDBgl)
hist(as.dist(ibd_padres, upper=F, diag=F))


plantulasDBgl<-dartR::gi2gl(plantulasNDB)
ibd_plantulas<-dartR::gl.grm(plantulasDBgl)
hist(as.dist(ibd_plantulas, upper=F, diag=F))
hist(ibd_plantulas)


distgenDIFF <- dist.gene(padresDBloci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
plot.phylo(bionj(distgenDIFF), type="p", cex=0.7, show.tip.label = F)
tiplabels(padresDB@other$Madre_procedencia)

distgenDISS <- diss.dist(padresDB, percent = T, mat = T)
plot(bionj(distgenDISS), type="p", cex=0.7, show.tip.label=F)
tiplabels(padresDB@other$Madre_procedencia)

padresnames<-read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/diversity/muestras2.csv", header=T)
rownames(distgenDISS)

colnames(distgenDISS)<-rownames(distgenDISS)
  padresnames[match(rownames(distgenDISS), padresnames$Codigo_muestra), "ID_LIFE"]
colnames(distgenDISS)<-padresnames[match(rownames(distgenDISS), padresnames$Codigo_muestra), "ID_LIFE"]

ps<-propShared(padresDB)
rownames(ps)<-padresnames[match(rownames(ps), padresnames$Codigo_muestra), "ID_LIFE"]
colnames(ps)<-padresnames[match(rownames(ps), padresnames$Codigo_muestra), "ID_LIFE"]


corrplot(1-distgenDISS,  order = "hclust", hclust.method="ward.D2",
         addrect = 3,
         is.corr = F,
         cl.lim = c(0,1),
         col=col4(100),
         win.asp = 1)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F","#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061"))
#####3 fasstructure

gl2faststructure(padresDBgl,outfile="nebrodesispadre.fs", outpath = "/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/diversity",probar =FALSE)
rep(padresDBgl$ind.names, each=2)
write.table( rep(padresDBgl$ind.names, each=2),"/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/diversity/padresnames.txt")



# Adults Diversity ----

padreshet<-read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/colony/nebrodi118plantulas/nebrodi_parents118.txt")
locusname<-colnames(read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/colony/nebrodi118plantulas/nebrodi_markers118.txt", header=T))
ind_nebroDiv<-GENHET(padreshet,estimfreq="T",locname=locusname)
ind_nebroDiv<- as.data.frame(ind_nebroDiv)
ind_nebroDiv2<- as.data.frame(apply(ind_nebroDiv, 2,FUN=as.numeric))
ind_nebroDiv2$sampleid<-ind_nebroDiv$sampleid

ind_nebroDiv2$sampleid<-padresnames[match(ind_nebroDiv2$sampleid, padresnames$Codigo_muestra), "ID_LIFE"]


ind_nebroDiv2<-ind_nebroDiv2[order(ind_nebroDiv2$sampleid, decreasing = F),]

ind_nebroDiv2$ibr<-ibr$LH

plot(ind_nebroDiv2$PHt~ind_nebroDiv2$ibr)

# Clustering Adults DAPC Figure S2 ----

kneb<-find.clusters(padresDB) # Usando n.pca = 30,n.clust = 3)

neb_dapc<-dapc(padresDB,pop = kneb$grp, n.pca = 12, n.da = 2)
# se usan 13 PC tras mirar 
# neb_dapc<-dapc(padresDB,pop = kneb$grp, n.pca = 13, n.da = 2)
# temp <- optim.a.score(neb_dapc, n.sim=1000, n=30)
rownames(neb_dapc$tab)<-padresnames[match(rownames(neb_dapc$tab), padresnames$Codigo_muestra), "ID_LIFE"]
names(neb_dapc$grp)<-padresnames[match(names(neb_dapc$grp), padresnames$Codigo_muestra), "ID_LIFE"]

v<-names(sort(neb_dapc$grp))

v2<-indorder

colfunk<-funky(6)

compoplot(neb_dapc,border=T,
          col=c(colfunk[1], colfunk[3],colfunk[2]),show.lab=T, legend=F,
          subset=v2)








# Ritland coancestry Figure S2 ----

reldata <- readgenotypedata("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/colony/nebrodi118plantulas/nebrodi_parents118.txt")
reldata$gdata$V1<-padresnames[match(reldata$gdata$V1, padresnames$Codigo_muestra), "ID_LIFE"]

reldata$gdata
rxy<-coancestry(reldata$gdata, ritland = 2, trioml=2, allow.inbreeding =T, error.rates = 0.01)
rela<-rxy$relatedness
rela<-rela[order(rela$ritland,decreasing = F),]

ibr<-rxy$inbreeding
ibr<-ibr[order(ibr$LH),]

barplot(ibr$LH, names.arg = ibr$ind.id, ylim=c(-0.2,0.88), las=2, 
        space = 0.5, ylab="Inbreeding", col=col4(30))

barplot(rela$ritland, names.arg = rela$group, ylim=c(-1,1), las=2, 
        space = 0.5, ylab="Rxy", col=col4(30))



convert_to_distance_matrix <- function(df) {
  colnames(df) <- c("id1", "id2", "distance")
  all_ids <- unique(c(df$id1, df$id2))
  distance_matrix <- matrix(NA, nrow = length(all_ids), ncol = length(all_ids))
  rownames(distance_matrix) <- colnames(distance_matrix) <- all_ids
  
  for (i in 1:nrow(df)) {
    distance_matrix[df[i, "id1"], df[i, "id2"]] <- df[i, "distance"]
    distance_matrix[df[i, "id2"], df[i, "id1"]] <- df[i, "distance"]
  }
  
  return(distance_matrix)
}

dfr<-convert_to_distance_matrix(rela[c(2,3,9)])

col3 <- colorRampPalette(c("red", "white", "blue")) 

# Figure S2
layout(matrix(c(1,2,2),nrow=3))

compoplot(neb_dapc,border=T,
          col=c(colfunk[3], colfunk[2],colfunk[1]),show.lab=T, legend=F,
          subset=v2)

c<-corrplot(dfr,  order = "hclust", hclust.method = "ward.D2",
         addrect = 3,
         is.corr = F,
         # cl.lim = c(-1,1),
         col = COL2('RdBu', 100),
         #col=col3(100),
         win.asp = 1,
         diag = F,
         cl.ratio=0.2,
         tl.col = funky(6)[c(rep(2,12),rep(3,10),rep(1,8))])

indorder<-c$corrPos$yName[1:30]
indorder<-indorder[c(1,30,2:29)]



# Graph plot ----
# No funciona

grafo<-gengraph(padresDB,cutoff = 75,ngrp = 1)
set.seed(3)
plot.igraph(grafo[[1]],vertex.color = funky(6)[neb_dapc$grp],
            edge.arrow.size=0, 
            vertex.label.cex=0.7, 
            vertex.label.family="Helvetica",
            vertex.label.font=2,
            vertex.label.dist=0,
            vertex.shape="circle", 
            vertex.size=8, 
            vertex.label.color="grey15", 
            edge.label.cex=0.01,
            edge.label.color="white",
            edge.width=0.5,
            edge.color="grey85",
            asp=0.5,
            #  ylim=c(-0.1,0.8),xlim=c(0,0.75),
            margin=-0,
            rescale=T)






pDist<-pairDistPlot(padresDB, neb_dapc$grp)
pDist$vi

# Tabla diversidad ----

ind_nebroDiv2<-ind_nebroDiv2[order(ind_nebroDiv2$PHt, decreasing = T),]
ind_nebroDiv2[,2:7]<-round(ind_nebroDiv2[,2:7],digits = 2)
rownames(ind_nebroDiv2)<-ind_nebroDiv2$sampleid
ind_nebroDiv2<-ind_nebroDiv2[,-1]

colnames(ind_nebroDiv2)[6]<-"INBR"
ind_nebroDiv2[,-2] %>%
  kbl(align='ccccccc') %>%
  kable_paper(full_width = F) %>%
  #row_spec(0, angle = -15) %>%
  kable_styling(bootstrap_options = "condensed", font_size = 10, position = "center") %>%
  column_spec(2, color = "white",background = spec_color(ind_nebroDiv2$PHt, option="C", direction=1)) %>%
  column_spec(3, color = "white",background = spec_color(ind_nebroDiv2$Hs_exp, option="C", direction=1)) %>%
  column_spec(4, color = "white",background = spec_color(ind_nebroDiv2$IR, option="C", direction=1)) %>%
  column_spec(5, color = "white",background = spec_color(ind_nebroDiv2$HL, option="C", direction=1)) %>%
  column_spec(6, color = "white",background = spec_color(ind_nebroDiv2$INBR, option="C", direction=1)) %>%
  footnote(general_title = "", general = c("PHt : proportion of heterozygous loci (PHt) in an individual.",
          "Hs_exp : standardized heterozygosity based on the mean expected heterozygosity (Coltman 1999).",
"IR : internal relatedness (Amos 2001).",
"HL: homozygosity by locus (Aparicio 2006).",
"INBR: inbreeding cofficient.")) %>%
 # as_image(width = 8) %>%
  save_kable("test.pdf", density=1200)



############Graf

layout(matrix(c(1,2),nrow=2))
compoplot(neb_dapc,border=T,
          col=funky(6),show.lab=T, legend=F,
          subset=v)
text(0,1.4,labels = "A)", font = 2, cex=1.2)
set.seed(3)
plot.igraph(grafo[[1]],vertex.color = funky(6)[neb_dapc$grp],
            edge.arrow.size=0, 
            vertex.label.cex=0.7, 
            vertex.label.family="Helvetica",
            vertex.label.font=2,
            vertex.label.dist=0,
            vertex.shape="circle", 
            vertex.size=8, 
            vertex.label.color="grey15", 
            edge.label.cex=0.01,
            edge.label.color="white",
            edge.width=0.5,
            edge.color="grey85",
            asp=0.5,
            #  ylim=c(-0.1,0.8),xlim=c(0,0.75),
            margin=-0,
            rescale=T)
text(-3,2,labels = "B)", font = 2, cex=1.2)



# Mantel test -------------------------------------------------------------





geo_nebrodi_adults<- read.table("Abies_nebrodi/openarray_data/diversity/coordenadas_nebrodi.txt", header=T)

source("/home/fbalao/Datos/R/script/RMantel2.R")

geoDist_nebrodi_adults<-geodist(geo_nebrodi_adults)
rownames(distgenDISS)<-c("01M","06M","08M","10M","11M","16M","18M","20M","21M",
                         "22M","27M","29M","02M","04M","07M","09M","12M","13M","14M",
                         "15M","17M","19M","23M","24M","25M","26M","28M","30M",
                         "31M","32M")
distgenDISS_2<-sort_dist_mat(distgenDISS, by_rows = T)

GeneD <- as.dist(distgenDISS_2)
GeoD<- as.dist(geoDist_nebrodi_adults)


manteltodo <- mantel(GeneD~GeoD, mrank=F, nperm=1000)
correlogram<- mgram(GeneD,GeoD, nperm=1000,mrank=T)
print(correlogram$mgram)
plot(correlogram,pval = 0.05, xlab = "Distance", ylab = "Mantel r", main=" Correlogram")

plot(mantel.rtest(GeoD,GeneD, nrepet = 9999), main = "Mantel's test")

# Structure ---------------------------------------------------------------
# Partiendo de nebrodf3 pasamos a otros formatos

source(here("scripts","diversidad","genind2structrure.R"))

genind2structure(nebrodf4, file="nebrodf4_structure.txt", pops=FALSE)

