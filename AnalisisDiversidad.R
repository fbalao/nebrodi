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
##################################################################

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


distgenDIFF <- dist.gene(padresDBloci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
plot.phylo(bionj(distgenDIFF), type="p", cex=0.7, show.tip.label = F)
tiplabels(padresDB@other$Madre_procedencia)

distgenDISS <- diss.dist(padresDB, percent = T, mat = T)
plot(bionj(distgenDISS), type="p", cex=0.7, show.tip.label=F)
tiplabels(padresDB@other$Madre_procedencia)

padresnames<-read.table("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_data/diversity/muestras2.csv", header=T)
rownames(distgenDISS)

rownames(distgenDISS)<-padresnames[match(rownames(distgenDISS), padresnames$Codigo_muestra), "ID_LIFE"]
colnames(distgenDISS)<-padresnames[match(rownames(distgenDISS), padresnames$Codigo_muestra), "ID_LIFE"]

ps<-propShared(padresDB)
rownames(ps)<-padresnames[match(rownames(ps), padresnames$Codigo_muestra), "ID_LIFE"]
colnames(ps)<-padresnames[match(rownames(ps), padresnames$Codigo_muestra), "ID_LIFE"]


corrplot(1-distgenDISS,  order = "hclust", hclust.method="ward.D2",
         addrect = 3,
         is.corr = F,
         cl.lim = c(0,1),
         col=col2(100),
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


##### clusters

kneb<-find.clusters(padresDB)
neb_dapc<-dapc(padresDB,pop = kneb$grp)
rownames(neb_dapc$tab)<-padresnames[match(rownames(neb_dapc$tab), padresnames$Codigo_muestra), "ID_LIFE"]
names(neb_dapc$grp)<-padresnames[match(names(neb_dapc$grp), padresnames$Codigo_muestra), "ID_LIFE"]

v<-names(sort(neb_dapc$grp))

layout(matrix(c(1,2),nrow=2))
compoplot(neb_dapc,border=T,
          col=funky(6),show.lab=T, legend=F,
          subset=v)

temp <- optim.a.score(neb_dapc)

neb_dapc<-dapc(padresDB,pop = kneb$grp)


rownames(padresDB@tab)<-padresnames[match(rownames(padresDB@tab), padresnames$Codigo_muestra), "ID_LIFE"]


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
            



pairDistPlot(padresDB, neb_dapc$grp)

###################3 Diversity

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

############3333333

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


corrplot(x.dist,  order = "hclust", hclust.method = "ward.D2",
         addrect = 3,
         is.corr = F,
      # cl.lim = c(-1,1),
         col=col3(100),
         win.asp = 1,
         cl.ratio=0.2,
        tl.col = funky(6)[c(rep(2,12),rep(3,10),rep(1,8))])


  dfr <- reshape(rela[c(2,3,9)], direction="wide", idvar="ind2.id", timevar="ind1.id")

dfr<-xtabs(ritland~.,rela[c(2,3,9)])
d <- as.dist(dfr)
attr(d, "Labels") <- c("1M",dfr[-29, 1])
##############################333 Tabla diversidad
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
