library(vcfR)
library(hierfstat)
library(poppr)
library(ggplot2)
library(factoextra)
library(RColorBrewer)
#import vcf into geneind object

vcf_nebrodi <- read.vcfR("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_design/hybrids_filter/nebrodi_hybrids_snps_final.vcf", verbose = FALSE )

nebrodi_genind <- vcfR2genind(vcf_nebrodi)
pop<-c(rep("neb",5), rep("cep",15), rep("alb",15))
nebrodi_genind@pop<-as.factor(pop)

private.table<-private_alleles(nebrodi_genind, count.alleles=T, drop=T, report="data.frame")

private.table[private.table$population=="alb",3 ]<-private.table[private.table$population=="alb",3 ]/15
private.table[private.table$population=="cep",3 ]<-private.table[private.table$population=="cep",3 ]/15
private.table[private.table$population=="neb",3 ]<-private.table[private.table$population=="neb",3 ]/5

ggplot(private.table) + geom_tile(aes(x = population, y = allele, fill = count))

relative_frequencies <- pop.freq(nebrodi_genind)

library('dartR')
nebrodi_genl<-gi2gl(nebrodi_genind, parallel = TRUE)
gl.report.pa(nebrodi_genl)
gl.fixed.diff(nebrodi_genl, tloc=0.05, test=TRUE, delta=0.02, reps=100, v=5, pb=T)


pc <- gl.pcoa(nebrodi_genl, nfactors=5)
gl.pcoa.plot(pc, nebrodi_genl, xaxis=1, yaxis=2)


#Elegimos los alelos 10 con mayor carga de cada unos de los 2  PC principales en el PCA
# Ordenamos de manera decreciente lso alelos por carga en valor absoluto
pc1order<-order(abs(pc$loadings[,1]),decreasing = T)
pc1_10<-pc1order[1:20]
pc2order<-order(abs(pc$loadings[,2]),decreasing = T)
pc2_10<-pc2order[1:20]


pc1order5max<-order(pc$loadings[,1], decreasing = T)[1:5]
pc1order5min<-order(pc$loadings[,1], decreasing = F)[1:5]
pc2order5max<-order(pc$loadings[,2], decreasing = T)[1:5]
pc2order5min<-order(pc$loadings[,2], decreasing = F)[1:5]





selectednebrodi<-nebrodi_genl[,c(pc1_10,pc2_10)]
pc10 <- gl.pcoa(selectednebrodi, nfactors=2)
gl.pcoa.plot(pc10, selectednebrodi, pop.labels="pop", xaxis=1, yaxis=2, )

freqgenl<-function(x){apply(as.matrix(x),2, tapply, pop(x), function(e) mean(e)/ploidy(x)[1])}

sel_freq<-freqgenl(selectednebrodi)


selectednebrodi2<-nebrodi_genl[,c(pc1order5max,pc1order5min,pc2order5max,pc1order5min)]
pcmaxmin <- gl.pcoa(selectednebrodi2, nfactors=2)
gl.pcoa.plot(pcmaxmin, selectednebrodi2, pop.labels="pop", xaxis=1, yaxis=2)


#################################### Esto es que se hizo finalmente


dapc_neb<-dapc(nebrodi_genl,n.pca = 35, n.da = 2)
scatter(dapc_neb, 1, 2,posi.da="bottomright", bg="white")
assignplot(dapc_neb, cex=0.5)
compoplot(dapc_neb, lab="",ncol=1, xlab="individuals", col=rainbow(3))


contrib <- loadingplot(dapc_neb$var.contr, axis=1, threshold=quantile(dapc_neb$var.contr,0.99))
contrib2 <- loadingplot(dapc_neb$var.contr, axis=2, threshold=quantile(dapc_neb$var.contr,0.99))

dapc_sel1<-order(dapc_neb$var.contr[,1], decreasing = T)[1:11]
dapc_sel2<-order(dapc_neb$var.contr[,2], decreasing = T)[1:10]

dapc_neb$var.contr[order(dapc_neb$var.contr[,1], decreasing = T)[1:10],1]
dapc_neb$var.contr[order(dapc_neb$var.contr[,2], decreasing = T)[1:10],2]

# el loci #542 tiene alta carga en el eje 1 y 2, por lo tanto cogemos uno más del eje 1
# para que haya 20
selectednebrodi_dapc_old<-nebrodi_genl[,c(dapc_sel1[-1],dapc_sel2)]

# Actualizado 12-1-2021, se añaden 7 loci más de hibridación tras fallar 7 loci del primer lote
# de arrays por ello se cogen 15 del eje 1 y 13 del eje 2 (27= 28-1 común)

dapc_sel1<-order(dapc_neb$var.contr[,1], decreasing = T)[1:15]
dapc_sel2<-order(dapc_neb$var.contr[,2], decreasing = T)[1:13]
selectednebrodi_dapc<-nebrodi_genl[,c(dapc_sel1[-1],dapc_sel2)]

newloci<-unique(selectednebrodi_dapc$loc.names[! selectednebrodi_dapc$loc.names%in% selectednebrodi_dapc_old$loc.names])
pc_dapc <- gl.pcoa(selectednebrodi_dapc, nfactors=2)

col1=brewer.pal(5,"Spectral")[c(3,2,4)]



gl.pcoa.plot(pc_dapc, selectednebrodi_dapc, pop.labels="pop", xaxis=1, yaxis=2,
             pt.color=col1, pop.labels=c())

#### Figure S1

x<-pc_dapc$scores[,1]
y<-pc_dapc$scores[,2]


circ<-rep("grey45",35)

bgcol<-col1[pop(selectednebrodi_dapc)]



plot(x,y, bg=bgcol, pch=21, cex=2, col=circ,
     xlab="PC1 (42.5%)", ylab="PC2 (20.8%)", cex.axis=1.3, cex.lab=1.4)

legend(0.2,-1.2,legend = c(expression(italic("A. alba")), 
                        expression(italic("A. cephalonica")),
                        expression(italic("A. nebrodensis"))),
       pch=21, cex=1.2, col="grey35", pt.bg=col1, pt.cex=1.5,
       box.lwd = 0, border = "white", bg="transparent", box.col = "white")




selectednebrodi_dapc$loc.names

write.table(newloci, file="/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_design/hybrids_filter/nebrodi_hybrids_7snps_newdataset2021.txt", col.names = F, row.names = F, quote = F)


# Selection for paternity test ---------------------------------------------

#import vcf into geneind object

vcf_nebrodip <- read.vcfR("/home/fbalao/Datos/ARTICULO/Abies_nebrodi/openarray_design/nebrodi_filter/clean_intergenicSNPs_nebrodi.vcf", verbose = FALSE )

nebrodip_genind <- vcfR2genind(vcf_nebrodip)

pca1 <- dudi.pca(nebrodip_genind$tab,scannf=FALSE,scale=FALSE)# 

fviz_pca_ind(pca1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


fviz_pca_var(pca1,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,    # Avoid text overlapping
             select.var=list(contrib=50)
)cont<-get_pca_var(pca1)
cont$contrib

even_indexes<-seq(2,536,2)
c1e<-pca1$c1[even_indexes,]

top50pc1<-rownames(c1e[order(c1e[,1], decreasing = T)[1:73],])
top50pc2<-rownames(c1e[order(c1e[,2], decreasing = T)[1:68],])

selected<-unique(c(top50pc1,top50pc2))

new<-nebrodip_genind$tab[,selected]
pcanew <- dudi.pca(new,scannf=FALSE,scale=FALSE)# 

fviz_pca_ind(pcanew,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

snps100<-substr(colnames(new),1,nchar(colnames(new))-2)
write.table(snps100, file="nebrodi_105snps_newfinaldata.txt", col.names = F, row.names = F, quote = F)

