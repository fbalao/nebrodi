#Format genind to colony2


db3<-codes[match(rownames(nebrodf3@tab), codes$Codigo_muestra),]
campo<-db3[db3$Vivero_campo==0,]

M10<-campo[campo$Madre_procedencia==10, ]

known_mother<-list()
for (i in 1: length(madres_names)){
  madres_names[i]
  genotypes<-campo[campo$Madre_procedencia==madres_names[i], ]
  known_mother[[i]]<-genotypes[order(genotypes$Plantula_madre),1]
}

lapply(known_mother, write, here("data","formats","nebrodi_KnownMatSibs.txt"), append=TRUE, ncolumns=1000, sep="\t")

markers<-rbind(as.character(nebrodf3@loc.fac), rep("0", length(nebrodf3@loc.fac)),rep("0", length(nebrodf3@loc.fac)), rep("0.01", length(nebrodf3@loc.fac) ))
markers<-markers[,1:(dim(markers)[2])%% 2 != 0]
write.table(markers,file=here("data","formats","nebrodi_markers.txt"), sep="\t", col.names = F, row.names = F, quote = F)
