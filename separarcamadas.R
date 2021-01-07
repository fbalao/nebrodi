# Madres y progenie
madres_names<-names(which(table(nebrodf2[nebrodf2@other$Vivero_campo==0]@other$Madre_procedencia)>2))
non_corrected_madres_names<-names(which(table(nebrodf2[nebrodf2@other$Vivero_campo==0]@other$Madre_procedencia)<=2))


NA_mother_replacer<-function(genindobj){
M<-as.data.frame(genindobj@tab)
M_mother<- which(genindobj@other$Plantula_madre==0)
ch<-M[M_mother,]
nach<-which(is.na(ch))
ch[nach]<-apply(M[,nach],2, median, na.rm=T)
newM1.1<-ch
ch[which(ch==2)]<-1
ch[which(ch==0)]<-NA

M2 <- replace_na(M[-M_mother,], as.list(ch))
M3<- rbind(newM1.1,M2)
M4<-as.genind(M3)
M4@other<-genindobj@other
M4
}

Mcorrected<-list()
for (i in 1: length(madres_names)){
madres_names[i]
  GENO<-nebrodf2[nebrodf2@other$Madre_procedencia==madres_names[i] & nebrodf2@other$Vivero_campo==0]
  Mcorrected[[madres_names[i]]]<-NA_mother_replacer(GENO)
}

Nocorrected<-list()
for (i in 1: length(non_corrected_madres_names)){
  non_corrected_madres_names[i]
  Nocorrected[[non_corrected_madres_names[i]]]<-nebrodf2[nebrodf2@other$Madre_procedencia==non_corrected_madres_names[i] & nebrodf2@other$Vivero_campo==0]
}


  nebrodf3<-repool(c(Mcorrected,Nocorrected))
  
  db2<-codes[match(rownames(nebrodf3@tab), codes$Codigo_muestra),]
  
  nebrodf3@other$Plantula_madre<-as.factor(db2$Plantula_madre)
  nebrodf3@other$Vivero_campo<-as.factor(db2$Vivero_campo)
  nebrodf3@other$Madre_procedencia<-as.factor(db2$Madre_procedencia)


  
  # se eliminan 10 muestras más
  ty<-function(x){sum(is.na(x))}
  
  y<-matrix(nrow = 21, ncol=2)
  for (z in 1:21){
    i<-seq(0,1,0.05)[z]
    db<-missingno(nebrodf3, type = "loci", cutoff = i, quiet = F, freq = FALSE)
    y[z,1]<-nlevels(db@loc.fac)
    y[z,2]<-mean(apply(db@tab,1,ty)/(y[z,1]*2))
  }
  
  cuttoff<-cbind(seq(0,1,0.05),y)
  cuttNA<-cbind(seq(0,1,0.05),y[,2])
  plot( cuttoff,type="b", ylab="N loci", xlab="cut off call rate", ylim=c(0,120))
  plot( cuttNA,type="b", ylab="NA %", xlab="cut off call rate")
  abline(v=0.25, lty=2)
  
  
  
  nebrodf3<-missingno(nebrodf3, type = "geno", cutoff = 0.30, quiet = FALSE, freq = FALSE)
  nebrodf3<-missingno(nebrodf3, type = "loci", cutoff = 0.30, quiet = FALSE, freq = FALSE)
  nebrodf3<-informloci(nebrodf3, MAF=0)
