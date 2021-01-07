library(tidyr)
library(dplyr)

df<-data.frame(x=c(0,1,1,1,0), y=c(NA,NA,1,1,0), z=c(NA,0,0,0,NA))
df<-as.data.frame(t(df))
df

#Elige las columas impares de la madre
ch<-df[1,1:(dim(df)[2])%% 2 != 0]

#Cambia los NA  de las columnas impares por los valores de la madre
df %>% replace_na(as.list(ch))
