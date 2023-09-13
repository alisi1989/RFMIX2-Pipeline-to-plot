#leggiamo il file
msp=read.table("/Users/alessandrolisi/Desktop/MyData2/RFMIX2_pipeline_master/Output/Mozabite1_ind1_allchr.msp.tsv", header=F, sep="\t")

#prendiamo solo le prime 8 colonne
msp2=msp[,1:8]

#eliminiamo la colonna 4,5 e 6
msp2=msp2[,-4]
msp2=msp2[,-4]
msp2=msp2[,-4]

#rinominiamo le colonne
names(msp2)=c("chm", "spos", "epos", "ind1", "ind1")

#creiamo 2 dataframe
ch1=msp2[,-5]
ch2=msp2[,-4]



# Sostituisci i valori nella colonna "ind1" con i corrispondenti valori delle ancestry
ch1$ind1 <- ifelse(ch1$ind1 == 0, "ANC0",
                   ifelse(ch1$ind1 == 1, "ANC1",
                   ifelse(ch1$ind1 == 2, "ANC2",
                   ifelse(ch1$ind1 == 3, "ANC3",
                   ifelse(ch1$ind1 == 4, "ANC4", ch1$ind1)))))
                   
# Sostituisci i valori nella colonna "ind1" con i corrispondenti valori delle ancestry
ch2$ind1 <- ifelse(ch2$ind1 == 0, "ANC0",
                   ifelse(ch2$ind1 == 1, "ANC1",
                   ifelse(ch2$ind1 == 2, "ANC2",
                   ifelse(ch2$ind1 == 3, "ANC3",
                   ifelse(ch2$ind1 == 4, "ANC4", ch2$ind1)))))

#salviamo i due file per singolo individuo
write.table(ch2, "/Users/alessandrolisi/Desktop/MyData2/RFMIX2_pipeline_master/Output/Mozabite1_ind1_hap2.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(ch1, "/Users/alessandrolisi/Desktop/MyData2/RFMIX2_pipeline_master/Output/Mozabite1_ind1_hap1.bed", quote=F, col.names=F, row.names=F, sep="\t")