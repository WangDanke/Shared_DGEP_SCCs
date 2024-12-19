##preprocess the metadata of ENA RNA-seq datasets
setwd("D:\\SCC_DGEP_shared\\SCC_ENA_seq")
getwd()
library(tidyr)
###
phe_GSE87410 = read.csv("GSE87410.csv", header = F)
phe_GSE87410 = phe_GSE87410[-c(1:29),]

phe_GSE87410 = phe_GSE87410 %>% t()
colnames(phe_GSE87410) = phe_GSE87410[1,]
phe_GSE87410 = phe_GSE87410[-1,]
rownames(phe_GSE87410) = phe_GSE87410[,2]
phe_GSE87410 = as.data.frame(phe_GSE87410)

pheno_GSE87410 = phe_GSE87410[,c("!Sample_title","ID_REF","!Sample_source_name_ch1",
                                 "!Sample_organism_ch1" ,"!Sample_characteristics_ch1","!Sample_library_strategy")]

library(dplyr)
library(tidyr)
library(stringr)
sample_title = str_split_fixed(pheno_GSE87410$`!Sample_title`,"_",3)
sample_title = as.data.frame(sample_title)
str(sample_title)
pheno_GSE87410 = cbind(pheno_GSE87410,sample_title)
str(pheno_GSE87410)
colnames(pheno_GSE87410) = c("Sample_title", "ID", "source", "sample_organism","Sample_type","organ")
save(pheno_GSE87410,file = "pheno_GSE87410.RData")

####
####cervical SCC

##GSE223804
phe_GSE223804 = read.csv("GSE223804_clinical.csv", header = F)
phe_GSE223804 = phe_GSE223804[-c(1:34),]
phe_GSE223804 = phe_GSE223804 %>% t()
colnames(phe_GSE223804) = phe_GSE223804[1,]
phe_GSE223804 = phe_GSE223804[-1,]
rownames(phe_GSE223804) = phe_GSE223804[,2]
phe_GSE223804 = as.data.frame(phe_GSE223804)

pheno_GSE223804 = phe_GSE223804[,c(1,46,10,11)]号
AA = as.data.frame(str_split_fixed(pheno_GSE223804$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE223804$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE223804 = cbind(pheno_GSE223804[,c(1,2)],AA[,2],BB[,2])
save(pheno_GSE223804,file = "pheno_GSE223804.RData")

###GSE144293
phe_GSE144293 = read.csv("GSE144293_clinical.csv", header = F)
phe_GSE144293 = phe_GSE144293[-c(1:32),]
phe_GSE144293 = phe_GSE144293 %>% t()
colnames(phe_GSE144293) = phe_GSE144293[1,]
phe_GSE144293 = phe_GSE144293[-1,]
rownames(phe_GSE144293) = phe_GSE144293[,2]
phe_GSE144293 = as.data.frame(phe_GSE144293)

pheno_GSE144293 = phe_GSE144293[,c(1,46,10,12)]
AA = as.data.frame(str_split_fixed(pheno_GSE144293$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE144293$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE144293 = cbind(pheno_GSE144293[,c(1,2)],AA[,2],BB[,2])
colnames(pheno_GSE144293) = c("Sample_title","ID","organ","tumor_type")
pheno_GSE144293$sample_type = "T"
save(pheno_GSE144293,file = "pheno_GSE144293.RData")

###
###cSCC

###GSE84293
phe_GSE84293 = read.csv("GSE84293_clinical.csv", header = F)
phe_GSE84293 = phe_GSE84293[-c(1:40),]
phe_GSE84293 = phe_GSE84293 %>% t()
colnames(phe_GSE84293) = phe_GSE84293[1,]
phe_GSE84293 = phe_GSE84293[-1,]
rownames(phe_GSE84293) = phe_GSE84293[,2]
phe_GSE84293 = as.data.frame(phe_GSE84293)

colnames(phe_GSE84293)
pheno_GSE84293 = phe_GSE84293[,c(1,2,8,10)]
AA = as.data.frame(str_split_fixed(pheno_GSE84293$`!Sample_characteristics_ch1`,":",2))

pheno_GSE84293 = cbind(pheno_GSE84293[,c(1,2,3)],AA[,2])
colnames(pheno_GSE84293) = c("Sample_title","ID","organ","tumor_type")

save(pheno_GSE84293,file = "pheno_GSE84293.RData")

####GSE139505
phe_GSE139505 = read.csv("GSE139505_clinical.csv", header = F)
phe_GSE139505 = phe_GSE139505[-c(1:33),]
phe_GSE139505 = phe_GSE139505 %>% t()
colnames(phe_GSE139505) = phe_GSE139505[1,]
phe_GSE139505 = phe_GSE139505[-1,]
rownames(phe_GSE139505) = phe_GSE139505[,2]
phe_GSE139505 = as.data.frame(phe_GSE139505)

colnames(phe_GSE139505)
pheno_GSE139505 = phe_GSE139505[,c(1,2,11,12,13)]
pheno_GSE139505$`!Sample_characteristics_ch1` = "skin"
AA = as.data.frame(str_split_fixed(pheno_GSE139505$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE139505$`!Sample_characteristics_ch1.2`,":",2))
pheno_GSE139505 = cbind(pheno_GSE139505[,c(1,2,3)],AA[,2],BB[,2])
colnames(pheno_GSE139505)= c("Sample_title","ID","organ","gender","Sample_type")
pheno_GSE139505$Sample_type = c("N","N","N","N","N","N","N","T","T","T","T","T","T","T","T","T")

save(pheno_GSE139505,file = "pheno_GSE139505.RData")

####GSE191334
phe_GSE191334 = read.csv("GSE191334_clinical.csv", header = F)
phe_GSE191334 = phe_GSE191334[-c(1:29),]
phe_GSE191334 = phe_GSE191334 %>% t()
colnames(phe_GSE191334) = phe_GSE191334[1,]
phe_GSE191334 = phe_GSE191334[-1,]
rownames(phe_GSE191334) = phe_GSE191334[,2]
phe_GSE191334 = as.data.frame(phe_GSE191334)

colnames(phe_GSE191334)
pheno_GSE191334 = phe_GSE191334[,c(1,2,8)]
View(pheno_GSE191334)
colnames(pheno_GSE191334) = c("Sample_title","ID","organ")
pheno_GSE191334$sample_type = c("N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T")

save(pheno_GSE191334,file = "pheno_GSE191334.RData")

###
###ESCC

###GSE149612
phe_GSE149612 = read.csv("GSE149612_clinical.csv", header = F)
phe_GSE149612 = phe_GSE149612[-c(1:32),]
phe_GSE149612 = phe_GSE149612 %>% t()
colnames(phe_GSE149612) = phe_GSE149612[1,]
phe_GSE149612 = phe_GSE149612[-1,]
rownames(phe_GSE149612) = phe_GSE149612[,2]
phe_GSE149612 = as.data.frame(phe_GSE149612)

colnames(phe_GSE149612)
pheno_GSE149612 = phe_GSE149612[,c(1, 2, 8, 12,18)]

AA = as.data.frame(str_split_fixed(pheno_GSE149612$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE149612$`!Sample_characteristics_ch1.2`,":",2))
colnames(pheno_GSE149612) = c("Sample_title","ID", "Sample","stage","Sample_type")
pheno_GSE149612$organ = "esophagus"
save(pheno_GSE149612,file = "pheno_GSE149612.RData")


##GSE164158
phe_GSE164158 = read.csv("GSE164158_clinical.csv", header = F)
phe_GSE164158 = phe_GSE164158[-c(1:34),]
phe_GSE164158 = phe_GSE164158 %>% t()
colnames(phe_GSE164158) = phe_GSE164158[1,]
phe_GSE164158 = phe_GSE164158[-1,]
rownames(phe_GSE164158) = phe_GSE164158[,2]
phe_GSE164158 = as.data.frame(phe_GSE164158)

colnames(phe_GSE164158)
pheno_GSE164158 = phe_GSE164158[,c(1,2,8,10,13,14,15,16)]
AA = as.data.frame(str_split_fixed(pheno_GSE164158$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE164158$`!Sample_characteristics_ch1.2`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE164158$`!Sample_characteristics_ch1.3`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE164158$`!Sample_characteristics_ch1.4`,":",2))
pheno_GSE164158 = cbind(pheno_GSE164158[,c(1:3)],AA[,2],BB[,2],CC[,2],DD[,2])
colnames(pheno_GSE164158) = c("Sample_title","ID","Sample_type","age","gender","t_stage","n_stage")
pheno_GSE164158$Sample_type = c("N","N","N","N","N","N","N","N","T","T","T","T","T","T","T","T")
save(pheno_GSE164158,file = "pheno_GSE164158.RData")



###GSE167488
phe_GSE167488 = read.csv("GSE167488_clinical.csv", header = F)
phe_GSE167488 = phe_GSE167488[-c(1:27),]
phe_GSE167488 = phe_GSE167488 %>% t()
colnames(phe_GSE167488) = phe_GSE167488[1,]
phe_GSE167488 = phe_GSE167488[-1,]
rownames(phe_GSE167488) = phe_GSE167488[,2]
phe_GSE167488 = as.data.frame(phe_GSE167488)

colnames(phe_GSE167488)
pheno_GSE167488 = phe_GSE167488[,c(1,2,8,10,13,14,15,16)]
AA = as.data.frame(str_split_fixed(pheno_GSE167488$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE167488$`!Sample_characteristics_ch1.2`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE167488$`!Sample_characteristics_ch1.3`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE167488$`!Sample_characteristics_ch1.4`,":",2))
pheno_GSE167488 = cbind(pheno_GSE167488[,c(1:4)],AA[,2],BB[,2],CC[,2],DD[,2])
colnames(pheno_GSE167488) = c("Sample_title","ID","sample","sample_type","age","stage","smoke","grade")
pheno_GSE167488$sample_type = c("N","N","N","T","T","T","N","N","N","T","T","T")
pheno_GSE167488$smoke = "1"
save(pheno_GSE167488,file = "pheno_GSE167488.RData")


###GSE189830
phe_GSE189830 = read.csv("GSE189830_clinical.csv", header = F)
phe_GSE189830 = phe_GSE189830[-c(1:36),]
phe_GSE189830 = phe_GSE189830 %>% t()
colnames(phe_GSE189830) = phe_GSE189830[1,]
phe_GSE189830 = phe_GSE189830[-1,]
rownames(phe_GSE189830) = phe_GSE189830[,2]
phe_GSE189830 = as.data.frame(phe_GSE189830)

colnames(phe_GSE189830)
pheno_GSE189830 = phe_GSE189830[,c(1,2,11)]
AA = as.data.frame(str_split_fixed(pheno_GSE189830$`!Sample_characteristics_ch1`,":",2))

pheno_GSE189830 = cbind(pheno_GSE189830[,c(1:2)],AA[,2])
pheno_GSE189830$organ = "esophagus"
colnames(pheno_GSE189830) = c("Sample_title","ID","sample","organ")
pheno_GSE189830$sample_type = c("N","N","T","N","N","T","N","N","T","N","N","T")
pheno_GSE189830$stage = c("II","II","II","III","III","III","0","0","0","II","II","II")
save(pheno_GSE189830,file = "pheno_GSE189830.RData")


###GSE194116
phe_GSE194116 = read.csv("GSE194116_clinical.csv", header = F)
phe_GSE194116 = phe_GSE194116[-c(1:31),]
phe_GSE194116 = phe_GSE194116 %>% t()
colnames(phe_GSE194116) = phe_GSE194116[1,]
phe_GSE194116 = phe_GSE194116[-1,]
rownames(phe_GSE194116) = phe_GSE194116[,2]
phe_GSE194116 = as.data.frame(phe_GSE194116)

colnames(phe_GSE194116)
pheno_GSE194116 = phe_GSE194116[,c(1,2,8,12)]

colnames(pheno_GSE194116) = c("Sample_title","ID","organ","sample_type")
pheno_GSE194116$sample_type = c("T","N","T","N","T","N","T","N","T","N","T","N")

save(pheno_GSE194116,file = "pheno_GSE194116.RData")


####HNSC

###GSE182227
phe_GSE182227 = read.csv("GSE182227_clinical.csv", header = F)
phe_GSE182227 = phe_GSE182227[-c(1:28),]
phe_GSE182227 = phe_GSE182227 %>% t()
colnames(phe_GSE182227) = phe_GSE182227[1,]
phe_GSE182227 = phe_GSE182227[-1,]
rownames(phe_GSE182227) = phe_GSE182227[,2]
phe_GSE182227 = as.data.frame(phe_GSE182227)

colnames(phe_GSE182227)
pheno_GSE182227 = phe_GSE182227[,c(1,2,8,10,11)]
AA = as.data.frame(str_split_fixed(pheno_GSE182227$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE182227$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE182227 = cbind(pheno_GSE182227[,c(1:3)],AA[,2],BB[,2])
colnames(pheno_GSE182227) = c("Sample_title","ID","organ","sample_type-1","HPV_state")
pheno_GSE182227$sample_type_2 = ifelse(pheno_GSE182227$`sample_type-1`==" tumor","T",pheno_GSE182227$`sample_type-1`)
pheno_GSE182227$sample_type = ifelse(pheno_GSE182227$`sample_type_2`==" adjacent normal tissue","N",pheno_GSE182227$`sample_type_2`)
pheno_GSE182227 = pheno_GSE182227[,c(1,2,3,5,6)]
save(pheno_GSE182227,file = "pheno_GSE182227.RData")


###GSE205308
phe_GSE205308 = read.csv("GSE205308_clinical.csv", header = F)
phe_GSE205308 = phe_GSE205308[-c(1:28),]
phe_GSE205308 = phe_GSE205308 %>% t()
colnames(phe_GSE205308) = phe_GSE205308[1,]
phe_GSE205308 = phe_GSE205308[-1,]
rownames(phe_GSE205308) = phe_GSE205308[,2]
phe_GSE205308 = as.data.frame(phe_GSE205308)

colnames(phe_GSE205308)
pheno_GSE205308 = phe_GSE205308[,c(1,2,8,10,11)]
AA = as.data.frame(str_split_fixed(pheno_GSE205308$`!Sample_characteristics_ch1`,":",2))

pheno_GSE205308 = cbind(pheno_GSE205308[,c(1:3)],AA[,2])
colnames(pheno_GSE205308) = c("Sample_title","ID","organ","sample_type-1","HPV_state")
pheno_GSE205308$HPV_state = ifelse(pheno_GSE205308$`AA[, 2]`==" HPV_associated","HPV+","HPV-")
pheno_GSE205308$sample_type = "T"
pheno_GSE205308 = pheno_GSE205308[,c(1,2,3,5,6)]
colnames(pheno_GSE205308) = c("Sample_title","ID","organ","PHV_state","sample_type")
save(pheno_GSE205308,file = "pheno_GSE205308.RData")


###GSE208576
phe_GSE208576 = read.csv("GSE208576_clinical.csv", header = F)
phe_GSE208576 = phe_GSE208576[-c(1:35),]
phe_GSE208576 = phe_GSE208576 %>% t()
colnames(phe_GSE208576) = phe_GSE208576[1,]
phe_GSE208576 = phe_GSE208576[-1,]
rownames(phe_GSE208576) = phe_GSE208576[,2]
phe_GSE208576 = as.data.frame(phe_GSE208576)

colnames(phe_GSE208576)
pheno_GSE208576 = phe_GSE208576[,c(1,2,8,12,13,14,15,16,17,18,19,20,21,22,23)]

AA = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.1`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.2`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.3`,":",2))
EE = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.4`,":",2))
FF = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.5`,":",2))
GG = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.6`,":",2))
HH = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.7`,":",2))
II = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.8`,":",2))
JJ = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.9`,":",2))
KK = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.10`,":",2))
LL = as.data.frame(str_split_fixed(pheno_GSE208576$`!Sample_characteristics_ch1.11`,":",2))

pheno_GSE208576 = cbind(pheno_GSE208576[,c(1:3)],AA[,2],BB[,2],CC[,2],DD[,2],
                        EE[,2],FF[,2],GG[,2],HH[,2],II[,2],JJ[,2],KK[,2],LL[,2])
colnames(pheno_GSE208576) = c("Sample_title","ID","Sampele_type","age","gender","smoke","drink","organ",
                              "t_stage","grade","HPV_state_1","radiotherapy","chemotherapy","live_state","live.time_momth")
pheno_GSE208576$Sampele_type = "T"
pheno_GSE208576$HPV_state = ifelse(pheno_GSE208576$HPV_state_1==" Positive","HPV+(16)","HPV-(16)")
pheno_GSE208576 = pheno_GSE208576[,-c(11)]

save(pheno_GSE208576,file = "pheno_GSE208576.RData")


###GSE228393
phe_GSE228393 = read.csv("GSE228393_clinical.csv", header = F)
phe_GSE228393 = phe_GSE228393[-c(1:25),]
phe_GSE228393 = phe_GSE228393 %>% t()
colnames(phe_GSE228393) = phe_GSE228393[1,]
phe_GSE228393 = phe_GSE228393[-1,]
rownames(phe_GSE228393) = phe_GSE228393[,2]
phe_GSE228393 = as.data.frame(phe_GSE228393)

colnames(phe_GSE228393)
pheno_GSE228393 = phe_GSE228393[,c(1,2,8)]

colnames(pheno_GSE228393) = c("Sample_title","ID","organ")
pheno_GSE228393$sample_type = c("T","T","T","T","T","N","N","N")

save(pheno_GSE228393,file = "pheno_GSE228393.RData")



####LUSC 

###GSE226069
phe_GSE226069 = read.csv("GSE226069_clinical.csv", header = F)
phe_GSE226069 = phe_GSE226069[-c(1:29),]
phe_GSE226069 = phe_GSE226069 %>% t()
colnames(phe_GSE226069) = phe_GSE226069[1,]
phe_GSE226069 = phe_GSE226069[-1,]
rownames(phe_GSE226069) = phe_GSE226069[,2]
phe_GSE226069 = as.data.frame(phe_GSE226069)

colnames(phe_GSE226069)
pheno_GSE226069 = phe_GSE226069[,c(1,2,8,11,12)]
AA = as.data.frame(str_split_fixed(pheno_GSE226069$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE226069$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE226069 = cbind(pheno_GSE226069[,c(1:3)],AA[,2],BB[,2])
colnames(pheno_GSE226069) = c("Sample_type","ID","organ","tumor_type","trans_state")
pheno_GSE226069$sample_type = "T"
save(pheno_GSE226069,file = "pheno_GSE226069.RData")


###GSE159857
phe_GSE159857 = read.csv("GSE159857_clinical.csv", header = F)
phe_GSE159857 = phe_GSE159857[-c(1:31),]
phe_GSE159857 = phe_GSE159857 %>% t()
colnames(phe_GSE159857) = phe_GSE159857[1,]
phe_GSE159857 = phe_GSE159857[-1,]
rownames(phe_GSE159857) = phe_GSE159857[,2]
phe_GSE159857 = as.data.frame(phe_GSE159857)

colnames(phe_GSE159857)
pheno_GSE159857 = phe_GSE159857[,c(1,2,8,11,12)]
AA = as.data.frame(str_split_fixed(pheno_GSE159857$`!Sample_characteristics_ch1`,":",2))

pheno_GSE159857 = cbind(pheno_GSE159857[,c(1,2,4,5)],AA[,2])
pheno_GSE159857$sample_type = ifelse(pheno_GSE159857$`!Sample_characteristics_ch1.2`=="disease state: Tumor","T","N")
pheno_GSE159857 = pheno_GSE159857[,c(1,2,3,5,6)]
colnames(pheno_GSE159857) = c("Sample_type","ID","patient","tumor_type","sample_type")
pheno_GSE159857$organ = "Lung"
save(pheno_GSE159857,file = "pheno_GSE159857.RData")


###GSE179879
phe_GSE179879 = read.csv("GSE179879_clinical.csv", header = F)
phe_GSE179879 = phe_GSE179879[-c(1:43),]
phe_GSE179879 = phe_GSE179879 %>% t()
colnames(phe_GSE179879) = phe_GSE179879[1,]
phe_GSE179879 = phe_GSE179879[-1,]
rownames(phe_GSE179879) = phe_GSE179879[,2]
phe_GSE179879 = as.data.frame(phe_GSE179879)

colnames(phe_GSE179879)
pheno_GSE179879 = phe_GSE179879[,c(1,2,8,11,12)]
AA = as.data.frame(str_split_fixed(pheno_GSE179879$`!Sample_characteristics_ch1`,":",2))

pheno_GSE179879 = cbind(pheno_GSE179879[,c(1,2,3,5)],AA[,2])
pheno_GSE179879$sample_type = ifelse(pheno_GSE179879$`!Sample_characteristics_ch1.2`=="disease state: Tumor","T","N")
colnames(pheno_GSE179879) = c("Sample_title","ID","source","Sample_type_1","tumor_type")
pheno_GSE179879$sample_type = "T"
pheno_GSE179879$organ = "Lung"
save(pheno_GSE179879,file = "pheno_GSE179879.RData")

###GSE226069
phe_GSE226069 = read.csv("GSE226069_clinical.csv", header = F)
phe_GSE226069 = phe_GSE226069[-c(1:29),]
phe_GSE226069 = phe_GSE226069 %>% t()
colnames(phe_GSE226069) = phe_GSE226069[1,]
phe_GSE226069 = phe_GSE226069[-1,]
rownames(phe_GSE226069) = phe_GSE226069[,2]
phe_GSE226069 = as.data.frame(phe_GSE226069)

colnames(phe_GSE226069)
pheno_GSE226069 = phe_GSE226069[,c(1,2,8,11,12)]
pheno_GSE226069$tumor_type = ifelse(pheno_GSE226069$`!Sample_characteristics_ch1`=="disease state: Lung adenocarcinoma",
                                     "Adenocarcinoma","Squamous cell carcinoma")
pheno_GSE226069$trans_state = ifelse(pheno_GSE226069$`!Sample_characteristics_ch1.1`=="genotype: Transforming",
                                     "Transforming","Control")
pheno_GSE226069 = pheno_GSE226069[,c(1,2,3,7,8)]
colnames(pheno_GSE226069) = c("Sample_title","ID","organ","tumor_type","trans_state")
pheno_GSE226069$sample_type = "T"

save(pheno_GSE226069,file = "pheno_GSE226069.RData")


###GSE226070
phe_GSE226070 = read.csv("GSE226070-GPL20301_clinical.csv", header = F)
phe_GSE226070 = phe_GSE226070[-c(1:29),]
phe_GSE226070 = phe_GSE226070 %>% t()
colnames(phe_GSE226070) = phe_GSE226070[1,]
phe_GSE226070 = phe_GSE226070[-1,]
rownames(phe_GSE226070) = phe_GSE226070[,2]
phe_GSE226070 = as.data.frame(phe_GSE226070)

colnames(phe_GSE226070)
pheno_GSE226070 = phe_GSE226070[,c(1,2,8,11,12)]
pheno_GSE226070$tumor_type = ifelse(pheno_GSE226070$`!Sample_characteristics_ch1`=="disease state: Lung adenocarcinoma",
                                    "Adenocarcinoma","Squamous cell carcinoma")
pheno_GSE226070$trans_state = ifelse(pheno_GSE226070$`!Sample_characteristics_ch1.1`=="genotype: Transforming",
                                     "Transforming","Control")
pheno_GSE226070 = pheno_GSE226070[,c(1,2,3,6,7)]
colnames(pheno_GSE226070) = c("Sample_title","ID","organ","tumor_type","trans_state")
pheno_GSE226070$sample_type = "T"

save(pheno_GSE226070,file = "pheno_GSE226070.RData")


###GSE84339
phe_GSE84339 = read.csv("GSE84339_clinical.csv", header = F)
phe_GSE84339 = phe_GSE84339[-c(1:33),]
phe_GSE84339 = phe_GSE84339 %>% t()
colnames(phe_GSE84339) = phe_GSE84339[1,]
phe_GSE84339 = phe_GSE84339[-1,]
rownames(phe_GSE84339) = phe_GSE84339[,2]
phe_GSE84339 = as.data.frame(phe_GSE84339)

colnames(phe_GSE84339)
pheno_GSE84339 = phe_GSE84339[,c(1,2,8,11,12)]
AA = as.data.frame(str_split_fixed(pheno_GSE84339$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE84339$`!Sample_characteristics_ch1.1`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE84339$`!Sample_characteristics_ch1.2`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE84339$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE84339 = cbind(pheno_GSE84339[,c(1,2,3)],AA[,2],BB[,2],CC[,2],DD[,2])
colnames(pheno_GSE84339) = c("Sample_title","ID","organ","age","gender","pathologial_stage","smoke_states")
pheno_GSE84339$organ = "Lung"
pheno_GSE84339$tumor_type = "Squamous cell carcinoma"
pheno_GSE84339$sample_type = "T"
save(pheno_GSE84339,file = "pheno_GSE84339.RData")



####HNSC 

###GSE20116
phe_GSE20116 = read.csv("GSE20116_clinical.csv", header = F)
phe_GSE20116 = phe_GSE20116[-c(1:27),]
phe_GSE20116 = phe_GSE20116 %>% t()
colnames(phe_GSE20116) = phe_GSE20116[1,]
phe_GSE20116 = phe_GSE20116[-1,]
rownames(phe_GSE20116) = phe_GSE20116[,2]
phe_GSE20116 = as.data.frame(phe_GSE20116)

colnames(phe_GSE20116)
pheno_GSE20116 = phe_GSE20116[,c(1,2,8)]
colnames(pheno_GSE20116) = c("Sample_title","ID","organ")
pheno_GSE20116$sample_type = c("N","T","N","T","N","T")
pheno_GSE20116$organ = "Oral"
save(pheno_GSE20116,file = "pheno_GSE20116.RData")



###GSE176221
phe_GSE176221 = read.csv("GSE176221_clinical.csv", header = F)
phe_GSE176221 = phe_GSE176221[-c(1:27),]
phe_GSE176221 = phe_GSE176221 %>% t()
colnames(phe_GSE176221) = phe_GSE176221[1,]
phe_GSE176221 = phe_GSE176221[-1,]
rownames(phe_GSE176221) = phe_GSE176221[,2]
phe_GSE176221 = as.data.frame(phe_GSE176221)

colnames(phe_GSE176221)
pheno_GSE176221 = phe_GSE176221[,c(1,2,11,12,13,14,15)]
AA = as.data.frame(str_split_fixed(pheno_GSE176221$`!Sample_characteristics_ch1.4`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE176221$`!Sample_characteristics_ch1.1`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE176221$`!Sample_characteristics_ch1.2`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE176221$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE176221 = cbind(pheno_GSE176221[,c(1,2,3)],AA[,2],BB[,2],CC[,2],DD[,2])
colnames(pheno_GSE176221) = c("Sample_title","ID","source","grade","gender","age","t_states")
pheno_GSE176221$sample_type = "T"
pheno_GSE176221$organ = "Oral"
save(pheno_GSE176221,file = "pheno_GSE176221.RData")


###GSE184616
phe_GSE184616 = read.csv("GSE184616_clinical.csv", header = F)
phe_GSE184616 = phe_GSE184616[-c(1:31),]
phe_GSE184616 = phe_GSE184616 %>% t()
colnames(phe_GSE184616) = phe_GSE184616[1,]
phe_GSE184616 = phe_GSE184616[-1,]
rownames(phe_GSE184616) = phe_GSE184616[,2]
phe_GSE184616 = as.data.frame(phe_GSE184616)

colnames(phe_GSE184616)
pheno_GSE184616 = phe_GSE184616[,c(1,2,8,10,13,14,15,16)]
pheno_GSE184616$HPV_state = ifelse(pheno_GSE184616$`!Sample_characteristics_ch1`=="patient diagnosis: HPV-negative Oral Squamous Cell Carcinoma","HPV-","NA")
pheno_GSE184616 = pheno_GSE184616[,c(1:3,5:9)]
AA = as.data.frame(str_split_fixed(pheno_GSE184616$`!Sample_characteristics_ch1.4`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE184616$`!Sample_characteristics_ch1.1`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE184616$`!Sample_characteristics_ch1.2`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE184616$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE184616 = cbind(pheno_GSE184616[,c(1,2,3,4)],AA[,2],BB[,2],CC[,2],DD[,2])

pheno_GSE184616$sample_type = ifelse(pheno_GSE184616$`BB[, 2]`==" Adjacent Normal","N","T")
pheno_GSE184616 = pheno_GSE184616[,c(1:4,6:9)]
colnames(pheno_GSE184616) = c("Sample_title","ID","source","smoke_state","age","gender","HPV_state","sample_type")
save(pheno_GSE184616,file = "pheno_GSE184616.RData")

###GSE186775
phe_GSE186775 = read.csv("GSE186775_clinical.csv", header = F)
phe_GSE186775 = phe_GSE186775[-c(1:27),]
phe_GSE186775 = phe_GSE186775 %>% t()
colnames(phe_GSE186775) = phe_GSE186775[1,]
phe_GSE186775 = phe_GSE186775[-1,]
rownames(phe_GSE186775) = phe_GSE186775[,2]
phe_GSE186775 = as.data.frame(phe_GSE186775)

colnames(phe_GSE186775)
pheno_GSE186775 = phe_GSE186775[,c(1,2,10)]
View(pheno_GSE186775)
pheno_GSE186775$sample_type = ifelse(pheno_GSE186775$`!Sample_characteristics_ch1` == "	tissue: OSCC tissue","T","N")
pheno_GSE186775$sample_type = ifelse(pheno_GSE186775$`!Sample_characteristics_ch1` == "tissue: OSCC tissue","T","N")
pheno_GSE186775 = pheno_GSE186775[,c(1,2,4)]
colnames(pheno_GSE186775) = c("Sample_title","ID","sample_type")
pheno_GSE186775$organ = "Oral"
save(pheno_GSE186775,file = "pheno_GSE186775.RData")

#
pheno_HNSC = rbind(pheno_GSE176221[,c("Sample_title","ID","sample_type","organ","ID_SRR")],
pheno_GSE184616[,c("Sample_title","ID","sample_type","organ","ID_SRR")],
pheno_GSE186775[,c("Sample_title","ID","sample_type","organ","ID_SRR")])



































































































