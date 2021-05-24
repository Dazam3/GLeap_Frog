#This Script serves to organize genoic coordinates within a single genome of reference and
#evaluate conservation of corridinates among data sets to allow for visualization within 
#PCA plots

##################################################################################
#Initiate required packages/dependencies
##################################################################################

setwd('/Volumes/Cocoon/WGA/Z_ATAC/BED_Files/') #Set working directory
set.seed(1234)
library(ggplot2)
library(ggfortify)

##################################################################################
#Read in data table generated using GLeap-Frog
##################################################################################

Input <- read.table('./RPKM_Master.txt', header = T)
Input[1] <- NULL
NfAl <- cbind(Input[,1:11], Input[,12:14])
NfOl <- cbind(Input[,1:11],Input[,15:24])
NfDr <- cbind(Input[,1:11],Input[,25:34])
Fall <- Input

Fall <- na.omit(Fall)
NfAl <- na.omit(NfAl)
NfOl <- na.omit(NfOl)
NfDr <- na.omit(NfDr)
# Optionally used to convert NA to zeros if user beleives all columns should be included in analysis
#Input2[is.na(Input2)] <- 0

Fall <- t(Fall)
NfAl <- t(NfAl)
NfOl <- t(NfOl)
NfDr <- t(NfDr)
counter <- c(ncol(Fall),ncol(NfAl),ncol(NfOl),ncol(NfDr))

##################################################################################
#Transform data to matrixies and eigan vectors
##################################################################################
test1 <- prcomp(Fall, scale.=T)
test2 <- prcomp(NfAl, scale.=T)
test3 <- prcomp(NfOl, scale.=T)
test4 <- prcomp(NfDr, scale.=T)
group11 <- c('Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur', 'Alim', 'Alim','Alim', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat','Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer')
group12 <- c('Nf_Dia6', 'Nf_Dia6', 'Nf_esc', 'Nf_esc', 'Nf_KvO', 'Nf_KvY', 'Nf_KvY', 'Nf_DiaM', 'Nf_DiaM', 'Nf_DiaM', 'Nf_KvO', 'Al_DiaM', 'Al_esc', 'Al_Kv', 'Ol_11', 'Ol_11', 'Ol_13', 'Ol_13', 'Ol_19', 'Ol_19', 'Ol_25', 'Ol_25', 'Ol_32', 'Ol_32', 'Dr_8s', 'Dr_8s', 'Dr_48', 'Dr_48', 'Dr_epb', 'Dr_epb', 'Dr_dom', 'Dr_dom', 'Dr_shl', 'Dr_shl')
group21 <- c('Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur', 'Alim', 'Alim','Alim')
group22 <- c('Nf_Dia6', 'Nf_Dia6', 'Nf_esc', 'Nf_esc', 'Nf_KvO', 'Nf_KvY', 'Nf_KvY', 'Nf_DiaM', 'Nf_DiaM', 'Nf_DiaM', 'Nf_KvO', 'Al_DiaM', 'Al_esc', 'Al_Kv')
group31 <- c('Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat', 'Olat')
group32 <- c('Nf_Dia6', 'Nf_Dia6', 'Nf_esc', 'Nf_esc', 'Nf_KvO', 'Nf_KvY', 'Nf_KvY', 'Nf_DiaM', 'Nf_DiaM', 'Nf_DiaM', 'Nf_KvO','Ol_11', 'Ol_11', 'Ol_13', 'Ol_13', 'Ol_19', 'Ol_19', 'Ol_25', 'Ol_25', 'Ol_32', 'Ol_32')
group41 <- c('Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur','Nfur', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer', 'Drer')
group42 <- c('Nf_Dia6', 'Nf_Dia6', 'Nf_esc', 'Nf_esc', 'Nf_KvO', 'Nf_KvY', 'Nf_KvY', 'Nf_DiaM', 'Nf_DiaM', 'Nf_DiaM', 'Nf_KvO', 'Dr_8s', 'Dr_8s', 'Dr_48', 'Dr_48', 'Dr_epb', 'Dr_epb', 'Dr_dom', 'Dr_dom', 'Dr_shl', 'Dr_shl')
            
Fall <- cbind(group11, group12, Fall)
NfAl <- cbind(group21, group22, NfAl)
NfOl <- cbind(group31, group32, NfOl)
NfDr <- cbind(group41, group42, NfDr)


##################################################################################
#Visulaize PCA and output graphs to file
##################################################################################
pdf('PCA_Report.pdf') #Name of ouput file generated
autoplot(test1, data=Fall, size = 4, colour='group11')
autoplot(test1, data=Fall, size = 4, colour='group12')
autoplot(test2, data=NfAl, size = 4, colour='group21')
autoplot(test2, data=NfAl, size = 4, colour='group22')
autoplot(test3, data=NfOl, size = 4, colour='group31')
autoplot(test3, data=NfOl, size = 4, colour='group32')
autoplot(test4, data=NfDr, size = 4, colour='group41')
autoplot(test4, data=NfDr, size = 4, colour='group42')
dev.off()

