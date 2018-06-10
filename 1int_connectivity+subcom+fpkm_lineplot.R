library(parallel)
library(snow)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(assertthat)
library(grid)
library(scales)
library(gridExtra)
library(Hmisc)
cl <- makeCluster(rep("localhost", 8), type = "SOCK")
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
table_up <- args[1]
table_down <- args[2]

time.start=proc.time()[[3]]

#downstream##########################################
#table_down <- "D:/shai/hiC_chip-seq/intron_data/CTCF_intron_chr19.table_down.txt"
table1<-read.table(table_down, sep = "\t", header = FALSE)
#table1 <- table1[1:20,]
vars <- colsplit(table1$V4, "_", c("ENST", "subcomp" ,"FPKM","connctivity"))
table1 <-append(table1,vars,4)
table1<- as.data.frame(table1)
table1$V4 <- NULL
table1 <- table1[order(table1$ENST),]#order by ENST
#subseting only the bigwig data 
bigwig_d <- table1[,11:ncol(table1)]
end_loci <- length(bigwig_d)
#upstream#################################################
#table_up <- "D:/shai/hiC_chip-seq/intron_data/CTCF_intron_chr19.table_up.txt"
table1u<-read.table(table_up,sep = "\t")
#table1u <- table1u[1:20,]
varsu <- colsplit(table1u$V4, "_", c("ENST", "subcomp" ,"FPKM","connctivity"))
table1u <-append(table1u,varsu,4)
table1u<- as.data.frame(table1u)
table1u$V4 <- NULL
table1u <- table1u[ order(table1u$ENST),]#order by ENST
#subseting only the bigwig data 
bigwig_u <- table1u[,c(11:ncol(table1u))]
start_loci <- length(bigwig_u)-1
###########################################################
full_table <- cbind(table1u,bigwig_d)
position = as.numeric(c(-start_loci:end_loci))

#retriving the name
split1 <- strsplit(table_down, "\\.")[[1]]
split2 <- strsplit(split1[1], "/")[[1]]
name <- split2[length(split2)]

#avarage pval per ENST#####################################
#av_vector <- (apply(full_table[,c(511:1511)], 1, mean, na.rm=TRUE))
#av_table <- cbind.data.frame(full_table[,4],av_vector)
#colnames(av_table) <- c("ENST",name)
#rownames(av_table) <- NULL
#write.table(av_table, file = paste0(name,"_av_table.txt"), sep = "\t" , row.names=FALSE)

#CONNECTIVITY######################
binsConnectivity = 4
#ordering by connectivity, low to high
tb <- full_table[order(full_table$connctivity),]
#get the bin ranges
x <- cut2(tb$connctivity, g=binsConnectivity)
ranges_tb <- cbind.data.frame(x, tb)
ranges_tb$x <- as.character(ranges_tb$x)
###########################################################################
#mean loop
#uni_tb <- unique(ranges_tb[,c(2,3,4,8,12:ncol(ranges_tb))]) #chr-start-end-subcomp

con_bin_table <- data.frame(matrix(NA, nrow = c(0:4), ncol = length(position)))
for (i in as.character(levels(x))){
  print(i)
  indx <- which(ranges_tb$x ==  i)#the row numbers for the bin
  mean_bin_vector <- t(as.matrix(apply(ranges_tb[indx, 12:ncol(ranges_tb)], 2, mean, na.rm=TRUE)))#the vector of the mean per nuc
  bin_mean_bin_vector <- cbind(i,mean_bin_vector)
  colnames(bin_mean_bin_vector) <- c("Connectivity",position)
  con_bin_table <- rbind(con_bin_table,bin_mean_bin_vector)
}

qtr1 <- ranges_tb[which(ranges_tb$x ==  as.character(levels(x)[1])),]
qtr4 <- ranges_tb[which(ranges_tb$x ==  as.character(levels(x)[4])),]
mean_qtr1 <- t(as.matrix(apply(qtr1[, 12:ncol(ranges_tb)], 2, mean, na.rm=TRUE)))
mean_qtr4 <- t(as.matrix(apply(qtr4[, 12:ncol(ranges_tb)], 2, mean, na.rm=TRUE)))
max_mean_qtr1 <- (apply(qtr1[,which.max(mean_qtr1)-500:which.max(mean_qtr1)+500], 1, mean, na.rm=TRUE))
max_mean_qtr4 <- (apply(qtr4[,which.max(mean_qtr4)-500:which.max(mean_qtr4)+500], 1, mean, na.rm=TRUE))
result <- t.test(x = max_mean_qtr1 ,y = max_mean_qtr4)
p_v <- as.character(scientific(result$p.value, digits = 3))
t_v <- as.character(substr(result$statistic,1,5))

out <- paste0("t.v:",t_v ,", p.v:",p_v)

##########################################################################
'
#for the statistics########################################
all_mean_vector <- (apply(sub_tb[,c(5575:6575)], 1, mean, na.rm=TRUE))

tb_max <- cbind.data.frame(tb$ENST ,tb$connctivity,all_mean_vector)
colnames(tb_max) <- c("ENST", "Connectivity","Avarage_Pval")
#write.table(tb_max, file = paste0(name,"_max_nuc.txt"),   sep = "\t" )

ranges <- cut2(all_mean_vector,g = 4)

a <- cbind.data.frame(ranges,all_mean_vector)
colnames(a) <- c("bin","val")
a$bin <- as.character(a$bin)

qtr1_indx <- which(a$bin ==  as.character(levels(ranges)[1]))
qtr2_indx <- which(a$bin ==  as.character(levels(ranges)[2]))
qtr3_indx <- which(a$bin ==  as.character(levels(ranges)[3]))
qtr4_indx <- which(a$bin ==  as.character(levels(ranges)[4]))

qtr1 <- a[qtr1_indx,2]
qtr2 <- a[qtr2_indx,2]
qtr3 <- a[qtr3_indx,2]
qtr4 <- a[qtr4_indx,2]

result <- wilcox.test(x = qtr3 ,y = qtr4)
p_v <- as.character(scientific(result$p.value, digits = 3))
out <- paste0(name,"\t",p_v)
#setwd("/home/shaidulberg/chipseq/Modifications")
#write.table(out, file = paste0("pvalues.txt") , sep = "\t" , append = TRUE , eol = "\n",row.names = FALSE, col.names = FALSE)

n <- max(length(qtr1), length(qtr2), length(qtr3) ,length(qtr4))
length(qtr1) <- n   
length(qtr2) <- n 
length(qtr3) <- n 
length(qtr4) <- n

qtr1_table <- cbind.data.frame("1Qt",qtr1)
colnames(qtr1_table) <- c("Qt", "Pval")
qtr2_table <- cbind.data.frame("2Qt",qtr2)
colnames(qtr2_table) <- c("Qt", "Pval")
qtr3_table <- cbind.data.frame("3Qt",qtr3)
colnames(qtr3_table) <- c("Qt", "Pval")
qtr4_table <- cbind.data.frame("4Qt",qtr4)
colnames(qtr4_table) <- c("Qt", "Pval")

box_table <- rbind.data.frame(qtr1_table,qtr2_table, qtr3_table, qtr4_table)
# compute lower and upper whiskers
ylim1 = boxplot.stats(box_table$Pval)$stats[c(1, 5)]
ylim1[2] <- ylim1[2]+4

box <- ggplot(box_table,aes(y=Pval)) + 
geom_boxplot(aes(x = Qt),notch=TRUE, outlier.colour = NA)+ 
ggtitle(paste0("p=", p_v )) +
theme(legend.position="none") +
coord_cartesian(ylim = ylim1*1.8)+ # scale y limits based on ylim1
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(size = rel(1)))
'
###########################################################

con_bin_table$Connectivity <- c("Low","Medium","High","Highest")

mCon <- as.data.frame(melt(con_bin_table,id.vars="Connectivity"))#melting the table so each loci is infront of its bin
colnames(mCon) <-c("Connectivity","Loci","Pval")
mCon$Connectivity <- na.omit(as.character(mCon$Connectivity))
mCon$Loci <- na.omit(as.numeric(as.character(mCon$Loci)))
mCon$Pval <- na.omit(as.numeric(mCon$Pval))

#forcing the order of the legended titles
mCon$O.E <- factor(mCon$Connectivity, levels= c("Highest","High","Medium","Low"), labels=c("Highest","High","Medium","Low"))
pCON <-ggplot(mCon, aes(x=Loci, y=Pval, group=Connectivity)) +
  geom_line(aes(color=Connectivity), size=0.7) +
  scale_x_continuous(breaks = c(-5000,-75,5000) , labels = c(-5000,"TSS",5000)) +
  #coord_cartesian(ylim = c(0, 25)) +
  scale_fill_gradient(low="darkgreen",high="green") +
  geom_segment(aes(x=-75,xend=75,y=0,yend=0),lwd=4,color="black")+
  labs(y=paste0("ChIP-seq Signal"))+
  geom_segment(aes(x=-start_loci,xend=end_loci,y=0,yend=0),lwd=1,color="black")+
  scale_colour_manual(name='Connectivity', values=c("Low"="black", "Medium"= 'blue3', "High" = 'dodgerblue', "Highest" = 'deepskyblue'), guide='legend') +
  # theme(legend.position="none")+ #no legend
  #theme(legend.justification = c("right", "top"),legend.position = c(.95, .95))+ #option for legend on the figure
  ggtitle(name, subtitle = out)

setwd("/home/shaidulberg/chipseq/Modifications/1intron_connectivity_figuers")
#setwd("D:/shai/hiC_chip-seq/intron_data")
png(file= paste0(name,"_1int_connectivity.png"),width=850,height=600,res = 150)
print(pCON)
#print(grid.arrange(p, box ,nrow=1, ncol=2,newpage = TRUE, widths=c(3,1)))
dev.off()

#### fpkm ###################################################################
FPKMtb <- full_table[order(full_table$FPKM),] # sorting from the smallest to the biggest
binVector <- data.frame(matrix(NA, nrow = length(FPKMtb$FPKM), ncol = 1))
FPKMtb <- append(FPKMtb,binVector,1)
FPKMtb <- as.data.frame(FPKMtb)
colnames(FPKMtb)[2] <- "bin"
binsFPKM = 4
bin_0 <- which(FPKMtb$FPKM==0)
FPKMtb$bin[bin_0] <- 0 #setting the 0 bin
x <- (length(FPKMtb$FPKM)-length(bin_0))%%binsFPKM
if(x!=0){
  FPKMtb <-  FPKMtb[-c(1:x),] 
}

part <- ((length(FPKMtb$FPKM)-length(bin_0))/binsFPKM)#getting the size of the bins
# giving a bin number to the designated part of the table, besides the zerom bin 
for (i in 1:binsFPKM){
  for (line in (1+(i-1)*part):((i*part))) {
    FPKMtb$bin[length(bin_0)+line] <- i
  }
}

fpkm_bin_table <- data.frame(matrix(NA, nrow = 0:binsFPKM, ncol = length(position)))
for (y in 0:3){   
  #n = as.character(length(which(uni_tb$binsFPKM==(y))))
  #write(paste0(y," n=",n),"/home/shaidulberg/chipseq/Modifications/exons_n.txt",append = TRUE)
  mean_bin_vector_fpkm <-  apply(FPKMtb[which(FPKMtb$bin==(y)),12:ncol(FPKMtb)], 2, mean , na.rm = TRUE)
  mean_bin_vector_fpkm <-  append(mean_bin_vector_fpkm,y,0) 
  fpkm_bin_table <- rbind(fpkm_bin_table,mean_bin_vector_fpkm)
  colnames(fpkm_bin_table) <-c("bin",position)
}
fpkm_bin_table$bin <- c("Non","Low","Medium","High")

qtr1 <- FPKMtb[which(FPKMtb$bin ==  0),]
qtr4 <- FPKMtb[which(FPKMtb$bin ==  4),]
mean_qtr1 <- t(as.matrix(apply(qtr1[, 12:ncol(qtr1)], 2, mean, na.rm=TRUE)))
mean_qtr4 <- t(as.matrix(apply(qtr4[, 12:ncol(qtr4)], 2, mean, na.rm=TRUE)))
max_mean_qtr1 <- (apply(qtr1[,which.max(mean_qtr1)-500:which.max(mean_qtr1)+500], 1, mean, na.rm=TRUE))
max_mean_qtr4 <- (apply(qtr4[,which.max(mean_qtr4)-500:which.max(mean_qtr4)+500], 1, mean, na.rm=TRUE))
result <- t.test(x = max_mean_qtr1 ,y = max_mean_qtr4)
p_v <- as.character(scientific(result$p.value, digits = 3))
t_v <- as.character(substr(result$statistic,1,5))

out <- paste0("t.v:",t_v ,", p.v:",p_v)

mFPKM <- melt(fpkm_bin_table,id.vars="bin")#melting the table so each loci is infront of its bin
colnames(mFPKM) <-c("FPKM","Loci","ChIP_seq_Signal")
mFPKM$FPKM <- factor(mFPKM$FPKM, levels= c("High","Medium","Low","Non"), labels=c("High","Medium","Low","Non"))
mFPKM$Loci <- as.numeric(as.character(mFPKM$Loci))
mFPKM$ChIP_seq_Signal <- as.numeric(mFPKM$ChIP_seq_Signal)

pFPKM <- ggplot(mFPKM, aes(x = Loci, y = ChIP_seq_Signal, group = FPKM)) +
  geom_line(aes(color = FPKM), size = 0.7) +
  scale_x_continuous(breaks = c(-5000,-75,5000) , labels = c(-5000,"TSS",5000)) +
  #coord_cartesian(ylim = c(0, 25)) +
  geom_segment(aes(x=-75,xend=75,y=0,yend=0),lwd=4,color="black")+
  labs(y=paste0("ChIP-seq Signal"))+
  geom_segment(aes(x=-start_loci,xend=end_loci,y=0,yend=0),lwd=1,color="black")+
  scale_colour_manual(name='FPKM', values=c("Non"='black', "Low"='saddlebrown', "Medium" = 'red3', "High" = "red1"), guide='legend') +
  #theme(legend.position="none")+ #no legend
  #theme(legend.justification = c("right", "top"),legend.position = c(.95, .95))+ #option for legend on the figure
  ggtitle(name, subtitle = out)

setwd("/home/shaidulberg/chipseq/Modifications/1intron_FPKM_figures")
#setwd("D:/shai/hiC_chip-seq/intron_data")
png(file= paste0(name,"_1intron_FPKM.png"),width=850,height=600,res = 150)
print(pFPKM)
#print(grid.arrange(p, box ,nrow=1, ncol=2,newpage = TRUE, widths=c(3,1)))
dev.off()
#########################################################################

#SUBCOMPARTMENT##########################################################
sub_bin_table <- data.frame(matrix(NA, nrow = c(0:6), ncol = length(position+1)))
binsSub <- na.omit(unique(full_table$subcomp))
for (i in binsSub){
  #n = as.character(length(which(subset_table[,1]==(i))))
  #write(paste0(i," n=",n),"/home/shaidulberg/chipseq/Modifications/exons_n.txt",append = TRUE)
  indx <- which(full_table$subcomp == (i))#the row numbers in the subcompartment
  mean_vector <- t(as.matrix(apply(full_table[indx,11:ncol(full_table)], 2, mean,na.rm=TRUE)))#the vector of the mean of the subcomp, per nuc
  bin_mean_vector <- cbind(i,mean_vector)
  colnames(bin_mean_vector) <-c("Subcomp",position)
  sub_bin_table <- rbind(sub_bin_table,bin_mean_vector)
}

mSub <- as.data.frame(melt(sub_bin_table,id.vars="Subcomp"))#melting the table so each loci is infront of its bin
colnames(mSub) <- c("Subcomp","Loci","Histone_Modification")
mSub$Loci <- as.numeric(as.character(mSub$Loci))
mSub$Histone_Modification <- as.numeric(as.character(mSub$Histone_Modification))

pSub <- ggplot(mSub, aes(x=Loci, y=Histone_Modification, group=Subcomp)) +
  geom_line(aes(color=Subcomp), size=0.7) +
  scale_x_continuous(breaks = c(-5000,-75,5000) , labels = c(-5000,"TSS",5000)) +
  #coord_cartesian(ylim = c(0, 25)) +
  labs(y=paste0("ChIP-seq Signal"))+
  geom_segment(aes(x=-75,xend=75,y=0,yend=0),lwd=4,color="black")+
  geom_segment(aes(x=-start_loci,xend=end_loci,y=0,yend=0),lwd=1,color="black")+
  scale_colour_manual(name= "Subcomp", values= c("A1"="red","A2"="blue","B1"="gold1","B2"="purple","B3"="green4","B4"="cyan"), guide='legend') +
  ggtitle(name)

setwd("/home/shaidulberg/chipseq/Modifications/1intron_subcom_figures")
#setwd("D:/shai/hiC_chip-seq/intron_data")
png(file= paste0(name,"_1intron_subcom.png"),width=850,height=600,res = 150)
print(pSub)
#print(grid.arrange(p, box ,nrow=1, ncol=2,newpage = TRUE, widths=c(3,1)))
dev.off()

stopCluster(cl)
time.end=proc.time()[[3]]
time.end-time.start
