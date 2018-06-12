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
#retriving the name
split1 <- strsplit(table_down, "\\.")[[1]]
split2 <- strsplit(split1[1], "/")[[1]]
name <- split2[length(split2)]

#downstream##########################################
table_down <- "D:/shai/hiC_chip-seq/intron_data/CTCF_intron_chr19.table_down.txt"
table1<-read.table(table_down, sep = "\t", header = FALSE)
table1 <- table1[1:200,]
vars <- colsplit(table1$V4, "_", c("ENST", "subcomp" ,"FPKM","connctivity"))
table1 <-append(table1,vars,4)
table1<- as.data.frame(table1)
table1$V4 <- NULL
table1 <- table1[order(table1$ENST),]#order by ENST
#subseting only the bigwig data 
bigwig_d <- table1[,11:ncol(table1)]
end_loci <- length(bigwig_d)
#upstream#################################################
table_up <- "D:/shai/hiC_chip-seq/intron_data/CTCF_intron_chr19.table_up.txt"
table1u<-read.table(table_up,sep = "\t")
table1u <- table1u[1:200,]
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


FPKMtb <- full_table[order(full_table$FPKM),] # sorting from the smallest to the biggest
ranges <- cut2(FPKMtb$FPKM, g=4)
ranges_FPKMtb <- cbind.data.frame(ranges, FPKMtb)
ranges_FPKMtb$ranges <- as.character(ranges_FPKMtb$ranges)
lowFPKMtb <- FPKMtb[c(which(ranges_FPKMtb$ranges ==  as.character(levels(ranges)[1]))),]
highFPKMtb <- FPKMtb[c(which(ranges_FPKMtb$ranges ==  as.character(levels(ranges)[4]))),]


#CONNECTIVITY - low FPKM######################
binsConnectivity = 4

figCreator <- function(table_input){
#ordering by connectivity, low to high
tb <- table_input[order(table_input$connctivity),]
#get the bin ranges
x <- cut2(tb$connctivity, g=binsConnectivity)
ranges_FPKMtb <- cbind.data.frame(x, tb)
ranges_FPKMtb$x <- as.character(ranges_FPKMtb$x)
###########################################################################
#mean loop
#uni_tb <- unique(ranges_tb[,c(2,3,4,8,12:ncol(ranges_tb))]) #chr-start-end-subcomp

con_bin_table <- data.frame(matrix(NA, nrow = c(0:4), ncol = length(position)))
for (i in as.character(levels(x))){
  print(i)
  indx <- which(ranges_FPKMtb$x == i)#the row numbers for the bin
  mean_bin_vector <- t(as.matrix(apply(ranges_FPKMtb[indx, c(12:ncol(ranges_FPKMtb))], 2, mean, na.rm=TRUE)))#the vector of the mean per nuc
  bin_mean_bin_vector <- cbind(i,mean_bin_vector)
  colnames(bin_mean_bin_vector) <- c("Connectivity",position)
  con_bin_table <- rbind(con_bin_table,bin_mean_bin_vector)
}

qtr1 <- ranges_FPKMtb[which(ranges_FPKMtb$x ==  as.character(levels(x)[1])),]
qtr4 <- ranges_FPKMtb[which(ranges_FPKMtb$x ==  as.character(levels(x)[4])),]
mean_qtr1 <- (apply(qtr1[,575:1575], 1, mean, na.rm=TRUE))
mean_qtr4 <- (apply(qtr4[,575:1575], 1, mean, na.rm=TRUE))
result <- t.test(x = mean_qtr1 ,y = mean_qtr4)
p_v <- as.character(scientific(result$p.value, digits = 3))
t_v <- as.character(substr(result$statistic,1,5))

out <- paste0("t.v:",t_v ,", p.v:",p_v)

con_bin_table$Connectivity <- c("Low","Medium","High","Highest")

mCon <- as.data.frame(melt(con_bin_table,id.vars="Connectivity"))#melting the table so each loci is infront of its bin
colnames(mCon) <-c("Connectivity","Loci","Pval")
mCon$Connectivity <- na.omit(as.character(mCon$Connectivity))
mCon$Loci <- na.omit(as.numeric(as.character(mCon$Loci)))
mCon$Pval <- na.omit(as.numeric(mCon$Pval))

#forcing the order of the legended titles
mCon$Connectivity <- factor(mCon$Connectivity, levels= c("Highest","High","Medium","Low"), labels=c("Highest","High","Medium","Low"))
pCON <- ggplot(mCon, aes(x=Loci, y=Pval, group=Connectivity)) +
  geom_line(aes(color=Connectivity), size=0.6) +
  scale_x_continuous(breaks = c(-5000,-75,5000) , labels = c(-5000,"TSS",5000)) +
  #coord_cartesian(ylim = c(0, 25)) +
  scale_fill_gradient(low="darkgreen",high="green") +
  geom_segment(aes(x=-75,xend=75,y=0,yend=0),lwd=4,color="black")+
  labs(y=paste0("ChIP-seq Signal"),subtitle = out)+
  geom_segment(aes(x=-start_loci,xend=end_loci,y=0,yend=0),lwd=1,color="black")+
  scale_colour_manual(name='Connectivity', values=c("Low"="black", "Medium"= 'blue3', "High" = 'dodgerblue', "Highest" = 'deepskyblue'), guide='legend') 
  #theme(legend.position="none") #no legend
  #theme(legend.justification = c("right", "top"),legend.position = c(.95, .95))+ #option for legend on the figure
  return(pCON)
}

pLOW <- figCreator(lowFPKMtb)
pHIGH <- figCreator(highFPKMtb)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pLOW)

#setwd("/home/shaidulberg/chipseq/Modifications/1intron_FPKM_figures")
setwd("D:/shai/hiC_chip-seq/intron_data")
png(file= paste0(name,"_1intron_H.L-FPKM.png"),width=1050,height=600,res = 170)
#print(pFPKM)
print(grid.arrange(pLOW+theme(legend.position="none")+labs(title = "Low Expression"),
                   pHIGH+theme(legend.position="none")+labs(title = "High Expression")+labs(y= NULL), mylegend,
                   nrow = 1,ncol=3,widths=c(3,3,1), newpage = TRUE,
                   top = textGrob(name ,gp=gpar(fontsize=20,font=2))))
dev.off()