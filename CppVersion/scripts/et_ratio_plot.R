
#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("lattice")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################
args<-commandArgs(trailingOnly=T);

run<-args[2];
time<-args[3];

png(paste(args[1],"et_ratio_",run,"_",time,".png",sep=""));
par(oma=c(1,2,1,1));
#layout(matrix(c(1:4),2,2,byrow = TRUE),widths=c(1,1,1,1),heights=c(1,1,1,1));

L<-125;

#read in cell_state matrix first
infile<-paste(args[1],"et_ratio_",run,"_",time,".csv",sep="");
print(paste("Reading file",infile));
data <- read.csv(infile, comment.char="#");
cluster <- as.numeric(data.matrix(data[,1]));
ratio <- sort(as.numeric(data.matrix(data[,2])),decreasing=TRUE);
plot(cluster,ratio,xlab="Cluster Rank",ylab="E/T ratio",ylim=c(0,1),col=c("blue"),main=paste('Episode',run,' Time =',time,'days'));

dev.off()
