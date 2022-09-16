
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
#Calculates plaque diameter in mm
calc_plaque_size<-function(dead_cells,site_scale){
   return(sqrt(dead_cells)*site_scale/1000);
}


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################
args<-commandArgs(trailingOnly=T);

run<-args[2];
current_time<-args[3];
max_log_vl<-as.numeric(args[8]);
max_time<-as.numeric(args[9]);

png(paste(args[1],"snapshot_",run,"_",current_time,".png",sep=""),width = 8, height = 8, units="in",res=200);
par(mar=c(4,4,2,6), xpd=TRUE, oma=c(1,2,0,1),cex=0.8);
layout(matrix(c(1:4),2,2,byrow = TRUE),widths=c(1,1,1,1),heights=c(1,1,1,1));

site_scale<-50; #how many microns is one site

L<-125;

#read in cell_state matrix first
cell_file<-paste(args[1],"cell_state_",run,"_",current_time,".csv",sep="");
print(paste("Reading file",cell_file));
data <- read.csv(cell_file, comment.char="#");
cell_state <- data.matrix(data[,1:ncol(data)]);

susc<-sum(cell_state%in%c(0,1,2));
infn<-sum(cell_state%in%c(3,4,5,6));
infp<-sum(cell_state%in%c(7,8,9,10));

temp_grid<-cell_state;
temp_grid[cell_state%in%c(0,1,2)]<-0;
temp_grid[cell_state%in%c(3)]<-1;
temp_grid[cell_state%in%c(4)]<-2;
temp_grid[cell_state%in%c(5)]<-3;
temp_grid[cell_state%in%c(6)]<-4;
temp_grid[cell_state%in%c(7)]<-5;
temp_grid[cell_state%in%c(8)]<-6;
temp_grid[cell_state%in%c(9)]<-7;
temp_grid[cell_state%in%c(10)]<-8;
temp_grid[cell_state%in%c(11)]<-9;
cols<-c('papayawhip','pink','lightblue','lavender','yellow','red','blue','purple','orange','black');                     
leg_text<-c("Susc","preA","preB","preA|B","preAB","A","B","A|B","AB","Dead")
num_breaks<-9

image(temp_grid,col=cols,breaks=c(-1:num_breaks),xaxt='n',yaxt='n',main='Cells', asp=1)
#leg_text<-c("empty","susc","pre-prod","prod inf","died","cyt killed","Trm killed")
legend("topright", inset=c(-0.4, 0),leg_text,fill=cols,bty='n',ncol=1)

#Panel 2: Virions
#read in free_virus matrix
#0,A,B,ABn,ABr
col_virus<-c('papayawhip','red','blue','purple','orange','green');
cols<-col_virus;
vl_file<-paste(args[1],"virus_stateA_",run,"_",current_time,".csv",sep="");
print(paste("Reading file",vl_file));
data <- read.csv(vl_file, comment.char="#");
free_virusA <- data.matrix(data[,1:ncol(data)]);

vl_file<-paste(args[1],"virus_stateB_",run,"_",current_time,".csv",sep="");
print(paste("Reading file",vl_file));
data <- read.csv(vl_file, comment.char="#");
free_virusB <- data.matrix(data[,1:ncol(data)]);

vl_file<-paste(args[1],"virus_stateAB_",run,"_",current_time,".csv",sep="");
print(paste("Reading file",vl_file));
data <- read.csv(vl_file, comment.char="#");
free_virusAB <- data.matrix(data[,1:ncol(data)]);

virionsA<-sum(free_virusA);
virionsB<-sum(free_virusB);
virionsAB<-sum(free_virusAB);
virions<-sum(virionsA,virionsB,virionsAB);
if(virions==0){
  log_virus<-0
  image(free_virusA,col='papayawhip',xaxt='n',yaxt='n',main='Virus', asp=1)
}else{
  temp_virus<-free_virusA+free_virusB+free_virusAB
  log_virus<-log(virions)/log(10);
  temp_virus[(free_virusAB > 0 & (free_virusA >  0 | free_virusB > 0))]<-5;
  temp_virus[(free_virusAB > 0 & free_virusA == 0 & free_virusB == 0)]<-4;
  temp_virus[(free_virusAB == 0 & free_virusA > 0) & free_virusB > 0]<-3;
  temp_virus[(free_virusAB == 0 & free_virusA == 0) & free_virusB > 0]<-2;
  temp_virus[(free_virusAB == 0 & free_virusA > 0) & free_virusB == 0]<-1;
  temp_virus[(free_virusAB == 0 & free_virusA == 0) & free_virusB == 0]<-0;
  
  # image(free_virus,col=c('grey',colorRampPalette(c('blue','white'))(10)),xaxt='n',yaxt='n',main='Virus');
  #p<-sample(rainbow(7,start=0,end=6))
  #image(log(free_virus[(L/4):(3*L/4),(L/4):(3*L/4)])/log(10),col=cols,breaks=c(-1:8),xaxt='n',yaxt='n',main='Virus');
  image(temp_virus,col=cols,breaks=c(-1:5),xaxt='n',yaxt='n',main='Virus', asp=1);
  leg_text<-c('None','A','B','A&B','AB','AB&(A|B)')
  legend("topright", inset=c(-0.5, 0),leg_text,fill=cols,bty='n',ncol=1)
}

#Panel 3: Time history: cells
#read in results.csv
results_file<-paste(args[1],"results.csv",sep="");
#run,time,susceptible,inf_nonprod,infn_A,infn_B,infn_ABn,infn_ABr,inf_prod,infp_A,infp_B,infp_ABn,infp_ABr,dead_cells,virions,A virions,B virions,AB virions,log virions,max log VL,max time,viral_cells,first sampled log VL,first sampled time,max sampled log VL,max sampled time,plaque size,total infected,log total virus
print(paste("Reading file",results_file));
data <- read.csv(results_file, comment.char="#");

time <- data[,2];
infn_A <- 100*data[,5]/(L*L);
infn_B <- 100*data[,6]/(L*L);
infn_ABn <- 100*data[,7]/(L*L);
infn_ABr <- 100*data[,8]/(L*L);
infp_A <- 100*data[,10]/(L*L);
infp_B <- 100*data[,11]/(L*L);
infp_ABn <- 100*data[,12]/(L*L);
infp_ABr <- 100*data[,13]/(L*L);
dead <- 100*data[,14]/(L*L);

cols<-c('pink','lightblue','lavender','yellow','red','blue','purple','orange','black');                     

plot(time,infn_A,type="l",xlab="days",ylab="% cells",xlim=c(0,5),ylim=c(0,100),col=cols[1],lty=1);
lines(time,infn_B,col=cols[2],lty=2);
lines(time,infn_ABn,col=cols[3],lty=3);
lines(time,infn_ABr,col=cols[4],lty=4);
lines(time,infp_A,col=cols[5],lty=5);
lines(time,infp_B,col=cols[6],lty=6);
lines(time,infp_ABn,col=cols[7],lty=7);
lines(time,infp_ABr,col=cols[8],lty=8);
lines(time,dead,col=cols[9],lty=9);

leg_text<-c("% preA","% preB","% preA|B","% preAB","% A","% B","% A|B","% AB","% Dead")
legend("topright", inset=c(-0.5, 0),leg_text,col=cols,lty=c(1:9),bty='n',ncol=1)

#Panel 3: Time history: virus
#read in results.csv
virusA <- data[,16]
for (i in 1:length(virusA)) if (virusA[i] != 0) virusA[i] = log(virusA[i])/log(10);

virusB <- data[,17]
for (i in 1:length(virusB)) if (virusB[i] != 0) virusB[i] = log(virusB[i])/log(10);

virusAB <- data[,18]
for (i in 1:length(virusAB)) if (virusAB[i] != 0) virusAB[i] = log(virusAB[i])/log(10);

cols<-c('red','blue','orange');                     
plot(time,virusA,type="l",xlab="days",ylab="log virus",xlim=c(0,5),ylim=c(0,8),col=cols[1],lty=1);
lines(time,virusB,col=cols[2],lty=2);
lines(time,virusAB,col=cols[3],lty=3);
leg_text<-c("A","B","AB")
legend("topright", inset=c(-0.3, 0),leg_text,col=cols,lty=c(1:3),bty='n',ncol=1)

#Timestamp
mtext(paste('Episode',run,' Time =',formatC(as.numeric(current_time),format='f',digits=1),'days'),bg="white",side=1,outer=T);

dev.off()

print(paste("Inf cells =",infn+infp,"Log10 VL =",format(round(log_virus,2)),sep=" "));
