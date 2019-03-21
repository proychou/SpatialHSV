
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

data <- read.csv("trial_data/hsv_cen_episodes.csv", comment.char="#")
#Episode,ID,duration,peak log VL,peak time,first pos VL,mean log VL
x=colnames(data)
#data<-data[1:100,]
base_peak<-sort(as.numeric(data[,4]),decreasing=TRUE)
base_peak[base_peak>7]<-7
base_peak_time<-sort(as.numeric(data[,5]),decreasing=TRUE)
base_peak_time[base_peak_time>5]<-5
base_days<-sort(as.numeric(data[,3]),decreasing=TRUE)
base_days[base_days>5]<-5
base_first<-sort(as.numeric(data[,6]),decreasing=TRUE)

avg_peak <- 0
avg_peak_time <- 0
avg_days <- 0

n<-nrow(data)
for (i in 1:n) {
    avg_peak = avg_peak + base_peak[i]
    avg_peak_time = avg_peak_time + base_peak_time[i]
    avg_days = avg_days + base_days[i]
}
avg_peak = avg_peak / n
avg_peak_time = avg_peak_time / n
avg_days = avg_days / n
print(paste("Base data log avg peak is ",avg_peak))
print(paste("Base data avg peak time is ",avg_peak_time))
print(paste("Base data avg duration is ",avg_days))

run_data <- read.csv("best_fine_tune_epis.csv", comment.char="#")
run_data<-run_data[1:n,]
all_days<-as.numeric(run_data[,15])
all_days[all_days>5]<-5
all_peak_vl<-as.numeric(run_data[,16])
all_peak_vl[all_peak_vl>7]<-7
all_peak_time<-as.numeric(run_data[,17])
all_peak_time[all_peak_time>5]<-5

all_days<-sort(all_days,decreasing=TRUE)
all_peak_vl<-sort(all_peak_vl,decreasing=TRUE)
all_peak_time<-sort(all_peak_time,decreasing=TRUE)

run_data <- read.csv("sd_cat3_epis_model16.csv", comment.char="#")
run_data<-run_data[1:n,]
van_days<-as.numeric(run_data[,15])
van_days[van_days>5]<-5
van_peak_vl<-as.numeric(run_data[,16])
van_peak_time<-as.numeric(run_data[,17])
van_peak_time[van_peak_time>5]<-5

van_days<-sort(van_days,decreasing=TRUE)
van_peak_vl<-sort(van_peak_vl,decreasing=TRUE)
van_peak_time<-sort(van_peak_time,decreasing=TRUE)

png("peak_vl_base16.png")
data3<-data.frame(base_peak,all_peak_vl,van_peak_vl)
plot(data3$base_peak,ylim=c(0,round(max(data3$base_peak,data3$all_peak_vl,data3$van_peak_vl))),xlab="Rank",ylab = "Log VL",type="l",pch=19,col="blue",lty=1)
lines(data3$all_peak_vl,type="l",pch=18,lty=2,col="green")
lines(data3$van_peak_vl,type="l",pch=19,lty=3,col="red")
legend(x=nrow(run_data)/2,y=max(data3$base_peak,data3$all_peak_vl,data3$van_peak_vl),legend=c("Cohort","Fine Tuned","Best 16"),col=c("blue","green","red"),lty=1:6,cex=0.8)
title("Peak Viral Load")
dev.off()

png("peak_time_base16.png")
data3<-data.frame(base_peak_time,all_peak_time,van_peak_time)
plot(data3$base_peak_time,ylim=c(0,round(max(data3$base_peak_time,data3$all_peak_time,data3$van_peak_time))),xlab="Rank",ylab = "Days",type="l",pch=19,col="blue",lty=1)
lines(data3$all_peak_time,type="l",pch=18,lty=2,col="green")
lines(data3$van_peak_time,type="l",pch=19,lty=3,col="red")
legend(x=nrow(run_data)/2,y=max(data3$base_peak_time,data3$all_peak_time,data3$van_peak_time),legend=c("Cohort","Fine Tuned","Best 16"),col=c("blue","green","red"),lty=1:6,cex=0.8)
title("Time to Peak Viral Load")
dev.off()

png("duration_base16.png")
data3<-data.frame(base_days,all_days,van_days)
plot(data3$base_days,ylim=c(0,round(max(data3$base_days,data3$all_days,data3$van_days))),xlab="Rank",ylab = "Days",type="l",pch=19,col="blue",lty=1)
lines(data3$all_days,type="l",pch=18,lty=2,col="green")
lines(data3$van_days,type="l",pch=19,lty=3,col="red")
legend(x=nrow(run_data)/2,y=max(data3$base_days,data3$all_days,data3$van_days),legend=c("Cohort","Fine Tuned","Best 16"),col=c("blue","green","red"),lty=1:6,cex=0.8)
title("Episode Duration")
dev.off()
