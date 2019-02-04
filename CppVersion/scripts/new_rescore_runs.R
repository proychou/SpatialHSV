
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

args<-commandArgs(trailingOnly=T);
beta<-as.numeric(args[1])
viral_prod<-as.numeric(args[2])
viral_diff<-as.numeric(args[3])
cyt_diff<-as.numeric(args[4])
cyt_uptake<-as.numeric(args[5])
beta_ic50<-as.numeric(args[6])
inf_ic50<-as.numeric(args[7])
prod_ic50<-as.numeric(args[8])
scale<-as.numeric(args[9])
cyt_trm_ic50<-as.numeric(args[10])
cyt_inf_death<-as.numeric(args[11])
hsv_fract<-as.numeric(args[14])

if (beta > 0) beta<-log(beta)/log(10)
if (viral_prod > 0) viral_prod<-log(viral_prod)/log(10)
if (viral_diff > 0) viral_diff<-log(viral_diff)/log(10)
if (cyt_diff > 0) cyt_diff<-log(cyt_diff)/log(10)
if (cyt_uptake > 0) cyt_uptake<-log(cyt_uptake)/log(10)
if (beta_ic50 > 0) beta_ic50<-log(beta_ic50)/log(10)
if (inf_ic50 > 0) inf_ic50<-log(inf_ic50)/log(10)
if (prod_ic50 > 0) prod_ic50<-log(prod_ic50)/log(10)
if (cyt_trm_ic50 > 0) cyt_trm_ic50<-log(cyt_trm_ic50)/log(10) 

model<-as.numeric(args[12])
k<-as.numeric(args[13])

write("model,k,log beta,log viral prod,log viral diff,log cyto diff,log cyto uptake,log inf ic50,log beta ic50,log prod ic50,log cyt_trm_ic50,cyt_inf_death,hsv_fract,peak VL RSS,peak VL AIC,peak time RSS,peak time AIC,duration RSS,duration AIC,cat2 RSS,cat2 AIC,cat3 RSS,cat3 AIC", file="scores.csv",append=FALSE)

cen_data <- read.csv("../../trial_data/hsv_cen_episodes.csv", comment.char="#")
#Episode,ID,duration,peak log VL,peak time,first pos VL,mean log VL
x=colnames(cen_data)
base_peak<-sort(as.numeric(cen_data[,4]),decreasing=TRUE)
base_peak[base_peak>7]<-7
base_peak_time<-sort(as.numeric(cen_data[,5]),decreasing=TRUE)
base_peak_time[base_peak_time>5]<-5
base_days<-sort(as.numeric(cen_data[,3]),decreasing=TRUE)
base_days[base_days>5]<-5
base_first<-sort(as.numeric(cen_data[,6]),decreasing=TRUE)

avg_peak <- 0
avg_peak_time <- 0
avg_days <- 0

n<-nrow(cen_data)
for (i in 1:n) {
    avg_peak = avg_peak + base_peak[i]
    avg_peak_time = avg_peak_time + base_peak_time[i]
    avg_days = avg_days + base_days[i]
}
avg_peak = avg_peak / n
avg_peak_time = avg_peak_time / n
avg_days = avg_days / n

sd_peak = sd(base_peak)
sd_peak_time = sd(base_peak_time)
sd_days = sd(base_days)
print(paste("Base data sd log peak is ",sd_peak))
print(paste("Base data sd peak time is ",sd_peak_time))
print(paste("Base data sd duration is ",sd_days))

print(paste("Censored data log avg peak is ",avg_peak))
print(paste("Censored data avg peak time is ",avg_peak_time))
print(paste("Censored data avg duration is ",avg_days))

run_data <- read.csv("all_episodes.csv", comment.char="#")
#infectivity,viral prod,viral diff,cyto diff,cyto uptake,beta ic50,inf ic50,beta ic50,scale_at_max,cyt_trm_ic50,cyt_inf_death,hsv_fract,Episode,run,duration,scaled peak logVL,peak time,first pos VL,mean VL,avg E/T ratio,3x3 E/T ratio,5x5 E/T ratio,7x7 E/T ratio,9x9 E/T ratio,11x11 E/T ratio,percent 0s,start HSV+,start BYST,nearest HSV+,nearest BYST,max cytokine,total infected,log total virus,cytokine kills,Tcell kills,Tcell act time
#1.377494e-01,100000,8.120941e+00,0,0,0,0,0,10,0,0,24,1,1,8.04656286031241,7.863,6.83799218782224,5.567,0.196,0,0.04,0.102,0.173,0.157,38.194,328,2500,2,3.162,282.884,16226,8.412,0,12322,0.25

scale_at_max<-mean(as.numeric(run_data[,9]))
scale_at_max<-1
days<-as.numeric(run_data[,15])
peak_vl<-as.numeric(run_data[,16])
peak_time<-as.numeric(run_data[,17])
first_vl<-as.numeric(run_data[,18])
#mean_vl<-as.numeric(run_data[,19])

days<-sort(days,decreasing=TRUE)
days[days>5]<-5
peak_vl<-sort(peak_vl,decreasing=TRUE)
peak_vl[peak_vl>7]<-7
peak_time<-sort(peak_time,decreasing=TRUE)
peak_time[peak_time>5]<-5
first_vl<-sort(first_vl,decreasing=TRUE)

peak_score <- 0
n<-min(nrow(cen_data),nrow(run_data))
for (i in 1:n) {
    peak_score = peak_score + ((peak_vl[i] - base_peak[i])/sd_peak)**2
}
peak_score=sqrt(peak_score)
print(paste("Model",model,"Weighted RSS Score from",n,"peak bins is",peak_score))
peak_AIC=n*log(peak_score/n)+2*k
print(paste("Model",model,"AIC Score from",n,"peak bins is",peak_AIC))

peak_time_score <- 0
for (i in 1:n) {
    peak_time_score = peak_time_score + ((peak_time[i] - base_peak_time[i])/sd_peak_time)**2
}
peak_time_score=sqrt(peak_time_score)
print(paste("Model",model,"RSS Score from",n,"peak time bins is",peak_time_score))
peak_time_AIC=n*log(peak_time_score/n)+2*k
print(paste("Model",model,"AIC Score from",n,"peak time bins is",peak_time_AIC))

days_score <- 0
for (i in 1:n) {
    days_score = days_score + ((days[i] - base_days[i])/sd_days)**2
}
days_score=sqrt(days_score)
print(paste("Model",model,"RSS Score from",n,"days bins is",days_score))
days_AIC=n*log(days_score/n)+2*k
print(paste("Model",model,"AIC Score from",n,"days bins is",days_AIC))

cat2_rss_score<-peak_score+peak_time_score
print(paste("Model",model,"2 Category RSS Score is",cat2_rss_score))

cat2_AIC_score<-n * log(peak_score/(n)) + n * log(peak_time_score/(n)) + 2 * k
print(paste("Model",model,"2 Category AIC Score is",cat2_AIC_score))

cat3_rss_score<-days_score+peak_score+peak_time_score
print(paste("Model",model,"3 Category RSS Score is",cat3_rss_score))

cat3_AIC_score<-n * log(peak_score/(n)) + n * log(peak_time_score/(n)) + n * log(days_score/(n)) + 2 * k
print(paste("Model",model,"3 Category AIC Score is",cat3_AIC_score))

write(paste(model,k,beta,viral_prod,viral_diff,cyt_diff,
	cyt_uptake,inf_ic50,beta_ic50,prod_ic50,cyt_trm_ic50,cyt_inf_death,hsv_fract,
	peak_score,peak_AIC, peak_time_score,peak_time_AIC,days_score,days_AIC,
	cat2_rss_score,cat2_AIC_score, cat3_rss_score,cat3_AIC_score,
	sep=","), file="scores.csv",append=TRUE)
