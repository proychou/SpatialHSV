
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

peak_aic<-vector(length=16)
peak_time_aic<-vector(length=16)
duration_aic<-vector(length=16)
cat3_aic<-vector(length=16)

for (model in 1:16)
{
    all_data <- read.csv(paste("sd_top_final_scores_model",model,".csv",sep=""), header=FALSE, comment.char="#")
#model,k,log beta,log viral prod,log viral diff,log cyto diff,log cyto uptake,log inf ic50,log beta ic50,log prod ic50,log cyt_trm_ic50,cyt_inf_death,hsv_fract,peak VL RSS,peak VL AIC,peak time RSS,peak time AIC,duration RSS,duration AIC,cat2 RSS,cat2 AIC,cat3 RSS,cat3 AIC
    rows<-nrow(all_data)
    if (model==1) 
	min_rows<-rows
    else
	min_rows<-min(min_rows,rows)

    if (model == 1) {
	peak_aic<-all_data[,15]
	peak_time_aic<-all_data[,17]
	duration_aic<-all_data[,19]
	cat3_aic<-all_data[,23]
    } else {
	if (model == 2) {
	    peak_aic<-cbind(peak_aic[1:min_rows],all_data[1:min_rows,15])
	    peak_time_aic<-cbind(peak_time_aic[1:min_rows],all_data[1:min_rows,17])
	    duration_aic<-cbind(duration_aic[1:min_rows],all_data[1:min_rows,19])
	    cat3_aic<-cbind(cat3_aic[1:min_rows],all_data[1:min_rows,23])
	} else {
	    peak_aic<-cbind(peak_aic[1:min_rows,],as.matrix(all_data[1:min_rows,15]))
	    peak_time_aic<-cbind(peak_time_aic[1:min_rows,],as.matrix(all_data[1:min_rows,17]))
	    duration_aic<-cbind(duration_aic[1:min_rows,],as.matrix(all_data[1:min_rows,19]))
	    cat3_aic<-cbind(cat3_aic[1:min_rows,],as.matrix(all_data[1:min_rows,23]))
	}
    }
}

x=seq.int(1,16,1);
x=c("none","A","B","C","D","AB","AC","AD","BC","BD","CD","ABC","ABD","ACD","BCD","all")
#1) BASELINE
#2) BASELINE + 1
#3) BASELINE + 2
#4) BASELINE + 3
#5) BASELINE + 4
#6) BASELINE + 1 + 2
#7) BASELINE + 1 + 3
#8) BASELINE + 1 + 4
#9) BASELINE + 2 + 3
#10) BASELINE + 2 + 4
#11) BASELINE + 3 + 4
#12) BASELINE + 1 + 2 + 3
#13) BASELINE + 1 + 2 + 4
#14) BASELINE + 1 + 3 + 4
#15) BASELINE + 2 + 3 + 4
#16) BASELINE + 1 + 2 + 3 + 4
par()
png("peak_aic.png");
boxplot(peak_aic,col=rainbow(16),beside=TRUE, ylim=c(-1000,0),names=x,xlab="Model",ylab = "AIC Score",cex.axis=0.8)
title("Peak Viral Load AIC Scores\n(Top 5% of fits)");
dev.off()

png("peak_time_aic.png");
boxplot(peak_time_aic,beside=TRUE,col=rainbow(16),ylim=c(-1000,0),names=x,xlab="Model",ylab = "AIC Score",cex.axis=0.8)
title("Time to Peak AIC Scores\n(Top 5% of fits)");
dev.off()

png("duration_aic.png");
boxplot(duration_aic,beside=TRUE,col=rainbow(16),ylim=c(-1000,0),names=x,xlab="Model",ylab = "AIC Score",cex.axis=0.8)
title("Duration AIC Scores\n(Top 5% of fits)");
dev.off()

png("cat3_aic.png");
boxplot(cat3_aic,beside=TRUE,col=rainbow(16),ylim=c(-1000,0),names=x,xlab="Model",ylab = "AIC Score",cex.axis=0.8)
title("All Category AIC Scores\n(Top 5% of fits)");
dev.off()
