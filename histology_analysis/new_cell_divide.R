library(ggplot2)
library(gplots)
set.seed(42)
do_plots = 0
data <- read.table("base_names.txt", header=FALSE,comment.char="#")
file_names<-data[,1]
p_cd8_coefs<-vector(length=length(file_names))
c_cd8_coefs<-vector(length=length(file_names))
p_combo_coefs<-vector(length=length(file_names))
c_combo_coefs<-vector(length=length(file_names))
r2_combo<-vector(length=length(file_names))

data <- read.csv("cropped_cells.csv", comment.char="#")
#ImageNumber,ObjectNumber,AreaShape_Area,AreaShape_Center_X,AreaShape_Center_Y,AreaShape_Compactness,AreaShape_Eccentricity,AreaShape_EulerNumber,AreaShape_Extent,AreaShape_FormFactor,AreaShape_MajorAxisLength,AreaShape_MaxFeretDiameter,AreaShape_MaximumRadius,AreaShape_MeanRadius,AreaShape_MedianRadius,AreaShape_MinFeretDiameter,AreaShape_MinorAxisLength,AreaShape_Orientation,AreaShape_Perimeter,AreaShape_Solidity,Location_Center_X,Location_Center_Y,Number_Object_Number

cluster_size<-100

N<-nrow(data)
images=as.numeric(data[,1])
xs=as.numeric(data[,4])
ys=as.numeric(data[,5])
curr_image=images[1]
m<-matrix(data=NA,nrow=N,ncol=11)
objs=0
pobjs=1
tot_clusters=0
#pdf("cell_clustering.pdf")
for (i in 1:N)
{
    if (images[i] == curr_image)
    {
	objs = objs+1
    }
    else
    {
	x<-cbind(xs[pobjs:objs],ys[pobjs:objs])
	colnames(x) <- c("x", "y")
	tot_objs = objs - pobjs + 1
	clusters=as.integer(tot_objs/cluster_size)
	(cl <- kmeans(x, clusters))
	for (j in 1:clusters)
	{
	    m[tot_clusters+j,1]=curr_image
	    m[tot_clusters+j,2]=j
	    m[tot_clusters+j,3]=cl$centers[j,1]
	    m[tot_clusters+j,4]=cl$centers[j,2]
	    m[tot_clusters+j,5]=cl$cluster[j]
	    m[tot_clusters+j,6]=cl$size[j]	#Nuclei
	    m[tot_clusters+j,7]=0		#CD4s
	    m[tot_clusters+j,8]=0		#CD8s
	    #m[tot_clusters+j,8]=0
	    m[tot_clusters+j,9]=j
	    m[tot_clusters+j,10]=max(x[,1])
	    m[tot_clusters+j,11]=max(x[,2])
	}
	tot_clusters=tot_clusters+clusters

	print(paste("Plotting ",tot_objs,"cells and",clusters,"clusters from file",file_names[curr_image]))
	print(paste("First index",pobjs,"at (",xs[pobjs],",",ys[pobjs],"), last index",objs,"at (",xs[objs],",",ys[objs],")"))
	pdf(paste(file_names[curr_image],"new_nuclei_clustering.pdf",sep="_"))
	plot(x, xlim=c(0,max(x[,1])),ylim=c(0,max(x[,2])),col = cl$cluster, main=paste(file_names[curr_image],"(",clusters,"Clusters )"))
	points(cl$centers, col=1:clusters,pch = 8, cex = 2)
	dev.off()
	#pobjs=objs+1
	pobjs=i
	objs = i
	curr_image=images[i]
    }
}
x<-cbind(xs[pobjs:objs],ys[pobjs:objs])
colnames(x) <- c("x", "y")
tot_objs = objs - pobjs + 1
clusters=as.integer(tot_objs/cluster_size)
(cl <- kmeans(x, clusters))
for (j in 1:clusters)
{
    m[tot_clusters+j,1]=curr_image
    m[tot_clusters+j,2]=j
    m[tot_clusters+j,3]=cl$centers[j,1]
    m[tot_clusters+j,4]=cl$centers[j,2]
    m[tot_clusters+j,5]=cl$cluster[j]
    m[tot_clusters+j,6]=cl$size[j]	#Nuclei
    m[tot_clusters+j,7]=0		#CD4s
    m[tot_clusters+j,8]=0		#CD8s
    m[tot_clusters+j,9]=j
    m[tot_clusters+j,10]=max(x[,1])
    m[tot_clusters+j,11]=max(x[,2])
}
tot_clusters=tot_clusters+clusters
print(paste("Plotting ",tot_objs,"cells and",clusters,"clusters from file",file_names[curr_image]))
print(paste("First index",pobjs,"at (",xs[pobjs],",",ys[pobjs],"), last index",objs,"at (",xs[objs],",",ys[objs],")"))
print(paste("First cluster is",cl$cluster[1]))
print(paste("First nuc cluster is",cl$cluster[1],"at",xs[pobjs],",",ys[pobjs]))
pdf(paste(file_names[curr_image],"new_nuclei_clustering.pdf",sep="_"))
plot(x, xlim=c(0,max(x[,1])),ylim=c(0,max(x[,2])),col = cl$cluster, main=paste(file_names[curr_image],"(",clusters,"Clusters )"))
points(cl$centers, col=1:clusters,pch = 8, cex = 2)
m<-m[1:tot_clusters,]
dev.off()

data2 <- read.csv("cropped_CD4s.csv", comment.char="#")
CD4s<-nrow(data2)
images=as.numeric(data2[,1])
cd4_obj_num=as.numeric(data2[,2])
cd4_xs=as.numeric(data2[,4])
cd4_ys=as.numeric(data2[,5])
cd4_clusters<-vector(length=CD4s)
cd4_colors<-vector(length=CD4s)
cd4_maxx<-vector(length=CD4s)
cd4_maxy<-vector(length=CD4s)
cd4_centers<-matrix(data=NA,nrow=CD4s,ncol=2)
curr_image=images[1]
objs=0
pobjs=1
rep_clusters=0
print(paste("Read in ",CD4s,"CD4s"))
#pdf("cd4_clustering.pdf")
for (i in 1:CD4s)
{
    cluster_set=0
    min_cluster_dist=100000
    for (j in 1:tot_clusters)
    {
	if (images[i] == m[j,1])
	{
	    cluster_dist = sqrt((cd4_xs[i]-m[j,3])*(cd4_xs[i]-m[j,3]) +
				(cd4_ys[i]-m[j,4])*(cd4_ys[i]-m[j,4]))
	    if (cluster_dist < min_cluster_dist)
	    {
		cd4_clusters[i]=j
		min_cluster_dist=cluster_dist
		cluster_set=1
	    }
	}
    }
    if (cluster_set==0)
    {
	cd4_clusters[i]=1
	print(paste("Failed to find color for object",cd4_obj_num[i],"at (",cd4_xs[i],",",cd4_ys[i],") in",file_names[images[i]]))
    }
    else
    {
	if(m[cd4_clusters[i],7]==0)
	{
	    rep_clusters=rep_clusters+1
	    #m[cd4_clusters[i],9]= rep_clusters
	}
	m[cd4_clusters[i],7]= m[cd4_clusters[i],7]+1
	cd4_centers[i,1]=as.integer(m[cd4_clusters[i],3])
	cd4_centers[i,2]=as.integer(m[cd4_clusters[i],4])
	print(paste("For object",cd4_obj_num[i],"at (",cd4_xs[i],",",cd4_ys[i],") cluster",cd4_clusters[i],"center (",cd4_centers[i,1],",",cd4_centers[i,2],") in",file_names[images[i]]))
	cd4_colors[i]=m[cd4_clusters[i],9]
	cd4_maxx[i]=m[cd4_clusters[i],10]
	cd4_maxy[i]=m[cd4_clusters[i],11]
    }

    if (images[i] == curr_image)
    {
	objs = objs+1
    }
    else
    {
	x<-cbind(cd4_xs[pobjs:objs],cd4_ys[pobjs:objs])
	colnames(x) <- c("x", "y")
	tot_objs = objs - pobjs + 1
	pdf(paste(file_names[curr_image],"new_CD4_clustering.pdf",sep="_"))
	print(paste("Plotting CD4s from file",file_names[curr_image]))
	print(paste("First index",pobjs,"at (",cd4_xs[pobjs],",",cd4_ys[pobjs],"), last index",objs,"at (",cd4_xs[objs],",",cd4_ys[objs],")"))
	plot(x, xlim=c(0,cd4_maxx[pobjs]),ylim=c(0,cd4_maxy[pobjs]),col = cd4_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD4s from",rep_clusters,"clusters)"))
	points(x=cd4_centers[pobjs:objs,1],y=cd4_centers[pobjs:objs,2], col=cd4_colors[pobjs:objs],pch = 8, cex = 2)
	dev.off()
	#pobjs=objs+1
	pobjs=i
	objs = i
	curr_image=images[i]
	rep_clusters=0
    }
}
x<-cbind(cd4_xs[pobjs:objs],cd4_ys[pobjs:objs])
colnames(x) <- c("x", "y")
tot_objs = objs - pobjs + 1
print(paste("Plotting CD4s from file",file_names[curr_image]))
print(paste("First index",pobjs,"at (",cd4_xs[pobjs],",",cd4_ys[pobjs],"), last index",objs,"at (",cd4_xs[objs],",",cd4_ys[objs],")"))
pdf(paste(file_names[curr_image],"new_CD4_clustering.pdf",sep="_"))
plot(x, xlim=c(0,cd4_maxx[pobjs]),ylim=c(0,cd4_maxy[pobjs]),col = cd4_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD4s from",rep_clusters,"clusters)"))
points(x=cd4_centers[pobjs:objs,1],y=cd4_centers[pobjs:objs,2], col=cd4_colors[pobjs:objs],pch = 8, cex = 2)
dev.off()

data2 <- read.csv("cropped_CD8s.csv", comment.char="#")
CD8s<-nrow(data2)
images=as.numeric(data2[,1])
cd8_obj_num=as.numeric(data2[,2])
cd8_xs=as.numeric(data2[,4])
cd8_ys=as.numeric(data2[,5])
cd8_clusters<-vector(length=CD8s)
cd8_colors<-vector(length=CD8s)
cd8_maxx<-vector(length=CD8s)
cd8_maxy<-vector(length=CD8s)
cd8_centers<-matrix(data=NA,nrow=CD8s,ncol=2)
curr_image=images[1]
objs=0
pobjs=1
rep_clusters=0
print(paste("Read in ",CD8s,"CD8s"))
#pdf("cd8_clustering.pdf")
for (i in 1:CD8s)
{
    cluster_set=0
    min_cluster_dist=100000
    for (j in 1:tot_clusters)
    {
	if (images[i] == m[j,1])
	{
	    cluster_dist = sqrt((cd8_xs[i]-m[j,3])*(cd8_xs[i]-m[j,3]) +
				(cd8_ys[i]-m[j,4])*(cd8_ys[i]-m[j,4]))
	    if (cluster_dist < min_cluster_dist)
	    {
		cd8_clusters[i]=j
		min_cluster_dist=cluster_dist
		cluster_set=1
	    }
	}
    }
    if (cluster_set==0)
    {
	cd8_clusters[i]=1
	print(paste("Failed to find color for object",cd8_obj_num[i],"at (",cd8_xs[i],",",cd8_ys[i],") in",file_names[images[i]]))
    }
    else
    {
	if(m[cd8_clusters[i],8]==0)
	{
	    rep_clusters=rep_clusters+1
	    #m[cd8_clusters[i],8]= rep_clusters
	}
	m[cd8_clusters[i],8]= m[cd8_clusters[i],8]+1
	cd8_centers[i,1]=as.integer(m[cd8_clusters[i],3])
	cd8_centers[i,2]=as.integer(m[cd8_clusters[i],4])
	print(paste("For object",cd8_obj_num[i],"at (",cd8_xs[i],",",cd8_ys[i],") cluster",cd8_clusters[i],"center (",cd8_centers[i,1],",",cd8_centers[i,2],") in",file_names[images[i]]))
	cd8_colors[i]=m[cd8_clusters[i],9]
	cd8_maxx[i]=m[cd8_clusters[i],10]
	cd8_maxy[i]=m[cd8_clusters[i],11]
    }

    if (images[i] == curr_image)
    {
	objs = objs+1
    }
    else
    {
	x<-cbind(cd8_xs[pobjs:objs],cd8_ys[pobjs:objs])
	colnames(x) <- c("x", "y")
	tot_objs = objs - pobjs + 1
	pdf(paste(file_names[curr_image],"new_CD8_clustering.pdf",sep="_"))
	print(paste("Plotting CD8s from file",file_names[curr_image]))
	print(paste("First index",pobjs,"at (",cd8_xs[pobjs],",",cd8_ys[pobjs],"), last index",objs,"at (",cd8_xs[objs],",",cd8_ys[objs],")"))
	plot(x, xlim=c(0,cd8_maxx[pobjs]),ylim=c(0,cd8_maxy[pobjs]),col = cd8_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD8s from",rep_clusters,"clusters)"))
	points(x=cd8_centers[pobjs:objs,1],y=cd8_centers[pobjs:objs,2], col=cd8_colors[pobjs:objs],pch = 8, cex = 2)
	dev.off()
	#pobjs=objs+1
	pobjs=i
	objs = i
	curr_image=images[i]
	rep_clusters=0
    }
}
x<-cbind(cd8_xs[pobjs:objs],cd8_ys[pobjs:objs])
colnames(x) <- c("x", "y")
tot_objs = objs - pobjs + 1
print(paste("Plotting CD8s from file",file_names[curr_image]))
print(paste("First index",pobjs,"at (",cd8_xs[pobjs],",",cd8_ys[pobjs],"), last index",objs,"at (",cd8_xs[objs],",",cd8_ys[objs],")"))
pdf(paste(file_names[curr_image],"new_CD8_clustering.pdf",sep="_"))
plot(x, xlim=c(0,cd8_maxx[pobjs]),ylim=c(0,cd8_maxy[pobjs]),col = cd8_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD8s from",rep_clusters,"clusters)"))
points(x=cd8_centers[pobjs:objs,1],y=cd8_centers[pobjs:objs,2], col=cd8_colors[pobjs:objs],pch = 8, cex = 2)
dev.off()

write("Image,cluster,CD4s,CD8s,Cells,E/T ratio (cd8s),E/T ratio (cd4s),E/T ratio (combo)", file="new_cluster_data.csv",append=FALSE)

clusts=0
pclusts=1
curr_image=1
for (j in 1:tot_clusters)
{
    et_ratio_cd8=1
    et_ratio_cd4=1
    et_ratio_both=1
    if (m[j,6] >= m[j,8])
	et_ratio_cd8=(m[j,8])/(m[j,6]);
    if (m[j,6] >= m[j,7])
	et_ratio_cd4=(m[j,7])/(m[j,6]);
    if (m[j,6] >= m[j,7]+m[j,8])
	et_ratio_both=(m[j,7]+m[j,8])/(m[j,6]);
    write(paste(m[j,1],m[j,2],m[j,7],m[j,8],m[j,6],et_ratio_cd8,et_ratio_cd4,et_ratio_both,sep=","),file="new_cluster_data.csv",append=TRUE)
    if (m[j,1] == curr_image)
    {
	clusts = clusts+1
    }
    else
    {
	effs<-(m[pclusts:clusts,7]+m[pclusts:clusts,8])
	y<-m[pclusts:clusts,8]/(m[pclusts:clusts,6]-effs)
	y_cd8<-sort(y,decreasing=TRUE)
	y<-m[pclusts:clusts,7]/(m[pclusts:clusts,6]-effs)
	y_cd4<-sort(y,decreasing=TRUE)
	y<-effs/(m[pclusts:clusts,6]-effs)
	y_combo<-sort(y,decreasing=TRUE)
	header="1"
	line1=y_cd8[1]
	line2=y_cd4[1]
	line3=y_combo[1]
	for (k in 2:length(y_cd8))
	{
	    header=paste(header,k,sep=",")
	    if (y_cd8[k] != 0)
	    {
		line1=paste(line1,y_cd8[k],sep=",")
	    } else {
		line1=paste(line1,"0.001",sep=",")
		y_cd8[k]=0.001
	    }
	    if (y_cd4[k] != 0)
	    {
		line2=paste(line2,y_cd4[k],sep=",")
	    } else {
		line2=paste(line2,"0.001",sep=",")
		y_cd4[k]=0.001
	    }
	    if (y_combo[k] != 0)
	    {
		line3=paste(line3,y_combo[k],sep=",")
	    } else {
		line3=paste(line3,"0.001",sep=",")
		y_combo[k]=0.001
	    }
	}
	write(header, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=FALSE)
	write(line1, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)
	write(line2, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)
	write(line3, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)
	pclusts=clusts+1

	ranks<-seq(1,length(y_combo))
	exponential.model <- lm(log(y_combo)~ ranks)
	exp_a<-exp(coef(exponential.model)[1])
	p_combo_coefs[curr_image]=exp_a
	exp_b<-coef(exponential.model)[2]
	c_combo_coefs[curr_image]=exp_b
	r2_combo[curr_image]<-summary(exponential.model)$r.squared
	equat<-paste(paste("y=",format(exp_a,digits=2),"*e^",format(exp_b,digits=2),"*x",sep=""),
	    paste("R^2=",format(r2_combo[curr_image],2),sep=""),sep="\n")
	exp_combo <- exp(predict(exponential.model,list(ranks=ranks)))

	ranks<-seq(1,length(y_cd8))
	exponential.model <- lm(log(y_cd8)~ ranks)
	exp_a<-exp(coef(exponential.model)[1])
	p_cd8_coefs[curr_image]=exp_a
	exp_b<-coef(exponential.model)[2]
	c_cd8_coefs[curr_image]=exp_b
	rsquared<-summary(exponential.model)$r.squared
	cd8_equat<-paste(paste("y=",format(exp_a,digits=2),"*e^",format(exp_b,digits=2),"*x",sep=""),
	    paste("R^2=",format(rsquared,2),sep=""),sep="\n")
	exp_cd8 <- exp(predict(exponential.model,list(ranks=ranks)))

	if (do_plots == 1){
	png(paste(file_names[curr_image],"et.png",sep="_"))
	plot(x=ranks, y=y_combo,xlab="Regional Ranking",ylab="E/T ratio",main=paste(file_names[curr_image],"E/T ratios\n(cluster size =",cluster_size,")"),col=c("blue"),type="l",lty=1)
	lines(x=ranks,y=y_cd8,col=c("red"),lty=2)
	lines(x=ranks,y=y_cd4,col=c("green"),lty=3)
	lines(x=ranks,y=exp_combo,col=c("black"),lty=4)
	legend(x=max(ranks)/2,y=max(y_cd8,y_cd4,y_combo),legend=c("CD4s+CD8s","CD8s","CD4s","Expon Fit"),col=c("blue","red","green","black"),lty=1:4,cex=0.8)
	text(x=max(ranks)/2,y=max(y_cd8,y_cd4,y_combo)/2,equat)
	dev.off()
	}

	curr_image=m[j,1]
    }
}
effs<-(m[pclusts:clusts,7]+m[pclusts:clusts,8])
y<-m[pclusts:clusts,8]/(m[pclusts:clusts,6]-effs)
y_cd8<-sort(y,decreasing=TRUE)
y<-m[pclusts:clusts,7]/(m[pclusts:clusts,6]-effs)
y_cd4<-sort(y,decreasing=TRUE)

effs<-(m[pclusts:clusts,7]+m[pclusts:clusts,8])
y<-effs/(m[pclusts:clusts,6]-effs)
y_combo<-sort(y,decreasing=TRUE)
header="1"
line1=y_cd8[1]
line2=y_cd4[1]
line3=y_combo[1]
for (k in 2:length(y_cd8))
{
    header=paste(header,k,sep=",")
    if (y_cd8[k] != 0)
    {
	line1=paste(line1,y_cd8[k],sep=",")
    } else {
	line1=paste(line1,"0.001",sep=",")
	y_cd8[k]=0.001
    }
    if (y_cd4[k] != 0)
    {
	line2=paste(line2,y_cd4[k],sep=",")
    } else {
	line2=paste(line2,"0.001",sep=",")
	y_cd4[k]=0.001
    }
    if (y_combo[k] != 0)
    {
	line3=paste(line3,y_combo[k],sep=",")
    } else {
	line3=paste(line3,"0.001",sep=",")
	y_combo[k]=0.001
    }
}
write(header, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=FALSE)
write(line1, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)
write(line2, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)
write(line3, file=paste(file_names[curr_image],"_new_ranked_ratios.csv",sep=""),append=TRUE)

#determine fit for last plot and add to png
ranks<-seq(1,length(y_cd8))
exponential.model <- lm(log(y_cd8)~ ranks)
exp_a<-exp(coef(exponential.model)[1])
p_cd8_coefs[curr_image]=exp_a
exp_b<-coef(exponential.model)[2]
c_cd8_coefs[curr_image]=exp_b
cd8_equat<-paste(paste("y=",format(exp_a,digits=2),"*e^",format(exp_b,digits=2),"*x",sep=""),
    paste("R^2=",format(rsquared,2),sep=""),sep="\n")
exp_cd8 <- exp(predict(exponential.model,list(ranks=ranks)))

ranks<-seq(1,length(y_combo))
exponential.model <- lm(log(y_combo)~ ranks)
exp_a<-exp(coef(exponential.model)[1])
p_combo_coefs[curr_image]=exp_a
exp_b<-coef(exponential.model)[2]
c_combo_coefs[curr_image]=exp_b
r2_combo[curr_image]<-summary(exponential.model)$r.squared
equat<-paste(paste("y=",format(exp_a,digits=2),"*e^",format(exp_b,digits=2),"*x",sep=""),
    paste("R^2=",format(r2_combo[curr_image],2),sep=""),sep="\n")
exp_combo <- exp(predict(exponential.model,list(ranks=ranks)))

if (do_plots == 1){
png(paste(file_names[curr_image],"et.png",sep="_"))
plot(x=ranks, y=y_combo,xlab="Regional Ranking",ylab="E/T ratio",main=paste(file_names[curr_image],"E/T ratios\n(cluster size =",cluster_size,")"),col=c("blue"),type="l",lty=1)
lines(x=ranks,y=y_cd8,col=c("red"),lty=2)
lines(x=ranks,y=y_cd4,col=c("green"),lty=3)
lines(x=ranks,y=exp_combo,col=c("black"),lty=4)
legend(x=max(ranks)/2,y=max(y_cd8,y_cd4,y_combo),legend=c("CD4s+CD8s","CD8s","CD4s","Expon Fit"),col=c("blue","red","green","black"),lty=1:4,cex=0.8)
text(x=max(ranks)/2,y=max(y_cd8,y_cd4,y_combo)/2,equat)
dev.off()
}

# file of ALL E/T ratios across all files
effs<-(m[1:tot_clusters,8])
y<-effs/(m[1:tot_clusters,6]-effs)
y_cd8<-sort(y,decreasing=TRUE)

effs<-(m[1:tot_clusters,7]+m[1:tot_clusters,8])
y<-effs/(m[1:tot_clusters,6]-effs)
y_combo<-sort(y,decreasing=TRUE)
header="1"
line1=y_cd8[1]
line2=y_combo[1]
for (k in 2:length(y_combo))
{
    header=paste(header,k,sep=",")
    if (y_cd8[k] != 0)
    {
	line1=paste(line1,y_cd8[k],sep=",")
    } else {
	line1=paste(line1,"0.001",sep=",")
    }
    if (y_combo[k] != 0)
    {
	line2=paste(line2,y_combo[k],sep=",")
    } else {
	line2=paste(line2,"0.001",sep=",")
    }
}
write(header, file="all_new_ranked_ratios.csv",append=FALSE)
write(line1, file="all_new_ranked_ratios.csv",append=TRUE)
write(line2, file="all_new_ranked_ratios.csv",append=TRUE)

#determine fit for CD8 only E/T ratio coefs (poly)
ranks<-seq(1,length(file_names))
sorted_p_cd8_coefs<-sort(p_cd8_coefs,decreasing=TRUE)
poly.model <- lm(sorted_p_cd8_coefs~ poly(ranks,3))
summary(poly.model)
rsquared<-summary(poly.model)$r.squared
coef_a<-poly.model$coefficient[4]
coef_b<-poly.model$coefficient[3]
coef_c<-poly.model$coefficient[2]
coef_d<-poly.model$coefficient[1]
equat<-paste(paste("y=",format(coef_a,digits=2),"rank^3+",format(coef_b,digits=2),"*rank^2+",format(coef_c,digits=2),"*rank+",format(coef_d,digits=2),sep=""),
    paste("R^2=",format(rsquared,2),sep=""),sep="\n")
fit_cd8s <- predict(poly.model,list(ranks=ranks))
paste("P coef cd8 equat:",equat)
sorted_p_cd8_coefs

if (do_plots == 1){
png("p_cd8_poly.png")
plot(x=ranks, y=sorted_p_cd8_coefs,xlab="Sample Ranking",ylab="P coef",main="Rank Order P Coef",col=c("blue"),type="l",lty=1)
lines(x=ranks,y=fit_cd8s,col=c("red"),lty=3)
legend(x=max(ranks)/2,y=max(sorted_p_cd8_coefs,fit_cd8s),legend=c("P coefs (CD8s)","Poly Fit"),col=c("blue","red"),lty=1:2,cex=0.8)
text(x=max(ranks)/2,y=max(p_cd8_coefs,fit_cd8s)/2,equat)
dev.off()
}
c_cd8_coefs
mean(c_cd8_coefs)
sd(c_cd8_coefs)

#determine fit for CD4+CD8 E/T ratio coefs
ranks<-seq(1,length(file_names))
sorted_p_combo_coefs<-sort(p_combo_coefs,decreasing=TRUE)
poly.model <- lm(sorted_p_combo_coefs~ poly(ranks,3))
summary(poly.model)
rsquared<-summary(poly.model)$r.squared
coef_a<-poly.model$coefficient[4]
coef_b<-poly.model$coefficient[3]
coef_c<-poly.model$coefficient[2]
coef_d<-poly.model$coefficient[1]
equat<-paste(paste("y=",format(coef_a,digits=2),"rank^3+",format(coef_b,digits=2),"*rank^2+",format(coef_c,digits=2),"*rank+",format(coef_d,digits=2),sep=""),
    paste("R^2=",format(rsquared,2),sep=""),sep="\n")
fit_combo <- predict(poly.model,list(ranks=ranks))
paste("P coef equat:",equat)
sorted_p_combo_coefs

if (do_plots == 1){
png("p_combo_poly.png")
plot(x=ranks, y=sorted_p_combo_coefs,xlab="Sample Ranking",ylab="P coef",main="Rank Order P Coef",col=c("blue"),type="l",lty=1)
lines(x=ranks,y=fit_combo,col=c("red"),lty=3)
legend(x=max(ranks)/2,y=max(sorted_p_combo_coefs,fit_combo),legend=c("R.O. P coefs","Poly Fit"),col=c("blue","red"),lty=1:2,cex=0.8)
text(x=max(ranks)/2,y=max(sorted_p_combo_coefs,fit_combo)/2,equat)
dev.off()
}

ranks<-seq(1,length(file_names))
p_combo_coefs<-sort(p_combo_coefs,decreasing=TRUE)
exponential.model <- lm(log(p_combo_coefs)~ ranks)
rsquared<-summary(exponential.model)$r.squared
exp_a<-exp(coef(exponential.model)[1])
exp_b<-coef(exponential.model)[2]
equat<-paste(paste("y=",format(exp_a,digits=2),"*e^",format(exp_b,digits=2),"*x",sep=""),
    paste("R^2=",format(rsquared,2),sep=""),sep="\n")
exp_combo <- exp(predict(exponential.model,list(ranks=ranks)))
paste("P coef exp equat:",equat)

if (do_plots == 1){
png("p_combo_exp.png")
plot(x=ranks, y=p_combo_coefs,xlab="Sample Ranking",ylab="P coef",main="Rank Order P Coef",col=c("blue"),type="l",lty=1)
lines(x=ranks,y=exp_combo,col=c("red"),lty=3)
legend(x=max(ranks)/2,y=max(p_combo_coefs,exp_combo),legend=c("R.O. P coefs","Expon Fit"),col=c("blue","red"),lty=1:2,cex=0.8)
text(x=max(ranks)/2,y=max(p_combo_coefs,exp_combo)/2,equat)
dev.off()
}
c_combo_coefs
mean(c_combo_coefs)
sd(c_combo_coefs)
r2_combo
median(r2_combo)
min(r2_combo)
max(r2_combo)
