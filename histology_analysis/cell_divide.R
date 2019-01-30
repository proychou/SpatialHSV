library(ggplot2)
library(gplots)

data <- read.table("base_names.txt", header=FALSE,comment.char="#")
file_names<-data[,1]

data <- read.csv("rough_cells.csv", comment.char="#")
#ImageNumber,ObjectNumber,AreaShape_Area,AreaShape_Center_X,AreaShape_Center_Y,AreaShape_Compactness,AreaShape_Eccentricity,AreaShape_EulerNumber,AreaShape_Extent,AreaShape_FormFactor,AreaShape_MajorAxisLength,AreaShape_MaxFeretDiameter,AreaShape_MaximumRadius,AreaShape_MeanRadius,AreaShape_MedianRadius,AreaShape_MinFeretDiameter,AreaShape_MinorAxisLength,AreaShape_Orientation,AreaShape_Perimeter,AreaShape_Solidity,Location_Center_X,Location_Center_Y,Number_Object_Number

N<-nrow(data)
images=as.numeric(data[,1])
xs=as.numeric(data[,4])
ys=as.numeric(data[,5])
curr_image=images[1]
m<-matrix(data=NA,nrow=N,ncol=8)
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
	clusters=as.integer(tot_objs/100)
	(cl <- kmeans(x, clusters))
	for (j in 1:clusters)
	{
	    m[tot_clusters+j,1]=curr_image
	    m[tot_clusters+j,2]=j
	    m[tot_clusters+j,3]=cl$centers[j,1]
	    m[tot_clusters+j,4]=cl$centers[j,2]
	    m[tot_clusters+j,5]=cl$cluster[j]
	    m[tot_clusters+j,6]=0
	    m[tot_clusters+j,7]=cl$size[j]
	    #m[tot_clusters+j,8]=0
	    m[tot_clusters+j,8]=j
	}
	tot_clusters=tot_clusters+j

	print(paste("Plotting ",tot_objs,"cells and",clusters,"clusters from file",file_names[curr_image]))
	pdf(paste(file_names[curr_image],"nuclei_clustering.pdf",sep="_"))
	plot(x, col = cl$cluster, main=paste(file_names[curr_image],"(",clusters,"Clusters )"))
	points(cl$centers, col=1:clusters,pch = 8, cex = 2)
	dev.off()
	pobjs=objs+1
	curr_image=images[i]
    }
}
x<-cbind(xs[pobjs:objs],ys[pobjs:objs])
colnames(x) <- c("x", "y")
tot_objs = objs - pobjs + 1
clusters=as.integer(tot_objs/100)
(cl <- kmeans(x, clusters))
for (j in 1:clusters)
{
    m[tot_clusters+j,1]=curr_image
    m[tot_clusters+j,2]=j
    m[tot_clusters+j,3]=cl$centers[j,1]
    m[tot_clusters+j,4]=cl$centers[j,2]
    m[tot_clusters+j,5]=cl$cluster[j]
    m[tot_clusters+j,6]=0
    m[tot_clusters+j,7]=cl$size[j]
    #m[tot_clusters+j,8]=0
    m[tot_clusters+j,8]=cl$cluster[j]
}
tot_clusters=tot_clusters+j
print(paste("Plotting ",tot_objs,"cells and",clusters,"clusters from file",file_names[curr_image]))
print(paste("First cluster is",cl$cluster[1]))
print(paste("First nuc cluster is",cl$cluster[1],"at",xs[pobjs],",",ys[pobjs]))
pdf(paste(file_names[curr_image],"nuclei_clustering.pdf",sep="_"))
plot(x, col = cl$cluster, main=paste(file_names[curr_image],"(",clusters,"Clusters )"))
points(cl$centers, col=1:clusters,pch = 8, cex = 2)
m<-m[1:tot_clusters,]
dev.off()

data2 <- read.csv("rough_CD8s.csv", comment.char="#")
CD8s<-nrow(data2)
images=as.numeric(data2[,1])
cd8_obj_num=as.numeric(data2[,2])
cd8_xs=as.numeric(data2[,4])
cd8_ys=as.numeric(data2[,5])
cd8_clusters<-vector(length=CD8s)
cd8_colors<-vector(length=CD8s)
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
	if(m[cd8_clusters[i],6]==0)
	{
	    rep_clusters=rep_clusters+1
	    #m[cd8_clusters[i],8]= rep_clusters
	}
	m[cd8_clusters[i],6]= m[cd8_clusters[i],6]+1
	cd8_centers[i,1]=as.integer(m[cd8_clusters[i],3])
	cd8_centers[i,2]=as.integer(m[cd8_clusters[i],4])
	print(paste("For object",cd8_obj_num[i],"at (",cd8_xs[i],",",cd8_ys[i],") cluster",cd8_clusters[i],"center (",cd8_centers[i,1],",",cd8_centers[i,2],") in",file_names[images[i]]))
	cd8_colors[i]=m[cd8_clusters[i],8]
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
	pdf(paste(file_names[curr_image],"CD8_clustering.pdf",sep="_"))
	print(paste("Plotting CD8s from file",file_names[curr_image]))
	print(paste("First CD8 cluster is",cd8_colors[pobjs],"at",cd8_xs[pobjs],",",cd8_ys[pobjs]))
	plot(x, col = cd8_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD8s from",rep_clusters,"clusters)"))
	points(x=cd8_centers[pobjs:objs,1],y=cd8_centers[pobjs:objs,2], col=cd8_colors[pobjs:objs],pch = 8, cex = 2)
	dev.off()
	pobjs=objs+1
	curr_image=images[i]
	rep_clusters=0
    }
}
x<-cbind(cd8_xs[pobjs:objs],cd8_ys[pobjs:objs])
colnames(x) <- c("x", "y")
tot_objs = objs - pobjs + 1
print(paste("Plotting CD8s from file",file_names[curr_image]))
pdf(paste(file_names[curr_image],"CD8_clustering.pdf",sep="_"))
plot(x, col = cd8_colors[pobjs:objs], main=paste(file_names[curr_image],"(",tot_objs,"CD8s from",rep_clusters,"clusters)"))
points(x=cd8_centers[pobjs:objs,1],y=cd8_centers[pobjs:objs,2], col=cd8_colors[pobjs:objs],pch = 8, cex = 2)
dev.off()

write("Image,cluster,CD8s,Cells,E/T ratio", file="cluster_data.csv",append=FALSE)

clusts=0
pclusts=1
curr_image=1
for (j in 1:tot_clusters)
{
    et_ratio=1
    if (m[j,7] != m[j,6])
	et_ratio=m[j,6]/(m[j,7]-m[j,6]);
    write(paste(m[j,1],m[j,2],m[j,6],m[j,7],et_ratio,sep=","),file="cluster_data.csv",append=TRUE)
    if (m[j,1] == curr_image)
    {
	clusts = clusts+1
    }
    else
    {
	y<-m[pclusts:clusts,6]/(m[pclusts:clusts,7]-m[pclusts:clusts,6])
	y2<-sort(y,decreasing=TRUE)
	header="1"
	line=y2[1]
	for (k in 2:length(y2))
	{
	    if (y2[k] != 0)
	    {
		header=paste(header,k,sep=",")
		line=paste(line,y2[k],sep=",")
	    }
	}
	write(header, file=paste(file_names[curr_image],"_ranked_ratios.csv",sep=""),append=FALSE)
	write(line, file=paste(file_names[curr_image],"_ranked_ratios.csv",sep=""),append=TRUE)
	pclusts=clusts+1
	curr_image=m[j,1]
    }
}
y<-m[pclusts:clusts,6]/(m[pclusts:clusts,7]-m[pclusts:clusts,6])
y2<-sort(y,decreasing=TRUE)
header="1"
line=y2[1]
for (k in 2:length(y2))
{
    if (y2[k] != 0)
    {
	header=paste(header,k,sep=",")
	line=paste(line,y2[k],sep=",")
    }
}
write(header, file=paste(file_names[curr_image],"_ranked_ratios.csv",sep=""),append=FALSE)
write(line, file=paste(file_names[curr_image],"_ranked_ratios.csv",sep=""),append=TRUE)

y<-m[,6]/(m[,7]-m[,6])
y2<-sort(y,decreasing=TRUE)
header="1"
line=y2[1]
for (k in 2:length(y2))
{
    if (y2[k] != 0)
    {
	header=paste(header,k,sep=",")
	line=paste(line,y2[k],sep=",")
    }
}
write(header, file="all_ranked_ratios.csv",append=FALSE)
write(line, file="all_ranked_ratios.csv",append=TRUE)
