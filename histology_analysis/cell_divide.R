#Analysis of histology images to cluster cells and CD8s
#Dave Swan, Jul 2018
#Edits Pavitra Roychoudhury
#Citation info: <tbd>

rm(list=ls())
library(ggplot2)
library(gplots)

data <- read.table("base_names.txt", header=FALSE,comment.char="#")
file_names<-data[,1];
file_names

#Columns are:
#ImageNumber,ObjectNumber,AreaShape_Area,AreaShape_Center_X,AreaShape_Center_Y,AreaShape_Compactness,AreaShape_Eccentricity, 
#AreaShape_EulerNumber,AreaShape_Extent,AreaShape_FormFactor,AreaShape_MajorAxisLength,AreaShape_MaxFeretDiameter,
#AreaShape_MaximumRadius,AreaShape_MeanRadius,AreaShape_MedianRadius,AreaShape_MinFeretDiameter,AreaShape_MinorAxisLength,
#AreaShape_Orientation,AreaShape_Perimeter,AreaShape_Solidity,Location_Center_X,Location_Center_Y,Number_Object_Number
data <- read.csv("rough_cells.csv", comment.char="#");
names(data);


#Nuclei clustering
nimages<-length(unique(data[,1]));
for(i in 1:nimages){
	tempdata<-subset(data,ImageNumber==i);
	df<-data.frame(xs=tempdata$AreaShape_Center_X,ys=tempdata$AreaShape_Center_Y)
	clusters<-floor(nrow(tempdata)/100)

	cl <- kmeans(df, clusters);
	
	clusterinfo<-data.frame(
		cluster=c(1:nrow(cl$centers)),
		cluster_x=cl$centers[,'xs'],
		cluster_y=cl$centers[,'ys'],
		cluster_size=cl$size,
		img_id=file_names[i],
		stringsAsFactors=F
	)
	write.csv(clusterinfo,paste0(file_names[i],'_clusterinfo.csv'),row.names=F)
	
	print(paste("Plotting ",nrow(tempdata),"cells and",clusters,"clusters from file",as.character(file_names[i])))
	pdf(paste(file_names[i],"nuclei_clustering_pavi.pdf",sep="_"))
	# svg(paste(file_names[i],"nuclei_clustering_pavi.svg",sep="_"),width=7,height=4)
	plot(df, col = cl$cluster, main=paste(file_names[i],"(",clusters,"Clusters )"))
	points(cl$centers, col=1:clusters,pch = 8, cex = 2)
	dev.off();
	
	tempdata$cluster<-cl$cluster;
	write.csv(tempdata,paste0(file_names[i],'_cells.csv'),row.names=F)
}

#Assign CD8s to a cluster
data2 <- read.csv("rough_CD8s.csv", comment.char="#");
head(data2)
nimages<-length(unique(data[,1]));

for(i in 1:nimages){
	cd8data<-subset(data2,ImageNumber==i);
	cd8data$nearest_cell_cluster<='';
	clusterdata<-read.csv(paste0(file_names[i],'_clusterinfo.csv'));
	for(j in 1: nrow(cd8data)){
		cd8data$nearest_cell_cluster[j]<-clusterdata$cluster[which.min(sqrt((clusterdata$cluster_x-cd8data$AreaShape_Center_X[j])^2+
																																					(clusterdata$cluster_y-cd8data$AreaShape_Center_Y[j])^2))]
	}
	pdf(paste(file_names[i],"CD8_clustering_pavi.pdf",sep="_"))
	# svg(paste(file_names[i],"CD8_clustering.svg",sep="_"),width=7,height=4)
	print(paste("Plotting CD8s from file",file_names[i]))
	plot(cbind(cd8data$AreaShape_Center_X,cd8data$AreaShape_Center_Y), col = 1:nrow(clusterdata), xlab='x',ylab='y',
			 main=paste(file_names[i],"(",nrow(cd8data),"CD8s from",length(unique(cd8data$nearest_cell_cluster)),"clusters)"))
	plot_data<-subset(clusterdata,cluster %in% cd8data$nearest_cell_cluster)
	points(cbind(plot_data$cluster_x,plot_data$cluster_y),col = 1:nrow(clusterdata), pch = 8, cex = 2)
	dev.off()
	
	#count cd8s in each cluster and add to file
	clusterdata$cd8_count<-unlist(lapply(clusterdata$cluster,function(x)sum(cd8data$nearest_cell_cluster==x)));
	clusterdata$et_ratio<-clusterdata$cd8_count/(clusterdata$cluster_size-clusterdata$cd8_count);
	write.csv(clusterdata,paste0(file_names[i],'_clusterinfo.csv'),row.names=F)
	write.csv(cd8data,paste0(file_names[i],'_cd8s.csv'),row.names=F)
}


#For the paper Fig 2B and C (14095 2wph)
#Here I'm trying to make the colours consistent across the cell and cd8 panel
library(RColorBrewer);
library(ggplot2)
n<-70; #need a nice big pallete to accommodate all the clusters
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)));
colours<-data.frame(cluster=1:70,colour=col_vector[1:70],stringsAsFactors=F);
# pie(rep(1,n), col=sample(col_vector, n))

ind<-which(file_names=='14095_2wph');
clusters_cells<-read.csv(paste0(file_names[ind],'_clusterinfo.csv'),stringsAsFactors=F);
celldata<-read.csv(paste0(file_names[ind],'_cells.csv'),stringsAsFactors=F);
cd8data<-read.csv(paste0(file_names[ind],'_cd8s.csv'),stringsAsFactors=F);
celldata$type<-'cell';
cd8data$type<-'CD8'
clusters_cd8s<-subset(clusters_cells,cluster%in%unique(cd8data$nearest_cell_cluster));
clusters_cells$type<-'cell';
clusters_cd8s$type<-'CD8'
names(cd8data)[names(cd8data)=='nearest_cell_cluster']<-'cluster';
plot_data<-rbind(celldata,cd8data);
cluster_data<-rbind(clusters_cells,clusters_cd8s);
rm(celldata,cd8data,clusters_cells,clusters_cd8s);

# pdf('14095_2wph_clustering_pavi_ggplot.pdf',width=7,height=4);
svg('14095_2wph_clustering_pavi_ggplot.svg',width=7,height=4);
ggplot(data=plot_data,aes(x=AreaShape_Center_X,y=AreaShape_Center_Y))+
	geom_point(aes(colour=factor(cluster,levels=1:70),shape=type,size=type))+
	scale_color_manual(breaks=colours$cluster,values=colours$colour,
										 guide=FALSE)+
  scale_shape_manual(values=c(17,1),guide=FALSE)+
  scale_size_manual(values=c(2,1),guide=FALSE)+
	geom_point(data=cluster_data,aes(x=cluster_x,y=cluster_y),pch='+',size=3)+
	theme_void()+
  facet_grid(type~.)
dev.off()

#Now plot E:T ratios
#impoort files, compute rank order

rank_order<-function(filename){
  temp<-read.csv(paste0(filename,'_clusterinfo.csv'),stringsAsFactors=F);
  temp$rank<-rank(-temp$et_ratio)
  return(temp)
}
plot_data<-do.call(rbind,lapply(file_names,rank_order))
head(plot_data);
# write.csv(plot_data,'./All_images_clusterinfo.csv',row.names=F)

#TBD: ggplot
# ggplot(plot_data)+
#   geom_histogram(aes(x=et_ratio))+
#   facet_wrap(~img_id)+
#   theme_bw
# pdf('ETratios_all.pdf')
svg('ETratios_all.svg',width=3,height=2)
ggplot(plot_data,aes(x=rank,y=et_ratio,group=img_id))+
  geom_line(aes(colour=img_id))+
  theme_bw()+
  ylab('E:T ratio')+
  theme(legend.position='none')
dev.off()
pdf('ETratios_all_log.pdf')
ggplot(subset(plot_data,et_ratio>0),aes(x=rank,y=et_ratio,group=img_id))+
  geom_line(aes(colour=img_id))+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position='none')
dev.off()

plot_data<-read.csv('All_images_clusterinfo.csv',stringsAsFactors=F)
# plot_data2<-
#   plot_data<-plot_data%>%
#   group_by(img_id)
plot_data$wk<-''
plot_data$wk[grepl('2wph',plot_data$img_id)]<-'2 wph'
plot_data$wk[grepl('8wph',plot_data$img_id)]<-'8 wph';
plot_data$img_num<-unlist(lapply(plot_data$img_id,function(x)
  strsplit(x,'_')[[1]][1]))
write.csv(plot_data,'All_images_clusterinfo.csv',row.names=F)

plot_data2<-plot_data;
plot_data2$wk<-'all'


svg('ETratios_splitweeks_plusall.svg',width=6,height=2)
ggplot(rbind(plot_data,plot_data2),aes(x=rank,y=et_ratio,group=img_id))+
  geom_line(aes(colour=img_id))+
  theme_bw()+
  ylab('E:T ratio')+
  theme(legend.position='none')+
  facet_grid(~wk)
dev.off()


## Fig 2D
#not sure if "all" is needed, josh wants log scale 
plot_data3<-plot_data;
plot_data3$et_ratio[plot_data3$et_ratio==0]<-0.005; #for visual display
svg('ETratios_splitweeks.svg',width=3,height=3.5)
# pdf('ETratios_splitweeks.pdf',width=4,height=2)
ggplot(plot_data3,aes(x=rank,y=et_ratio,group=img_num))+
  geom_line(aes(colour=img_num),size=1)+
	scale_y_log10(limits=c(min(plot_data3$et_ratio),2))+
  theme_bw()+
  ylab('E:T ratio')+
  theme(legend.position='none')+
  facet_grid(wk~.)+
	theme(strip.background = element_blank(),
				strip.text = element_blank())+
	geom_text(aes(label = wk), x = Inf, y = Inf, hjust = 1.5, vjust = 1.5) 
dev.off()


## Fig 2E
#just the 0s
#"To highlight that there are lots of zero regions, let's do a scatter plot with total E:T ratio 
# of the whole biopsy on the x-axis and proportion of regions with E:T=0 on the y-axis.  Presumably, 
# there will be some correlation. I'd do separate scatter plots for 2 & 8 weeks."
library(dplyr)
plot_data4<-group_by(plot_data,img_id,wk,img_num)%>%
	summarize(clusters=n(),
						et_ratio_all=sum(cd8_count)/(sum(cluster_size-cd8_count)),
						zeros_prop=sum(cd8_count==0)/clusters)
svg('ETratios_zeros.svg',width=3,height=3.5)
# pdf('ETratios_zeros.pdf',width=4,height=2)
ggplot(plot_data4,aes(x=et_ratio_all,y=zeros_prop))+
	geom_point()+
	theme_bw()+
	ylab('Prop regions with E:T=0')+
	xlab('E:T ratio of the whole biopsy')+
	theme(legend.position='none',
				axis.text.x=element_text(angle=45,hjust=1))+
	facet_grid(wk~.)+
	geom_text(aes(label = wk), x = Inf, y = Inf, hjust = 1.5, vjust = 1.5) +
	theme(strip.background = element_blank(),
				strip.text = element_blank())
# +
# 	geom_smooth(method = 'lm')
# + geom_smooth(method = 'lm', formula = y ~ poly(x,2), aes(colour = 'polynomial'), se= FALSE) +
# 	geom_smooth(method = 'nls', formula = y ~ a * log(x) +b, aes(colour = 'logarithmic'), se = FALSE, start = list(a=1,b=1)) +
# 	geom_smooth(method = 'nls', formula = y ~ a*exp(b *x), aes(colour = 'Exponential'), se = FALSE, start = list(a=1,b=1))
dev.off()


## Return distributions  used for fitting
library(broom)
library(dplyr)
plot_data<-read.csv('All_images_clusterinfo.csv',stringsAsFactors=F)
nlsfit<-do.call('rbind',lapply(file_names, function(f){
  tmp<-tidy(nls(et_ratio~ a*(exp(b*rank)),subset(plot_data,img_id==f),start=list(a=1,b=0)));
  tmp$img_id<-f;
  return(tmp)
}));
write.csv(nlsfit,'./Exp_fit_ETratios_all.csv',row.names=F)
#for statistic quoted in text: The exponential slopes of the rank order abundance
#curves (median: __, IQR: __, range: _) were strikingly similar across all specimens
median(subset(nlsfit,term=='b')$estimate)
IQR(subset(nlsfit,term=='b')$estimate)
range(subset(nlsfit,term=='b')$estimate)

## For the statistic quoted in text: 
# "noted ratios of TRM to epidermal cells ranging from 0 to 1 with a range of medians across all 20 biopsies of 0.01 to 0.24 (IQR: ). "
library(dplyr)
length(file_names); #there are actually 19 biopsies, not 20
stat1<-
	plot_data%>%
	group_by(img_id,wk,img_num)%>%
	summarise(n_clusters=n(),
						median_et_ratio=median(et_ratio))
IQR(stat1$median_et_ratio)
summary(stat1$median_et_ratio)
summary(plot_data$et_ratio)
nrow(plot_data)
sum(plot_data$et_ratio==0)
sum(plot_data$et_ratio==0)/nrow(plot_data)*100




#For a sentence in text about ET ratios across all biopsies 
library(dplyr)
all_data<-read.csv('All_images_clusterinfo.csv',stringsAsFactors=F)
head(all_data)

summ1<-all_data%>%
	group_by(img_id)%>%
	summarize(median_ET_allclusters=median(et_ratio),
						iqr_ET_allclusters=IQR(et_ratio),
						min_ET_all_clusters=min(et_ratio),
						max_ET_all_clusters=max(et_ratio),
						ET_whole_biopsy=sum(cd8_count)/(sum(cluster_size-cd8_count)));
write.csv(summ1,'./Summary_ETratios.csv',row.names=F)
median(summ1$ET_whole_biopsy)	
IQR(summ1$ET_whole_biopsy)
quantile(summ1$ET_whole_biopsy)
range(summ1$ET_whole_biopsy)
nrow(summ1)
