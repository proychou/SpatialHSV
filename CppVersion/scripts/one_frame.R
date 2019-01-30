
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
cyt_ic50<-7.406*1.25e-4; #default
args<-commandArgs(trailingOnly=T);

run<-args[2];
time<-args[3];
plot_cytokines<-args[4];
cyt_expiration<-as.numeric(args[5]);
cyt_ic50<-as.numeric(args[6]);
cyt_hill<-as.numeric(args[7]);
max_log_vl<-as.numeric(args[8]);
max_time<-as.numeric(args[9]);
panels<-as.numeric(args[10]);
cyto_act<-as.numeric(args[11]);

if (panels == 4) {
    png(paste(args[1],"snapshot_",run,"_",time,".png",sep=""),width = 6, height = 6, units="in",res=200);
    par(mar=c(1,1,2,1),oma=c(1,0,0,0),cex=0.5);
    layout(matrix(c(1:4),2,2,byrow = TRUE),widths=c(1,1,1,1),heights=c(1,1,1,1));
} else if (panels == 3) {
    png(paste(args[1],"snapshot_",run,"_",time,".png",sep=""),width = 6, height = 4, units="in",res=200);
    par(mar=c(1,1,2,1),oma=c(1,0,0,0),cex=0.5);
    layout(matrix(c(1:3),1,3,byrow = TRUE),widths=c(1,1,1),heights=c(1,1,1));
} else {
    png(paste(args[1],"snapshot_",run,"_",time,".png",sep=""),width = 6, height = 6, units="in",res=200);
    par(mar=c(1,1,2,1),oma=c(1,0,0,0),cex=0.5);
    layout(matrix(c(1:2),1,2,byrow = TRUE),widths=c(1,1),heights=c(1,1));
}

site_scale<-50; #how many microns is one site
cyt_disp_thresh<-7.406*1.25e-4;
#cyt_effect<-function(dose,xmid=cyt_ic50,scal=-7.458*1.25e-4) (1-((1+exp(xmid/scal))/(1+exp((xmid-dose)/scal))));
cyt_effect<-function(dose,xmid=cyt_ic50) {
    return ((dose)/(xmid+dose));
}

L<-125;

#read in cell_state matrix first
cell_file<-paste(args[1],"cell_state_",run,"_",time,".csv",sep="");
print(paste("Reading file",cell_file));
data <- read.csv(cell_file, comment.char="#");
cell_state <- data.matrix(data[,1:ncol(data)]);

#0=uninf, 1=inf_np, 2=inf_p, 3=dead,4=cyto killed,5=trm killed,6=extra col
uninf<-sum(cell_state%in%c(0,1,2));
inf_np<-sum(cell_state==3);
inf_p<-sum(cell_state==4);
nat_death<-sum(cell_state==5);
trm_death<-sum(cell_state==7);
dead_cells<-sum(cell_state>=5);
temp_grid<-cell_state;
temp_grid[cell_state%in%c(0,1,2)]<-0;
temp_grid[cell_state%in%c(3)]<-1;
temp_grid[cell_state%in%c(4)]<-2;
temp_grid[cell_state%in%c(5)]<-3;
if (panels > 3) {
  cyt_death<-sum(cell_state==6);
  cols<-c('papayawhip','lightgreen','red','black','gray','purple',        
			  'orange');                     
  temp_grid[cell_state%in%c(6)]<-4;
  temp_grid[cell_state%in%c(7)]<-5;
  leg_text<-c(paste(uninf),paste(inf_np),paste(inf_p),paste(nat_death),paste(cyt_death),paste(trm_death))
  num_breaks<-6
} else {
  cols<-c('papayawhip','lightgreen','red','black','purple',        
			  'orange');                     
  temp_grid[cell_state%in%c(7)]<-4;
  leg_text<-c(paste(uninf),paste(inf_np),paste(inf_p),paste(nat_death),paste(trm_death))
  num_breaks<-5
}

image(temp_grid[(L/4):(3*L/4),(L/4):(3*L/4)],col=cols,breaks=c(-1:num_breaks),xaxt='n',yaxt='n',main='Cells', asp=1)
#leg_text<-c("empty","susc","pre-prod","prod inf","died","cyt killed","Trm killed")
legend('topright',leg_text,fill=cols,ncol=1,bg='white')
legend('bottomleft',paste("lesion=",round(calc_plaque_size(dead_cells,site_scale),2),'mm'))

#Panel 2: Virions
#read in free_virus matrix
#0,<10^2, <10^3, <10^4, <10^5, <10^6, <10^7,<10^8 and an extra colour
col_virus<-c('white','purple','blue','green','yellow','orange','red','darkred','black');
cols<-col_virus;
vl_file<-paste(args[1],"virus_state_",run,"_",time,".csv",sep="");
print(paste("Reading file",vl_file));
data <- read.csv(vl_file, comment.char="#");
free_virus <- data.matrix(data[,1:ncol(data)]);

virions<-sum(free_virus);
if(virions==0){
  log_virus<-0
  image(free_virus[(L/4):(3*L/4),(L/4):(3*L/4)],col='white',xaxt='n',yaxt='n',main='Virus', asp=1)
}else{
  log_virus<-log(virions)/log(10);
  temp_virus<-log(free_virus)/log(10)
  temp_virus[free_virus ==0]<-0;
  temp_virus[free_virus > 0 & log(free_virus)/log(10) <2]<-1;
  temp_virus[log(free_virus)/log(10) >=2 & log(free_virus)/log(10) <2.25]<-2;
  temp_virus[log(free_virus)/log(10) >=2.25 & log(free_virus)/log(10) <2.5]<-3;
  temp_virus[log(free_virus)/log(10) >=2.5 & log(free_virus)/log(10) <2.75]<-4;
  temp_virus[log(free_virus)/log(10) >=2.75 & log(free_virus)/log(10) <3]<-5;
  temp_virus[log(free_virus)/log(10) >=3 & log(free_virus)/log(10) <3.25]<-6;
  temp_virus[log(free_virus)/log(10) >=3.25]<-7;
  
  # image(free_virus,col=c('grey',colorRampPalette(c('blue','white'))(10)),xaxt='n',yaxt='n',main='Virus');
  #p<-sample(rainbow(7,start=0,end=6))
  #image(log(free_virus[(L/4):(3*L/4),(L/4):(3*L/4)])/log(10),col=cols,breaks=c(-1:8),xaxt='n',yaxt='n',main='Virus');
  image(temp_virus[(L/4):(3*L/4),(L/4):(3*L/4)],col=cols,breaks=c(-1:8),xaxt='n',yaxt='n',main='Virus', asp=1);
  #leg_text<-c('0','<10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6','>=10^6')
  #legend('topright',leg_text,fill=cols,bty='n',ncol=1)
}
legend('bottomleft',paste("log VL=",format(round(log_virus,2)),", max=",format(round(max_log_vl,2))," at time=",format(round(max_time,2)),sep=""))

if (panels > 2) {
    #Panel 3: T cells +cytokine
    #read in t-cell matrix AND cytokine matrix
    #Vacant, PAT HSV, ACT HSV, PAT BYST, ACT BYST, dead T-cell and an extra colour
    tcell_file<-paste(args[1],"tcell_state_",run,"_",time,".csv",sep="");
    print(paste("Reading file",tcell_file));
    data <- read.csv(tcell_file, comment.char="#");
    t_cell_state <- data.matrix(data[,1:ncol(data)]);
    cyto_file<-paste(args[1],"cyto_state_",run,"_",time,".csv",sep="");
    print(paste("Reading file",cyto_file));
    data <- read.csv(cyto_file, comment.char="#");
    cytokine <- data.matrix(data[,1:ncol(data)]);

    no_tcell<-sum(t_cell_state%in%c(0,5));
    temp_tcells<-t_cell_state;
    temp_tcells[t_cell_state%in%c(0,5)]<-0;
    if (panels > 3 & cyto_act == 1) {
      col_tcells<-c('white','green','red','blue','purple','black');
      temp_tcells[t_cell_state==6&cytokine==0]<-1;
      temp_tcells[t_cell_state==7|(t_cell_state==6&cytokine>0)]<-2;
      temp_tcells[t_cell_state==8&cytokine==0]<-3;
      temp_tcells[t_cell_state==9|(t_cell_state==8&cytokine>0)]<-4;
      pat_hsv<-sum(t_cell_state==6&cytokine==0);
      act_hsv<-sum(t_cell_state==7|(t_cell_state==6&cytokine>0));
      pat_byst<-sum(t_cell_state==8&cytokine==0)
      act_byst<-sum(t_cell_state==9|(t_cell_state==8&cytokine>0));
      leg_text<-c(paste(no_tcell),paste(pat_hsv),paste(act_hsv),paste(pat_byst),paste(act_byst))
      num_breaks<-5
    } else {
      col_tcells<-c('white','green','red','blue','black');
      temp_tcells[t_cell_state==6]<-1;
      temp_tcells[t_cell_state==7]<-2;
      temp_tcells[t_cell_state==8]<-3;
      pat_hsv<-sum(t_cell_state==6);
      act_hsv<-sum(t_cell_state==7);
      pat_byst<-sum(t_cell_state==8)
      leg_text<-c(paste(no_tcell),paste(pat_hsv),paste(act_hsv),paste(pat_byst))
      num_breaks<-4
    }
    image(temp_tcells[(L/4):(3*L/4),(L/4):(3*L/4)],col=col_tcells,breaks=c(-1:num_breaks),xaxt='n',yaxt='n',main='T cells', asp=1);

    cols<-col_tcells;
    legend('topright',leg_text,fill=cols,bg='white');

    if (panels > 3) {
	#Panel 4: Cytokine
	#0=Vacant, 1=dead, 2=susc, 3=protected, 4=inf,5=dying,5=extra col
	temp_cells<-cell_state;
	if (plot_cytokines == 0) {
	    col_cyto_cells<-c('white','black','papayawhip','green','yellow','red','orange');
	    cols<-col_cyto_cells;
	    temp_cells[cell_state%in%c(5,6,7)]<-1; # dead
	    temp_cells[cell_state%in%c(1)&cyt_effect(cytokine)<0.5]<-2; # susc
	    temp_cells[cell_state%in%c(1)&cyt_effect(cytokine)>=0.5]<-3; # protected
	    temp_cells[cell_state%in%c(2)]<-4; # unprotected
	    temp_cells[cell_state%in%c(3,4)]<-5; # inf
	    image(temp_cells[(L/4):(3*L/4),(L/4):(3*L/4)],col=cols,breaks=c(-1:6),xaxt='n',yaxt='n',main='Cytokine effects', asp=1);
	    leg_text<-c("empty","dead cell","susc cell","protected","no protect","infected")
	} else {
	    tot_cytos<-sum(cytokine);
	    temp_cytos<-apply(cytokine,MARGIN=c(1,2),FUN=cyt_effect);
	    if (length(which(cell_state%in%c(0,3,4,5,6,7)))) {
	      temp_cytos[cell_state%in%c(0,3,4,5,6,7)]<- -temp_cytos[cell_state%in%c(0,3,4,5,6,7)]; # unprotected
	    }
	    temp_cytos<-0.5*(temp_cytos+1);
	    if (tot_cytos == 0) {
		cols<-"white";
	    } else {
		cols<-colorRampPalette(c("blue","white","white","red"))(100);
	    }

	    image(temp_cytos[(L/4):(3*L/4),(L/4):(3*L/4)],zlim=c(0,1),col=cols,xaxt='n',yaxt='n',main='Cytokine effects', asp=1);
	    #image(temp_cytos[(L/4):(3*L/4),(L/4):(3*L/4)],col=cols,breaks=c(-1:8),xaxt='n',yaxt='n',main='Cytokine effects');
	    #leg_text<-c("dead","cyt killed","trm killed","latent","prod inf","no effect","some effect",">ic50")
	    #legend('topright',leg_text,fill=cols,bty='n');
	}
	#types<-as.numeric(as.data.frame(table(temp_cells),stringsAsFactors=FALSE)[,1]);
	#image(temp_cells[(L/4):(3*L/4),(L/4):(3*L/4)],col=cols,breaks=c(-1:6),xaxt='n',yaxt='n',main='T cells');
	#legend('topright',leg_text[types+1],fill=cols[types+1],bty='n');
	#legend('topright',leg_text,fill=cols,bty='n');
    }
}

#Timestamp
mtext(paste('Episode',run,' Time =',time,'days'),side=1,outer=T);

dev.off()

print(paste("Inf cells =",inf_np+inf_p,"Log10 VL =",format(round(log_virus,2)),"ic50=",formatC(cyt_ic50,format='e',digits=2),sep=" "));
