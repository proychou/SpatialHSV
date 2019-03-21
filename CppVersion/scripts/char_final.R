#########################################################
### Reading in HSV data from simulation
#########################################################
#Rscript ../char_final.R $beta $viral_prod $viral_diff $cyt_diff $cyt_uptake $beta_ic50 $inf_ic50 $prod_ic50
args<-commandArgs(trailingOnly=T);
beta<-args[1]
viral_prod<-args[2]
viral_diff<-args[3]
cyt_diff<-args[4]
cyt_uptake<-args[5]
beta_ic50<-args[6]
inf_ic50<-args[7]
prod_ic50<-args[8]
old_scale<-as.numeric(args[9])
cyt_trm_ic50<-as.numeric(args[10]) 
cyt_inf_death<-as.numeric(args[11])
hsv_fract<-as.numeric(args[12])

data <- read.csv("results.csv", comment.char="#")
#episode,time,prodeptible,inf_nonprod,inf_prod,dead_cells,pat hsv Trm,act hsv Trm,pat byst Trm,act byst Trm,virions,log virions,max log VL,max time,viral_cells,cytokines,cyto_cells,gauss offset,avg E/T ratio,3x3 E/T ratio,5x5 E/T ratio,7x7 E/T ratio,9x9 E/T ratio,11x11 E/T ratio,percent 0s,closest HSV+,closest Tcell,total infected,log total virus,cytokine kills,Tcell kills,first beta antigen,first beta cycto,first sampled log VL,first sampled time,max sampled log VL,max sampled time

print (nrow(data))
write("sample,run,time,log VL", file="all_samples.csv",append=FALSE)

p_ids<-data[,1]
col_time<-data[,2]
log_virus<-data[,12]
tot_cyto<-data[,16]
startHSV<-data[,7]
actHSV<-data[,8]
startTcell<-data[,9]
et_ratio<-data[,19]
et_3x3<-data[,20]
et_5x5<-data[,21]
et_7x7<-data[,22]
et_9x9<-data[,23]
et_11x11<-data[,24]
perc_0s<-data[,25]
near_HSV<-data[,26]
near_Tcell<-data[,27]
tot_inf<-data[,28]
log_Vt<-data[,29]
cyto_kills<-data[,30]
beta_kills<-data[,31]

id_index<-vector(length=nrow(data))
epis<-vector(length=nrow(data))
avg_dur<-vector(length=nrow(data))
vl_max<-vector(length=nrow(data))
cyto_max<-vector(length=nrow(data))
first_pos<-vector(length=nrow(data))
avg_shed<-vector(length=nrow(data))
epi_rate<-vector(length=nrow(data))
vl_means<-vector(length=nrow(data))
vl_mean2<-vector(length=nrow(data))
vl_vars<-vector(length=nrow(data))
vl_measures<-vector(length=nrow(data))
tot_time<-vector(length=nrow(data))
all_epis_id<-vector(length=nrow(data))
all_epis_dur<-vector(length=nrow(data))
all_epis_max<-vector(length=nrow(data))
all_epis_cyto_max<-vector(length=nrow(data))
all_epis_max_time<-vector(length=nrow(data))
all_epis_max_scale<-vector(length=nrow(data))
#all_epis_mean_VL<-vector(length=nrow(data))
all_epis_firsts<-vector(length=nrow(data))
all_et_ratio<-vector(length=nrow(data))
all_et_3x3<-vector(length=nrow(data))
all_et_5x5<-vector(length=nrow(data))
all_et_7x7<-vector(length=nrow(data))
all_et_9x9<-vector(length=nrow(data))
all_et_11x11<-vector(length=nrow(data))
all_perc_0s<-vector(length=nrow(data))
all_startHSV<-vector(length=nrow(data))
all_startTcell<-vector(length=nrow(data))
all_near_HSV<-vector(length=nrow(data))
all_near_Tcell<-vector(length=nrow(data))
all_tot_inf<-vector(length=nrow(data))
all_log_Vt<-vector(length=nrow(data))
all_cyto_kills<-vector(length=nrow(data))
all_beta_kills<-vector(length=nrow(data))
all_beta_act_times<-vector(length=nrow(data))
total_epis=0
vls_this_pat=0
beta_act_time=-1

pat_cnt=0

last_id=0

last_log_vl=0
next_last_log_vl=0
epi_first_pos=0

first_time=0
last_time=0
next_last_time=0

next_last_i=0
last_i=0
first_i = 1

epi_start=0
cyto_max_this_epi=0
max_this_epi=0
max_time_this_epi=0
#mean_epi_VL=0
num_epi_VL=0
start_recording=0
sample=1

for (i in 1:nrow(data))
{
    if (is.na(log_virus[i]))
	log_virus[i]=0

    if (p_ids[i]!=last_id)
    { 
	# here we count on-going shedding as episode (unlike cohort data)
	# since we may not allow enough time (or immune syst) to resolve
	if (pat_cnt > 0 && last_log_vl != 0 && epi_start>0)
	{
	  print(paste("ID:",last_id,"had",epis[pat_cnt],"episodes"))
	  epis[pat_cnt]=epis[pat_cnt]+1  
	  #if (epis[pat_cnt]==0)
	  #    epis[pat_cnt]=1
	  avg_dur[pat_cnt]=(avg_dur[pat_cnt]*(epis[pat_cnt]-1) + last_time-epi_start+0.25*runif(1)+0.25*runif(1)) / epis[pat_cnt]
	  if (avg_dur[pat_cnt]==0)
		avg_dur[pat_cnt]=0.25*runif(1)+0.25*runif(1)
	  #first_pos[pat_cnt]=(first_pos[pat_cnt]*(epis[pat_cnt]+1) + epi_first_pos) / epis[pat_cnt]
	  if (max_this_epi > vl_max[pat_cnt])
	      vl_max[pat_cnt]=max_this_epi

	  if (cyto_max_this_epi > cyto_max[pat_cnt])
	      cyto_max[pat_cnt]=cyto_max_this_epi

           total_epis=total_epis+1  
           all_epis_dur[total_epis]=last_time-epi_start+0.25*runif(1)+0.25*runif(1)
           all_epis_cyto_max[total_epis]=cyto_max_this_epi
           all_epis_max[total_epis]=max_this_epi
           all_epis_max_time[total_epis]=max_time_this_epi
           all_epis_max_scale[total_epis]=scale_at_max
	   #all_epis_mean_VL[total_epis]=mean_epi_VL
           all_epis_firsts[total_epis]=epi_first_pos
           all_epis_id[total_epis]=i-1
	   all_et_ratio[total_epis]=et_ratio[i-1]
	   all_et_3x3[total_epis]=et_3x3[i-1]
	   all_et_5x5[total_epis]=et_5x5[i-1]
	   all_et_7x7[total_epis]=et_7x7[i-1]
	   all_et_9x9[total_epis]=et_9x9[i-1]
	   all_et_11x11[total_epis]=et_11x11[i-1]
	   all_perc_0s[total_epis]=perc_0s[i-1]
	   all_near_HSV[total_epis]=near_HSV[i-1]
	   all_near_Tcell[total_epis]=near_Tcell[i-1]
	   all_tot_inf[total_epis]=tot_inf[i-1]
	   all_log_Vt[total_epis]=log_Vt[i-1]
	   all_cyto_kills[total_epis]=cyto_kills[i-1]
	   all_beta_kills[total_epis]=beta_kills[i-1]
	   all_beta_act_times[total_epis]=beta_act_time
	   print(paste("Old ID: epi=",epis[pat_cnt],"end=",last_time,"1st=",all_epis_firsts[total_epis],"max=",all_epis_max[total_epis],"at t=",all_epis_max_time[total_epis],"dur=",all_epis_dur[total_epis]))
	}
	if (pat_cnt > 0 && epi_start > 0)
	{
	    first_pos[pat_cnt]=first_pos[pat_cnt] / epis[pat_cnt]
	    if (last_time > first_time)
	    {
		avg_shed[pat_cnt]=100*avg_shed[pat_cnt]/ (last_time - first_time)
		if (vl_measures[pat_cnt] > 0)
		    vl_means[pat_cnt]=vl_means[pat_cnt]/ vl_measures[pat_cnt]
		print (paste("Patient",pat_cnt,"first_i=",first_i,"i=",i-1))
		print(log_virus[first_i])
		print(log_virus[first_i:(i-1)])
		vl_mean2[pat_cnt]=mean(log_virus[first_i-1+which(log_virus[first_i:(i-1)]>0)])
		vl_vars[pat_cnt]=var(log_virus[first_i-1+which(log_virus[first_i:(i-1)]>0)])
		cis<-t.test(log_virus[first_i:(i-1)])
		cis$conf.int
	    }
	    else
	    {
		avg_shed[pat_cnt]=100
		vl_means[pat_cnt]=last_log_vl
		vl_mean2[pat_cnt]=last_log_vl
		vl_measures[pat_cnt]=1
		vl_vars[pat_cnt]=0
	    }
	    epi_rate[pat_cnt]=0
	    if (last_time > first_time)
		epi_rate[pat_cnt]=epis[pat_cnt]*365/(last_time - first_time)
	    tot_time[pat_cnt]=(last_time-first_time)
	    if (tot_time[pat_cnt]==0)
		tot_time[pat_cnt]=1
	
	    print (paste("Patient",pat_cnt,"avg shed=",avg_shed[pat_cnt],"epi rate=",epi_rate[pat_cnt],"avg first=",first_pos[pat_cnt],"avg dur=",avg_dur[pat_cnt],"avg log VL=",vl_means[pat_cnt],"nonzeros=",vl_measures[pat_cnt],"means2=",vl_mean2[pat_cnt],"tot_time=",tot_time[pat_cnt]))
	    
	}
	pat_cnt = pat_cnt+1
	beta_act_time=-1

#	reinitialize for the new patient...
	scale<-tot_inf[i]
	scale<-1
	if (scale > 10)
	    scale = 10
	if (log_virus[i]> 0)
	    #scaled_log_virus=log(scale*(10**log_virus[i]))/log(10)
	    scaled_log_virus=log_virus[i]
	else
	    scaled_log_virus=0

	cyto_max_this_epi=tot_cyto[i]
	max_this_epi=scaled_log_virus
	scale_at_max=scale
	max_time_this_epi=0
        epi_first_pos = scaled_log_virus
	avg_shed[pat_cnt]=0
	vl_means[pat_cnt]=0
        vl_measures[pat_cnt]=0
	id_index[pat_cnt]=i
	epis[pat_cnt]=0
	last_id = p_ids[i]
	first_i=i
        all_startHSV[total_epis+1]=startHSV[i]
        all_startTcell[total_epis+1]=startTcell[i]
        epi_start=0

	#mean_epi_VL=0
	print(paste("For new ID:",p_ids[i],"first_i=",first_i,"scaled first VL=",scaled_log_virus))
	write(paste(i,p_ids[i],col_time[i],scaled_log_virus,sep=","),file="hsv_samples.csv",append=TRUE)

	#if (log_virus[i] >0)
	if (scaled_log_virus >=2)
	{
	   epi_start=col_time[i]
	   first_pos[pat_cnt] = 10**epi_first_pos
	   #epis[pat_cnt]=epis[pat_cnt]+1  
	   print(paste("New ID:",last_id,"start=",epi_start,"scaled log VL=",scaled_log_virus))
	   vl_means[pat_cnt]=scaled_log_virus
	   vl_measures[pat_cnt]= 1
	   #mean_epi_VL=log_virus[i]
	   num_epi_VL=1
	}
	next_last_log_vl=0
	last_log_vl=scaled_log_virus
	next_last_time=0
	next_last_i=0
	last_i=i
	last_time=col_time[i]
	first_time=col_time[i]
	start_recording=1
    } 
    else
    {
	scale<-tot_inf[i]
	scale<-1
	if (scale > 10)
	    scale = 10
	if (log_virus[i]> 0)
	    #scaled_log_virus=log(scale*(10**log_virus[i]))/log(10)
	    scaled_log_virus=log_virus[i]
	else
	    scaled_log_virus=0

	if (start_recording)
	{
	    if (beta_act_time < 0 && actHSV[i] > 0) {
		beta_act_time=col_time[i]
		print(paste("HSV activation at t=",col_time[i]))
	    }
	    if (scaled_log_virus >=2 && last_log_vl < 2)
	    {
	       epi_start=col_time[i]
	       max_this_epi=scaled_log_virus
	       scale_at_max=scale
	       max_time_this_epi=0.25*runif(1)
	       print(paste("Same ID: start=",epi_start,"log VL=",log_virus[i]))
	       #epis[pat_cnt]=epis[pat_cnt]+1  
	       if (i != first_i)
		   avg_shed[pat_cnt]= avg_shed[pat_cnt] + (col_time[i] -last_time)/2
	       epi_first_pos=scaled_log_virus
	       first_pos[pat_cnt]=first_pos[pat_cnt] + 10**epi_first_pos
	       #mean_epi_VL=10**log_virus[i]
	       num_epi_VL=1
	    }
	    else if (scaled_log_virus >max_this_epi)
	    {
		max_this_epi=scaled_log_virus
	        scale_at_max=scale
		max_time_this_epi=col_time[i]-epi_start+0.25*runif(1)
	    }
	    print(paste("Checking at t=",col_time[i],"scale=",scale,"scaled log virus=",scaled_log_virus,"max=",max_this_epi,"at t=",max_time_this_epi,"epi start=",epi_start))

	    if (tot_cyto[i] >cyto_max_this_epi)
		cyto_max_this_epi=tot_cyto[i]

	    if (scaled_log_virus ==0 && last_log_vl != 0)
	    #if (log_virus[i] >=2 && last_log_vl < 2)
	    {
	       epis[pat_cnt]=epis[pat_cnt]+1
	       #mean_epi_VL=log(mean_epi_VL/num_epi_VL)/log(10)
	       num_epi_VL=1
	       avg_dur[pat_cnt]=(avg_dur[pat_cnt]*(epis[pat_cnt]-1) + (col_time[i]-epi_start+0.25*runif(1)+0.25*runif(1))) / epis[pat_cnt]
           all_epis_dur[total_epis]=last_time-epi_start+0.25*runif(1)+0.25*runif(1)
	      if (avg_dur[pat_cnt]==0)
		    avg_dur[pat_cnt]=0.25*runif(1)+0.25*runif(1)
	       #first_pos[pat_cnt]=(first_pos[pat_cnt]*(epis[pat_cnt]-1) + epi_first_pos) / epis[pat_cnt]
	      if (max_this_epi > vl_max[pat_cnt])
		  vl_max[pat_cnt]=max_this_epi
	      if (tot_cyto[i] >cyto_max_this_epi)
		  cyto_max_this_epi=tot_cyto[i]

	       avg_shed[pat_cnt]= avg_shed[pat_cnt] + (col_time[i] - last_time)/2
	       total_epis=total_epis+1  
	       all_epis_dur[total_epis]=(col_time[i]-epi_start+0.25*runif(1)+0.25*runif(1))
	       all_epis_cyto_max[total_epis]=cyto_max_this_epi
	       all_epis_max[total_epis]=max_this_epi
	       all_epis_max_time[total_epis]=max_time_this_epi
	       all_epis_max_scale[total_epis]=scale_at_max
	       #all_epis_mean_VL[total_epis]=mean_epi_VL
	       all_epis_firsts[total_epis]=epi_first_pos
	       all_epis_id[total_epis]=i
	       all_et_ratio[total_epis]=et_ratio[i]
	       all_et_3x3[total_epis]=et_3x3[i]
	       all_et_5x5[total_epis]=et_5x5[i]
	       all_et_7x7[total_epis]=et_7x7[i]
	       all_et_9x9[total_epis]=et_9x9[i]
	       all_et_11x11[total_epis]=et_11x11[i]
	       all_perc_0s[total_epis]=perc_0s[i]
	       all_near_HSV[total_epis]=near_HSV[i]
	       all_near_Tcell[total_epis]=near_Tcell[i]
	       all_tot_inf[total_epis]=tot_inf[i]
	       all_log_Vt[total_epis]=log_Vt[i]
	       all_cyto_kills[total_epis]=cyto_kills[i]
	       all_beta_kills[total_epis]=beta_kills[i]
	       all_beta_act_times[total_epis]=beta_act_time
	       max_this_epi=scaled_log_virus;
	       scale_at_max=scale
	       cyto_max_this_epi=tot_cyto[i];

	       print(paste("Curr ID: epi=",epis[pat_cnt],"start=",epi_start,"end=",col_time[i],"1st=",all_epis_firsts[total_epis],"max=",all_epis_max[total_epis],"at t=",all_epis_max_time[total_epis],"dur=",all_epis_dur[total_epis]))
	    }
	    
	    if (scaled_log_virus !=0 && last_log_vl != 0 &&  i != first_i)
	    #if (log_virus[i] >=2 && last_log_vl >= 2 &&  i != first_i)
	    {
		avg_shed[pat_cnt]= avg_shed[pat_cnt] + (col_time[i]- last_time)
	        #mean_epi_VL=mean_epi_VL+10**log_virus[i]
	        num_epi_VL=num_epi_VL+1
	    }
	    if (scaled_log_virus != 0)
	    #if (log_virus[i] >= 2)
	    {
		vl_means[pat_cnt]= vl_means[pat_cnt] + scaled_log_virus
		vl_measures[pat_cnt]= vl_measures[pat_cnt] + 1
	    }
	    this_id=p_ids[i]
	    write(paste(i,this_id,col_time[i],scaled_log_virus,sep=","),file="hsv_samples.csv",append=TRUE)
	}
	next_last_log_vl=last_log_vl
	last_log_vl=scaled_log_virus
	next_last_time=last_time
	last_time=col_time[i]
    }
}
#if (last_log_vl != 0 && start_recording > 0)
if (last_log_vl >= 2 && start_recording > 0)
{
  epis[pat_cnt]=epis[pat_cnt]+1

  avg_dur[pat_cnt]=(avg_dur[pat_cnt]*(epis[pat_cnt]-1) + (last_time-epi_start+0.25*runif(1)+0.25*runif(1))) / epis[pat_cnt]
  if (avg_dur[pat_cnt]==0)
	avg_dur[pat_cnt]=0.25*runif(1)+0.25*runif(1)
  #first_pos[pat_cnt]=(first_pos[pat_cnt]*(epis[pat_cnt]-1) + epi_first_pos) / epis[pat_cnt]
  if (max_this_epi > vl_max[pat_cnt])
      vl_max[pat_cnt]=max_this_epi

  if (tot_cyto[i] >cyto_max_this_epi)
      cyto_max_this_epi=tot_cyto[i]
   #mean_epi_VL=log(mean_epi_VL/num_epi_VL)/log(10)

   total_epis=total_epis+1  
   all_epis_dur[total_epis]=(last_time-epi_start+0.25*runif(1)+0.25*runif(1))
   print(paste("Last ID: end=",last_time,"epis=",total_epis,"dur=",all_epis_dur[total_epis]))
   all_epis_cyto_max[total_epis]=cyto_max_this_epi
   all_epis_max[total_epis]=max_this_epi
   all_epis_max_time[total_epis]=max_time_this_epi
   all_epis_max_scale[total_epis]=scale_at_max
   #all_epis_mean_VL[total_epis]=mean_epi_VL
   all_epis_firsts[total_epis]=epi_first_pos
   all_epis_id[total_epis]= nrow(data)
   all_et_ratio[total_epis]=et_ratio[nrow(data)]
   all_et_3x3[total_epis]=et_3x3[nrow(data)]
   all_et_5x5[total_epis]=et_5x5[nrow(data)]
   all_et_7x7[total_epis]=et_7x7[nrow(data)]
   all_et_9x9[total_epis]=et_9x9[nrow(data)]
   all_et_11x11[total_epis]=et_11x11[nrow(data)]
   all_perc_0s[total_epis]=perc_0s[nrow(data)]
   all_near_HSV[total_epis]=near_HSV[nrow(data)]
   all_near_Tcell[total_epis]=near_Tcell[nrow(data)]
   all_tot_inf[total_epis]=tot_inf[nrow(data)]
   all_log_Vt[total_epis]=log_Vt[nrow(data)]
   all_cyto_kills[total_epis]=cyto_kills[nrow(data)]
   all_beta_kills[total_epis]=beta_kills[nrow(data)]
   all_beta_act_times[total_epis]=beta_act_time
}

last_shed=avg_shed[pat_cnt]
avg_shed[pat_cnt]=100
last_pat_vl=vl_means[pat_cnt]
vl_means[pat_cnt]=last_log_vl
vl_mean2[pat_cnt]=last_log_vl
vl_vars[pat_cnt]=0
epi_rate[pat_cnt]=0
first_pos[pat_cnt]=first_pos[pat_cnt] / epis[pat_cnt]
if (last_time > first_time)
{
    epi_rate[pat_cnt]=epis[pat_cnt]*365/(last_time - first_time)
    avg_shed[pat_cnt]=100*last_shed/ (last_time - first_time)
    if (vl_measures[pat_cnt] > 0)
	vl_means[pat_cnt]=last_pat_vl/ vl_measures[pat_cnt]
    vl_mean2[pat_cnt]=mean(log_virus[first_i-1+which(log_virus[first_i:nrow(data)]>0)])
    vl_vars[pat_cnt]=var(log_virus[first_i-1+which(log_virus[first_i:nrow(data)]>0)])
    print (paste("Last Patient",pat_cnt,"last time",last_time,"first time",first_time,"shed time=",last_shed,"shed rate=",avg_shed[pat_cnt],"epi rate=",epi_rate[pat_cnt]))
    cis<-t.test(log_virus[first_i:nrow(data)])
    cis$conf.int
}
tot_time[pat_cnt]=(last_time-first_time)
if (tot_time[pat_cnt]==0)
    tot_time[pat_cnt]=1

print (paste("Last Patient",pat_cnt,"last time",last_time,"first time",first_time,"avg shed=",avg_shed[pat_cnt],"avg vl=",vl_means[pat_cnt]))

print (paste(pat_cnt," runs"))
write("Runs,total time,episodes,avg days,peak log VL,avg shed,avg log VL,log VL mean,log VL variance,avg first pos,epi rate,max cytokine",file="run_summary_nz.csv",append=FALSE)
for (i in 1:pat_cnt)
{
    if (vl_means[i] > 0)
    {
	this_id=p_ids[id_index[i]]
	write(paste(this_id,tot_time[i],epis[i],avg_dur[i],vl_max[i],avg_shed[i],vl_means[i],vl_mean2[i],vl_vars[i],log(first_pos[i])/log(10),epi_rate[i],cyto_max[i],sep=","),file="run_summary_nz.csv",append=TRUE)
    }
}
print (paste(total_epis," episodes"))
write("infectivity,viral prod,viral diff,cyto diff,cyto uptake,beta ic50,inf ic50,beta ic50,scale_at_max,cyt_trm_ic50,cyt_inf_death,hsv_fract,Episode,run,duration,scaled peak logVL,peak time,first pos VL,avg E/T ratio,3x3 E/T ratio,5x5 E/T ratio,7x7 E/T ratio,9x9 E/T ratio,11x11 E/T ratio,percent 0s,start HSV+,start BYST,nearest HSV+,nearest BYST,max cytokine,total infected,log total virus,cytokine kills,Tcell kills,Tcell act time", file="all_episodes.csv",append=FALSE)
for (i in 1:total_epis)
{
    write(paste(beta,viral_prod,viral_diff,cyt_diff,cyt_uptake,beta_ic50,inf_ic50,prod_ic50,all_epis_max_scale[i],cyt_trm_ic50,cyt_inf_death,hsv_fract,i,p_ids[all_epis_id[i]],all_epis_dur[i],all_epis_max[i],all_epis_max_time[i],all_epis_firsts[i],all_et_ratio[i],all_et_3x3[i],all_et_5x5[i],all_et_7x7[i],all_et_9x9[i],all_et_11x11[i],all_perc_0s[i],all_startHSV[i],all_startTcell[i],all_near_HSV[i],all_near_Tcell[i],all_epis_cyto_max[i], 
   all_tot_inf[i], all_log_Vt[i], all_cyto_kills[i],
   all_beta_kills[i],all_beta_act_times[i],sep=","),file="all_episodes.csv",append=TRUE)
}
