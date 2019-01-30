# Spatial model of HSV in skin (2D)
# If running for the first time on mac, may need to install X11 for visualization
# Pavitra Roychoudhury
# Fred Hutchinson Cancer Research Center
# Latest update: 14-May-17


# INITIALIZATION: lines to run when this file is sourced
rm(list=ls()); 
library(RColorBrewer);

# 0.  MAIN FUNCTION-only used if there is no wrapper file
Spatial_main<-function(visualize='off'){
	options(error=recover);
	settings<-get_default_settings();
	settings$place_virus<-list(type='infected producer',num=1);
	parameters<-get_default_parameters();
	
	#modify settings and parameters here or use a wrapper, e.g. Spatial_wrapper.R
	settings$visualize<-visualize;
	
	#initialize and run sim
	settings<-init_grid(settings,parameters);
	simulation(settings,parameters);
}

## FUNCTIONS: all the functions that are needed for the simulation are below this line

# 1a.  Get settings for the simulation
get_default_settings<-function(){
	settings<-list(L=200);                          #Specify grid dimensions
	settings$dirname<-'Results/';   								#Directory to store results
	settings$out_fname<-'results.out';							#File to store cell and virus counts
	settings$sim_length<-5;											#Length of simulation in hours
	settings$save_grids<-TRUE;          #Save grids at the end of the sim? 
	
	#Visualization
	settings$visualize<-'off';                         #Options: 'screen','pdf','off','gif'
	settings$col_cells<-c('grey95','papayawhip','orange','turquoise','peachpuff4',        #0=Vacant, 1=uninf, 2=inf_np, 3=inf_p, 3=dead,
	                      'peachpuff2','orange','turquoise4');                     #5=extra col, 6=extra col, 7=extra col 
	settings$col_tcells<-c('grey95','forestgreen','deeppink','black',
	                       'green','khaki1','blue');#Vacant, TRM, TEM, dead cell, T-cell + cytokine, no T-cell+cytokine, and an extra colour
	settings$refr_freq<-25;												#Refresh frequency (hrs)
	settings$site_scale<-50; #how many microns is one site
	
	#Place virus
	#settings$place_virus<-list(type='single point',num=1);				#Start with single virus at the center of grid
	#settings$place_virus<-list(type='random',num=10);					#Start with num viruses, randomly scattered
	#settings$place_virus<-list(type='none');							#No free virus at start of simulation
	settings$place_virus<-list(type='infected producer',num=1);
	
	#TRM 
	settings$trm_init<-'gaussian'; #random or gaussian
	settings$trm_move_method<-'random'; #'random', 'directed' or 'levy'
	
	return(settings);
}

# 1b. Get default parameters 
get_default_parameters<-function(){
	parameters<-list(); 
	
	#Rates (per site or cell per day)
	parameters$max_rate<-25;						  #Sets time resolution, should be greater than or equal to sum of all rates for a given state
	
	#Uninfected cells
	parameters$uninf_cell_death<-0;	#Rate of uninfected cell death 
	parameters$infectivity_free<-0.1;				#Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
	parameters$cell_div<-0.077; 					#Rate of generation of new susceptible cells (repop of empty sites)
	
	#Infected cells
	parameters$inf_cell_death<-1.25;			#Rate of infected cell death due to infection
	parameters$viral_lag<-8;						  #Rate at which inf non-producer becomes viral producing cell (per cell per day)
	parameters$vir_prod_rate_mean<-100000;	#Rate at which an infected cell produces free virus (per cell per day, mean and sd)
	parameters$vir_prod_rate_sd<-10000;
	
	#Virus
	parameters$diff_rate_virus<-12; 						  #Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
	parameters$freevirus_decay<-8.8;		#Rate of decay of free virus (per day)
	parameters$neur_drip<-0;						  #Rate of arrival of virus from neurons (per site per day)

	#T-cells
	parameters$tem_trafficking<-0.01;				#Rate of arrival of TEM from LN (per site per day)
	parameters$trm_decay<-0;	#Death rate of TRMs and TEMs in absence of Ag (per site per day, default = 1E-3)
	parameters$trm_conv_rate<-0.05;					#Rate at which TEM becomes TRM in the absence of infection (per site per day, default = 1E-4)
	parameters$trm_dbl<-1;				#Proliferation (doubling) rate of TRMs/TEMs (per cell per day)
	parameters$tem_delay<-1;      #Delay in the arrival of TEM (per day, default = 1 = 24h delay)
	parameters$tcell_killing_rate<-96;      #Rate of target killing by T-cells (both aTRM and TEM)
	parameters$trm_motility<-5;					  #T-cell motility (jumps/movements per day)
	parameters$trm_levy_alpha<-2.5;       #Levy exponent (should be between 1 and 3)
	parameters$trm_levy_max<-10;          #Maximum jump per timestep (sites)

	#Cytokine-related parameters
	parameters$cyt_secretion_rate_tcell_mean<-0.05;     #Rate of cytokine secretion by T-cells, mean and sd (pg/cell/day)
	parameters$cyt_secretion_rate_tcell_sd<-0.005; 
	parameters$cyt_secretion_rate_infcell_mean<-0.005;     #Rate of cytokine secretion by infected cells, mean and sd (pg/cell/day)
	parameters$cyt_secretion_rate_infcell_sd<-0.0005; 
	parameters$cyt_diff_rate<-65.5;      #Cytokine diffusion/spread rate in sites/day
	parameters$cyt_decay_rate<-2.16;       #Rate of clearance of cytokine (/day)
	parameters$cyt_ic50<-7.406*1.25e-4;
	parameters$cyt_effect<-function(dose,xmid=parameters$cyt_ic50,scal=-7.458*1.25e-4) (1-((1+exp(xmid/scal))/(1+exp((xmid-dose)/scal))));             #Effect of cytokine-based activation given conc of cytokine at site (pg/ml). xmid and scal were estimated from Jia's experiments (see folder IFNG_jia) and multiplied to convert from units/ml to pg/cell
	parameters$cyt_effect_infectivity<-0; #Do cytokines affect infectivity: yes (1) or no (0)
	parameters$cyt_effect_virprod<-0; #Do cytokines affect virus production rate: yes (1) or no (0)
	parameters$cyt_effect_lag<-0; #Max effect of cytokines on lag (multiplier, remember that lag here is a rate)
	parameters$cyt_effect_uninfcelldeath<-48; #Max death rate of uninfected cells in the presence of cytokine, 0 means no effect
	parameters$cyt_effect_infncelldeath<-48; #Max death rate of pre-productively infected cells cells in the presence of cytokine, 0 means no effect
	parameters$cyt_effect_infpcelldeath<-48; #Max death rate of infected cells in the presence of cytokine, 0 means no effect
	
	#Probabilities, fractions, etc
	parameters$host_density<-1;           #Fraction of sites occupied by susceptible cells at start	
	parameters$diff_range<-3;						  #n x n site, default n = 3, i.e. immediate Moore neighbourhood. 
	parameters$tcell_density<-0.1;				#Density (fraction of total sites) of T-cells at the start of the simulation
	parameters$tcell_detection_range<-3;	#Spatial range of detection of antigen from infected cell (n x n), default n = 3

	#Other
	parameters$pat_trm_activatedby<-'infp';		#what activates a T-cell: 'inf','trm', 'act_tem' or 'any'
	parameters$tcell_killing_target<-'infp';	#What are the T-cells killing? 'infp','inf','all'. 'inf' means both early and late
	
	return(parameters);
}



# 2. Initialize grid and other arrays
init_grid<-function(settings,parameters){
	print('Starting simulation.. hold on to your hat!');

	#Compute simulation length in number of updates
	settings$tot_updates<-settings$sim_length*parameters$max_rate;
	
	#Create arrays
	cell_state<<-array(data=as.integer(0),dim=c(settings$L,settings$L));
	t_cell_state<<-array(data=as.integer(0),dim=c(settings$L,settings$L));
	free_virus<<-array(data=as.integer(0),dim=c(settings$L,settings$L));											
	cytokine<<-array(data=0,dim=c(settings$L,settings$L));
	
	#Directories
	print(paste('Results directory =',settings$dirname));
	dir.create(settings$dirname,showWarnings=FALSE);
	save(parameters,file=paste(settings$dirname,'/simulationparameters',sep=''));
	save(settings,file=paste(settings$dirname,'/simulationsettings',sep=''));
	if(settings$visualize=='gif')dir.create(paste(settings$dirname,'/imgs',sep=''),showWarnings=FALSE)
	
	#Dynamic variables - this stores time elapsed and population counts
	dynamic<<-list(num_updates=0,time=0,
				   susceptible_cells=0,infected_nonprodcells=0,infected_prodcells=0,dead_cells=0,
				   virions=0,trm=0,tem=0,dead_tcells=0,cytokine=0,max_vl=0,max_vl_time=0,viral_cells=0,cyto_cells=0);
	output_results(settings,dynamic,T); #write header
	
	#Load T-cell movement function based on setting
	set_trm_move_function(settings$trm_move_method);
	
	#Visualization setup
	settings$cyt_disp_thresh<-parameters$cyt_ic50; #Above what conc should cyt be displayed in the t-cell panel
	if(settings$visualize!='off'){
	  if(settings$visualize=='pdf'){
	    pdf(file=paste(settings$dirname,'/Results.pdf',sep=''),
	        width=8,height=8); 
	  }else if(settings$visualize=='screen'){
	    quartz(width=6,height=6);   
	  }else if(settings$visualize=='gif'){
	    # library(animation)
	    # png(file=paste(settings$dirname,'/imgs/Results.png',sep=''))
	  }else{
	    stop('Not a valid setting for "visualize"')
	  }   	
	  par(mar=c(1,1,2,1),oma=c(1,0,0,0));
	  layout(matrix(c(1:4),2,2,byrow = TRUE),widths=c(1,1,1,1),heights=c(1,1,1,1));
  }
	gc(); return(settings)
}

# 2. SIMULATION
simulation<-function(settings,parameters){
	
	gaussian_mask_diffusion<-compute_gaussian(parameters$diff_range); 
	 
	place_host(parameters); #Place host cells
	if(settings$place_virus$type!='none') place_virus(settings); #Place virus	
	place_tcells(settings,parameters); #Place T-cells
	
	dynamic$num_updates<<-0;
	compute_totals(dynamic,parameters);
	output_results(settings,dynamic,F);
	if(settings$visualize!='off'){
		visualize(settings); 
	}
	next_screen_grab<-dynamic$num_updates+settings$refr_freq;
	
	while(dynamic$num_updates<settings$tot_updates){
		
		#Generate random order of sites to update
		update_order<-sample(seq(1,settings$L^2),settings$L^2,replace=TRUE)
		row_order<-((update_order-1)%%settings$L)+1;
		col_order<-floor((update_order-1)/settings$L)+1;
		dice_cells_all<-runif(settings$L^2,min=0,max=1);
		dice_v_all<-runif(settings$L^2,min=0,max=1);
		dice_t_all<-runif(settings$L^2,min=0,max=1);
		dice_cyt_all<-runif(settings$L^2,min=0,max=1);


		#Update sites in random order (same as picking a random site to update one at a time)
		for(rep in 1:settings$L^2){
			i<-row_order[rep];
			j<-col_order[rep];
			dice_cells<-dice_cells_all[rep];
			dice_v<-dice_v_all[rep];
			dice_t<-dice_t_all[rep];
			dice_cyt<-dice_cyt_all[rep];

			
			##Any site on grid
			##------------------
			##Viral drip from neuron 
			if(dice_cells<=parameters$neur_drip/parameters$max_rate){
				neur_drip(i,j);
			
			}else{
			
			  ##Cytokine layer
			  ## 0: no cytokine, 1: cytokine present
			  ##--------------
			  #State 1: cytokine present
			  if(cytokine[i,j]>0){
			    #Event 0a: Decay
			    cyt_decay(i,j,parameters);
			    
			    #Diffusion of cytokine
			    if(dice_cyt<=parameters$cyt_diff_rate/parameters$max_rate){
			      diffusion_cyt(settings,i,j,gaussian_mask_diffusion,parameters);}	
			  }
			  
			  
				##T-cells - these are antigen-specific cells
			  ## 0: no T-cell present
			  ## 1: TRM
			  ## 2: TEM
			  ## 3: dead T-cell
				##--------------------------------------------
				#State 0 or 3: no T-cell present or occupied by dead T-cell
				if(t_cell_state[i,j]==0||t_cell_state[i,j]==3){
			
					#TEM arrives from LN (trafficking) if inf cell in detection nbhd
					if(dynamic$time>(1/parameters$tem_delay)&&(
					   (detect_inf(settings,parameters$tcell_detection_range,i,j,what='inf')&&
					    dice_t<=parameters$tem_trafficking/parameters$max_rate)||
					   (detect_inf(settings,parameters$tcell_detection_range,i,j,what='t_cell')&&
					    dice_t<=parameters$cyt_effect(dose=cytokine[i,j])*parameters$tem_trafficking/parameters$max_rate))){
						t_cell_state[i,j]<<-2;
					}
			
			
				#State 1: TRM/TEM present
				}else if(t_cell_state[i,j]%in%c(1,2)){
				  
			    #Case a: Antigen present
				  if(detect_inf(settings,parameters$tcell_detection_range,i,j,
				                what=parameters$pat_trm_activatedby)){
				    
				    #Event 1a.0: produce cytokine 
				    cytokine[i,j]<<-cytokine[i,j]+max(0,rnorm(1,mean=parameters$cyt_secretion_rate_tcell_mean/parameters$max_rate,
				                                              sd=parameters$cyt_secretion_rate_tcell_sd/parameters$max_rate));
				  
				    #Event 1a.1: proliferate
				    if(dice_t<=parameters$trm_dbl/parameters$max_rate){
				      tcell_div(i,j,settings);
				    }
				    
				  #Case b: Antigen absent
				  }else{
				    #Event 1b.0: produce cytokine depending on amount of cytokine at site
				    # cytokine[i,j]<<-cytokine[i,j]+parameters$cyt_effect(dose=cytokine[i,j])*
				    #   max(0,rnorm(1,mean=parameters$cyt_secretion_rate_tcell_mean/parameters$max_rate,
				    #               sd=parameters$cyt_secretion_rate_tcell_sd/parameters$max_rate));
				    
				    #Event 1b.1: decay/die
				    if(dice_t<=parameters$trm_decay/parameters$max_rate){
				      t_cell_state[i,j]<<-3;
				      
				    #Event 1b.2: move (only TRM)
				    }else if(t_cell_state[i,j]==1 && 
				             dice_t<=(parameters$trm_decay+parameters$trm_motility)/parameters$max_rate){
				      move_tcell(i,j,settings,parameters);
				    }
				  }
				}
							
				##Host cells
			  ## 0: vacant site
			  ## 1: susceptible keratinocyte
			  ## 2: infected, not producing virus
			  ## 3: infected, producing virus
			  ## 4: dead keratinocyte
				##--------------------------------------------
				#State 0 or 4: vacant site or occupied by dead cell
				if(cell_state[i,j]==0||cell_state[i,j]==4){
				
					#Event 0.1: occupation by susceptible
					if(dice_cells<=(parameters$neur_drip+parameters$cell_div)/parameters$max_rate){
						cell_state[i,j]<<-1;
					}
				
			
				#State 1: uninfected cell 
				}else if(cell_state[i,j]==1){
			
					#Event 1.1: infection due to free virus at the site
					if(free_virus[i,j]>0){
					  #with cytokine
					  if((parameters$cyt_effect_infectivity==1&&dice_cells<=(parameters$neur_drip+
					                    (1-parameters$cyt_effect(dose=cytokine[i,j]))*
					                    parameters$infectivity_free*free_virus[i,j])/parameters$max_rate)){
					    cell_state[i,j]<<-2; #becomes an infected cell
					    free_virus[i,j]<<-free_virus[i,j]-1;
					    cytokine[i,j]<<-0;
					  #without cytokine
					  }else if(parameters$cyt_effect_infectivity==0&&
					           dice_cells<=(parameters$neur_drip+parameters$infectivity_free*free_virus[i,j])/
					      parameters$max_rate){
					    cell_state[i,j]<<-2; #becomes an infected cell
					    free_virus[i,j]<<-free_virus[i,j]-1;	 
					  }

					#Event 1.2: apoptosis
					}else if(parameters$cyt_effect_uninfcelldeath==0&&
					         dice_cells<=(parameters$neur_drip+parameters$uninf_cell_death)/parameters$max_rate){
					  cell_state[i,j]<<-4;
					}else if(parameters$cyt_effect_uninfcelldeath>0&&
					         dice_cells<=(parameters$neur_drip+(parameters$cyt_effect(
					           dose=cytokine[i,j]))*parameters$cyt_effect_uninfcelldeath)/parameters$max_rate){
					  cell_state[i,j]<<-4; 
					  cytokine[i,j]<<-0;
					} 
			
				#State 2: infected cell (non-virus producing, not presenting Ag)
				}else if(cell_state[i,j]==2){					
				  
				  #Event 2.1: become a virus-producing cell
				  if(parameters$cyt_effect_lag==0&&
				      dice_cells<=(parameters$neur_drip+parameters$viral_lag)/parameters$max_rate){
				    cell_state[i,j]<<-3; #becomes virus-producing infected cell
				  
				  }else if(parameters$cyt_effect_lag>0&&
				      dice_cells<=(parameters$neur_drip+(1-parameters$cyt_effect(dose=cytokine[i,j]))
				                   *parameters$viral_lag)/parameters$max_rate){
				    cell_state[i,j]<<-3; #becomes virus-producing infected cell
				    cytokine[i,j]<<-0;
						
					#Event 2.2: dies 
				  }else if(parameters$cyt_effect_infncelldeath==0&&dice_cells<=(
				    parameters$neur_drip+parameters$viral_lag+parameters$uninf_cell_death)/parameters$max_rate){
				    cell_state[i,j]<<-4; #death of infected cell
				    
				  }else if(parameters$cyt_effect_infncelldeath>0&&dice_cells<=(parameters$neur_drip+(
				             parameters$cyt_effect(dose=cytokine[i,j]))*parameters$cyt_effect_infncelldeath)
				           /parameters$max_rate){
				    cell_state[i,j]<<-4; #death of infected cell		
				    cytokine[i,j]<<-0;
				    
				  #Event 2.3: killed by a T-cell
				  }else if(parameters$tcell_killing_target=='inf'){
				    num_t_cells<-detect_inf(settings,3,i,j,what='t_cell',
				                            return_pos='sum');
				    
				    if(num_t_cells>0&&dice_cells<=(parameters$neur_drip+parameters$viral_lag+parameters$uninf_cell_death+
				                    parameters$tcell_killing_rate*num_t_cells)/parameters$max_rate){
				      cell_state[i,j]<<-4; 
				    }
				  }
				  
				#State 3: infected cell (producing virus, presenting Ag)
				}else if(cell_state[i,j]==3){
				  
					num_t_cells<-detect_inf(settings,3,i,j,what='t_cell',
					                        return_pos='sum');
					
					#Event 3.0: secrete cytokine
					cytokine[i,j]<<-cytokine[i,j]+
					  max(0,rnorm(1,mean=parameters$cyt_secretion_rate_infcell_mean/parameters$max_rate,
					                                                sd=parameters$cyt_secretion_rate_infcell_sd/parameters$max_rate));

					#Event 3.1: die
					if(parameters$cyt_effect_infpcelldeath==0&&
					   dice_cells<=(parameters$neur_drip+parameters$inf_cell_death)/parameters$max_rate){
					  cell_state[i,j]<<-4; #death of infected cell
					  
					}else if(parameters$cyt_effect_infpcelldeath>0&&dice_cells<=(parameters$neur_drip+(
					     parameters$cyt_effect(dose=cytokine[i,j]))*parameters$cyt_effect_infpcelldeath)
					     /parameters$max_rate){
						cell_state[i,j]<<-4; #death of infected cell	
						cytokine[i,j]<<-0;
					
					#Event 3.2: killed by a T-cell
					}else if(num_t_cells>0&&dice_cells<=(parameters$neur_drip+parameters$inf_cell_death+
					                      parameters$tcell_killing_rate*num_t_cells)/parameters$max_rate){
					  cell_state[i,j]<<-4; #death of infected cell
					
					#Event 3.3: produce free virus
					}else{
					  if(parameters$cyt_effect_virprod==0){
					    free_virus[i,j]<<-free_virus[i,j]+max(0,round(rnorm(1,mean=parameters$vir_prod_rate_mean/parameters$max_rate,
					                                                        sd=parameters$vir_prod_rate_sd/parameters$max_rate)));
					  }else{
					    free_virus[i,j]<<-free_virus[i,j]+(1-parameters$cyt_effect(dose=cytokine[i,j]))*
					      max(0,round(rnorm(1,mean=parameters$vir_prod_rate_mean/parameters$max_rate,
					                        sd=parameters$vir_prod_rate_sd/parameters$max_rate)));
					    cytokine[i,j]<<-0;
					  }
					}
				}
			
				##Cell-free virus
				##--------------------------------------------
				if(free_virus[i,j]>0){

				  #Event 0a: Decay
					viral_decay(i,j,parameters,dice_v);
				
					#Event 0b: Diffusion of remaining free virus at site
					if(dice_v<=parameters$diff_rate_virus/parameters$max_rate){
						diffusion_virus(settings,i,j,gaussian_mask_diffusion,parameters);}			
				
				}
			}			
			
		}
		
		dynamic$num_updates<<-dynamic$num_updates+1;
		compute_totals(dynamic,parameters);
		
		if(settings$visualize!='off'&&dynamic$num_updates==next_screen_grab){
			visualize(settings);
			next_screen_grab<-dynamic$num_updates+settings$refr_freq;
		}
		output_results(settings,results,F);
		if(dynamic$susceptible_cells==0 ||
		   (dynamic$virions==0&&parameters$neur_drip==0&&
		   dynamic$infected_nonprodcells==0&&dynamic$infected_prodcells==0)){
		  #quit if there is no possibility of virus or targets exhausted
		  break()
		}
	}
	
	end_simulation(settings,parameters);
}




#PLACE HOST cells
place_host<-function(parameters){
  if(parameters$host_density==1){
    cell_state[,]<<-1;
  }else{
    locs<-sample(seq(1,length(cell_state)),length(cell_state),replace=FALSE);
    locs<-locs[runif(length(locs),min=0,max=1)<parameters$host_density];
    cell_state[locs]<<-1;
    rm(locs);
  }
  gc();
}

#PLACE VIRUS
place_virus<-function(settings){
	if(settings$place_virus$num==0){
	  return(NULL); #do nothing and exit
	}else{
	  if(settings$place_virus$type=='single point'){ #place virus at the center of the grid
	    free_virus[settings$L/2,settings$L/2]<<-settings$place_virus$num; 
	    
	  }else if(settings$place_virus$type=='random'){ #virus spread randomly on grid
	    locs<-sample(seq(1,(settings$L^2)),settings$place_virus$num,replace=TRUE);
	    counts<-as.data.frame(table(locs),stringsAsFactors=F);
	    free_virus[as.integer(counts$locs)]<<-counts$Freq;
	    
	  }else if(settings$place_virus$type=='infected producer'){ #place infected cell(s) instead of virus
	    if(settings$place_virus$num==1){ #at the center of the grid
	      cell_state[settings$L/2,settings$L/2]<<-3; 
	    }else{ #randomly on the grid
	      locs<-sample(seq(1,(settings$L^2)),min(settings$place_virus$num,settings$L^2),replace=TRUE);
	      cell_state[locs]<<-3;
	    }
	  }
	}
}

#PLACE T-CELLS
place_tcells<-function(settings,parameters){
	if(parameters$tcell_density>0){
	  if(settings$trm_init=='random'){
	    locs<-sample(seq(1,length(t_cell_state)),length(t_cell_state),replace=FALSE);
	    locs<-locs[runif(length(locs),min=0,max=1)<parameters$tcell_density];
	    t_cell_state[locs]<<-1; rm(locs);
	  }else if(settings$trm_init=='gaussian'){
	    gausskern<-compute_gaussian(settings$L,sig=settings$L/3,amp=parameters$tcell_density);
	    gausskern<-parameters$tcell_density*gausskern/max(gausskern);
	    t_cell_state[runif(length(gausskern),0,1)<gausskern]<<-1;
	  }
	}
}

#NEURONAL DRIP
neur_drip<-function(i,j){
	free_virus[i,j]<<-free_virus[i,j]+1;
}

#Calculates plaque diameter in mm
calc_plaque_size<-function(settings){
  with(dynamic,
       return(sqrt(dead_cells)*settings$site_scale/1000));
}

#VISUALIZE
visualize<-function(settings){
  if(settings$visualize=='gif'){
    png(file=paste(settings$dirname,'/imgs/',
                   round(dynamic$time,digits=1),'.png',sep='')) ;
    par(mar=c(1,1,2,1),oma=c(1,0,0,0));
    layout(matrix(c(1:4),2,2,byrow = TRUE),widths=c(1,1,1,1),heights=c(1,1,1,1));
  }
  
  #Panel 1: Host cells 
  cols<-settings$col_cells;
  types<-as.numeric(as.data.frame(table(cell_state),stringsAsFactors=FALSE)[,1]);
  image(cell_state,col=cols,breaks=c(-1:7),xaxt='n',yaxt='n',main='Cells')
  leg_text<-c('0: vacant','1: uninf','2: inf_np','3: inf_p','4: dead')
  legend('topright',leg_text[types+1],fill=cols[types+1],bty='n',ncol=1)
  legend('bottomleft',paste(round(calc_plaque_size(settings),2),'mm'))
  
  #Panel 2: Virus
  if(sum(free_virus)==0){
    image(free_virus,col='grey95',xaxt='n',yaxt='n',main='Virus')
  }else{
    # image(free_virus,col=c('grey',colorRampPalette(c('blue','white'))(10)),xaxt='n',yaxt='n',main='Virus');
    image(free_virus,col=c('grey95',brewer.pal(9,'Blues')),xaxt='n',yaxt='n',main='Virus');
  }
  legend('bottomleft',paste(dynamic$virions))
  
  #Panel 3: T cells +cytokine
  cols<-settings$col_tcells;
  temp_grid<-t_cell_state;
  temp_grid[t_cell_state%in%c(1,2)&cytokine>settings$cyt_disp_thresh]<-4;
  temp_grid[t_cell_state%in%c(0,3)&cytokine>settings$cyt_disp_thresh]<-5;
  types<-as.numeric(as.data.frame(table(temp_grid),stringsAsFactors=FALSE)[,1]);
  image(temp_grid,col=cols,breaks=c(-1:6),xaxt='n',yaxt='n',main='T cells');
  leg_text<-c('Vacant','TRM','TEM','dead','Cyt+T','Cyt-T');
  legend('topright',leg_text[types+1],fill=cols[types+1],bty='n');
  
  #Panel 4: Cytokine
  if(sum(cytokine)==0){
    image(cytokine,col='grey95',xaxt='n',yaxt='n',main='Cytokine')
  }else{
    # image(cytokine,col=c('grey',colorRampPalette(c('Red','white'))(10)),xaxt='n',yaxt='n',main='Cytokine');
    image(cytokine*10,col=c('grey95',brewer.pal(9,'RdPu')),xaxt='n',yaxt='n',main='Cytokine');
  }
  
  #Timestamp
  mtext(paste('Time =',round(dynamic$time,digits=1),'h'),side=1,outer=T);
  
  if(settings$visualize=='gif'){dev.off()}
}

#COMPUTE GAUSSIAN kernel - only run once, at the beginning of simulation
compute_gaussian<-function(mask_range,sig=1.2,amp=1){
	temp<-(mask_range-1)/2;
	y<-matrix(rep(seq(-temp,temp),mask_range),mask_range); x<-t(y);
	gaussian_mask<-(amp/(2*pi*(sig^2)))*exp(-(x^2+y^2)/(2*sig^2));
	gaussian_mask<-gaussian_mask/sum(gaussian_mask); #this makes everything sum to 1
	return(gaussian_mask);
	rm(sig,y,x,temp); gc();
}

#DIFFUSION of free virus
diffusion_virus<-function(settings,i,j,gaussian_mask,parameters){
  if(free_virus[i,j]==1){
		nbr<-pick_nbr(settings,i,j);
		if(!is.null(nbr)){
			free_virus[i,j]<<-0;
			free_virus[nbr[1],nbr[2]]<<-free_virus[nbr[1],nbr[2]]+1;
		}
	}else{
		burst<-gaussburst(free_virus[i,j],gaussian_mask); #nxn array with num of phage in burst at each site
		free_virus[i,j]<<-0; #virus at focal site [i,j] will be spread across diffusion range
		distrib_progeny(settings,i,j,burst,parameters$diff_range,parameters); 
	}
}

#DIFFUSION of cytokine
diffusion_cyt<-function(settings,i,j,gaussian_mask,parameters){
  burst<-cytokine[i,j]*gaussian_mask;
  cytokine[i,j]<<-0;
  imin<-max(1,i-1);imax<-min(i+1,settings$L); #deal with edges
  jmin<-max(1,j-1);jmax<-min(j+1,settings$L);
  imin1<-1-(i-imin)+1;imax1<-1+(imax-i)+1;
  jmin1<-1-(j-jmin)+1;jmax1<-1+(jmax-j)+1;
  cytokine[imin:imax,jmin:jmax]<<-cytokine[imin:imax,jmin:jmax]+burst[imin1:imax1,jmin1:jmax1]; #boundary is "porous" - see notes
}


#GAUSSIAN BURST: generates a burst according to discrete approximation of a gaussian distribution
gaussburst<-function(burst_size,gaussian_mask){
	burst<-round(burst_size*gaussian_mask);
	if(sum(burst)>burst_size){ 
    exc<-sum(burst)-burst_size;
    tmp<-sample(seq(1,length(burst))[c(burst>0)],exc,replace=TRUE);
    counts<-as.data.frame(table(tmp),stringsAsFactors=F);
    burst[as.integer(counts$tmp)]<-burst[as.integer(counts$tmp)]-as.integer(counts$Freq);
    burst[burst<0]<-0;
	}else if(sum(burst)<burst_size){
		def<-burst_size-sum(burst);
		tmp<-sample(seq(1,length(burst)),def,replace=TRUE);
		counts<-as.data.frame(table(tmp),stringsAsFactors=F);
		burst[as.integer(counts$tmp)]<-burst[as.integer(counts$tmp)]+as.integer(counts$Freq);
	}
	return(burst);
}

#DISTRIBUTE PROGENY in free virus array for burst or diffusion
distrib_progeny<-function(settings,i,j,burst,dist_range,parameters){
	temp2<-(dist_range-1)/2;
	imin<-max(1,i-temp2);imax<-min(i+temp2,settings$L); #deal with edges
	jmin<-max(1,j-temp2);jmax<-min(j+temp2,settings$L);
	imin1<-temp2-(i-imin)+1;imax1<-temp2+(imax-i)+1;
	jmin1<-temp2-(j-jmin)+1;jmax1<-temp2+(jmax-j)+1;
	free_virus[imin:imax,jmin:jmax]<<-free_virus[imin:imax,jmin:jmax]+burst[imin1:imax1,jmin1:jmax1]; #boundary is "porous" - see notes
}

#PICK NEIGHBOUR from immediate moore neighbourhood
pick_nbr<-function(settings,i,j){
	delta<-sample(c(-1,0,1),2,replace=TRUE);
	nbr<-c(i,j)+delta;
	if(any(nbr<1|nbr>settings$L)||delta==0){
		return(NULL); #outside grid
	}else{
		return(nbr);
	}
}

#DETECT stuff within detection range
#If return_pos is true, function will return the position of the detected object
#If there are multiple objects detected, one will be picked at random and its position 
#returned
#If return_sum is true, function will return the total number of "what"s in the nbhd
detect_inf<-function(settings,det_range,i,j,what,return_pos=F){
  if(what=='none'){
    return(TRUE);
  }else{
    imin<-max(1,i-(det_range-1)/2);imax<-min(i+(det_range-1)/2,settings$L);
    jmin<-max(1,j-(det_range-1)/2);jmax<-min(j+(det_range-1)/2,settings$L);
    
    if(return_pos==F){		
      if(what=='any'){
        det_inf<-any(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2)|
                       cell_state[imin:imax,jmin:jmax]%in%c(2,3));
      }else if(what=='inf'){
        det_inf<-any(cell_state[imin:imax,jmin:jmax]%in%c(2,3));
      }else if(what=='infp'){
        det_inf<-any(cell_state[imin:imax,jmin:jmax]==3);
      }else if(what=='infn'){
        det_inf<-any(cell_state[imin:imax,jmin:jmax]==2);
      }else if(what=='all_cells'){
        det_inf<-any(cell_state[imin:imax,jmin:jmax]!=0);
      }else if(what=='trm'){
        det_inf<-any(t_cell_state[imin:imax,jmin:jmax]==1);
      }else if(what=='tem'){
        det_inf<-any(t_cell_state[imin:imax,jmin:jmax]==2);
      }else if(what=='t_cell'){
        det_inf<-any(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2));
      }else if(what=='uninf'){
        det_inf<-any(cell_state[imin:imax,jmin:jmax]==1);
      }else if(what=='cytokine'){
        det_inf<-any(cytokine[imin:imax,jmin:jmax]>0);
      }
      return(det_inf);
      
    }else if(return_pos==T){
      if(what=='any'){
        det_inf<-which(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2)|
                         cell_state[imin:imax,jmin:jmax]%in%c(2,3),arr.ind=T);
      }else if(what=='inf'){
        det_inf<-which(cell_state[imin:imax,jmin:jmax]%in%c(2,3),arr.ind=T);
      }else if(what=='infp'){
        det_inf<-which(cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
      }else if(what=='infn'){
        det_inf<-which(cell_state[imin:imax,jmin:jmax]==2,arr.ind=T);
      }else if(what=='all cells'){
        det_inf<-which(cell_state[imin:imax,jmin:jmax]!=0,arr.ind=T);
      }else if(what=='trm'){
        det_inf<-which(t_cell_state[imin:imax,jmin:jmax]==1,arr.ind=T);
      }else if(what=='tem'){
        det_inf<-which(t_cell_state[imin:imax,jmin:jmax]==2,arr.ind=T);
      }else if(what=='t_cell'){
        det_inf<-which(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2),arr.ind=T);
      }else if(what=='uninf'){
        det_inf<-which(cell_state[imin:imax,jmin:jmax]==1,arr.ind=T);
      }else if(what=='cytokine'){
        det_inf<-which(cytokine[imin:imax,jmin:jmax]>0,arr.ind=T);
      }
      
      if(length(det_inf)==0){ 
        return(NULL);
      }else{
        if(nrow(det_inf)>1) det_inf<-det_inf[sample(seq(1,nrow(det_inf)),1),];
        det_inf<-det_inf+c(imin,jmin)-1;
        return(det_inf);
      }		
    }else if(return_pos=='sum'){
      if(what=='any'){
        det_inf<-sum(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2)|
                       cell_state[imin:imax,jmin:jmax]%in%c(2,3));
      }else if(what=='inf'){
        det_inf<-sum(cell_state[imin:imax,jmin:jmax]%in%c(2,3));
      }else if(what=='infp'){
        det_inf<-sum(cell_state[imin:imax,jmin:jmax]==3);
      }else if(what=='infn'){
        det_inf<-sum(cell_state[imin:imax,jmin:jmax]==2);
      }else if(what=='all cells'){
        det_inf<-sum(cell_state[imin:imax,jmin:jmax]!=0);
      }else if(what=='trm'){
        det_inf<-sum(t_cell_state[imin:imax,jmin:jmax]==1);
      }else if(what=='tem'){
        det_inf<-sum(t_cell_state[imin:imax,jmin:jmax]==2);
      }else if(what=='t_cell'){
        det_inf<-sum(t_cell_state[imin:imax,jmin:jmax]%in%c(1,2));
      }else if(what=='uninf'){
        det_inf<-sum(cell_state[imin:imax,jmin:jmax]==1);
      }else if(what=='cytokine'){
        det_inf<-sum(cytokine[imin:imax,jmin:jmax]);
      }
      return(det_inf);
      
    }else{
      stop('Not a valid option for return_pos in detect_inf().')
    }	
  }
}

#T-CELL MOVEMENT 
#Based on settings$trm_move_method, load appropriate TRM patrolling function
set_trm_move_function<-function(trm_move_method){
  if(trm_move_method=='random'){
    #Random walk
    move_tcell<<-function(i,j,settings,parameters){
      delta<-sample(c(-1,0,1),2,replace=TRUE);
      if(all(delta==0)){ #same spot as before, so no movement
        newpos<-NULL;
      }else{
        newpos<-c(i,j)+delta;
      }	
      if(!is.null(newpos)){
        #Check if outside grid or if occupied. If so, no move
        if(any(newpos<1|newpos>settings$L)||!(t_cell_state[newpos[1],newpos[2]] %in% c(0,3))){
          newpos<-NULL; 
          
        #Make the move
        }else{
          t_cell_state[newpos[1],newpos[2]]<<-t_cell_state[i,j];
          t_cell_state[i,j]<<-0;
        }
      }
    }
  }else if(trm_move_method=='directed'){
    move_tcell<<-function(i,j,settings,parameters){
      det_cyt<-detect_inf(settings,parameters$tcell_detection_range,i,j,what='cytokine',return_pos=T);
      det_inf<-detect_inf(settings,parameters$tcell_detection_range,i,j,what='infp',return_pos=T);
        
      if(is.null(det_cyt)&is.null(det_inf)){ #random if there's no cytokine or infected cell in nbhd
        delta<-sample(c(-1,0,1),2,replace=TRUE);
        if(all(delta==0)){ #same spot as before, so no movement
          newpos<-NULL;
        }else{
          newpos<-c(i,j)+delta;
        }
      }else{ #else, move towards cytokine and antigen (preference to infected cell, then to cytokine gradient)
        if(is.null(det_inf)){
          del<-det_cyt-c(i,j);
        }else{
          del<-det_inf-c(i,j);
        }
        del[del>0]<-1;
        del[del<0]<--1;
        newpos<-unique(rbind(c(0,del[2]),c(del[1],0),del)); 
        if(nrow(newpos)>1) newpos<-newpos[sample(seq(1:nrow(newpos)),1),]; #pick new pos at random
        newpos<-newpos+c(i,j);
      }
      
      if(!is.null(newpos)){
        
        #Check if outside grid or if occupied. If so, no move
        if(any(newpos<1|newpos>settings$L)||!(t_cell_state[newpos[1],newpos[2]] %in% c(0,3))){
          newpos<-NULL; 
          
          #Make the move
        }else{
          t_cell_state[newpos[1],newpos[2]]<<-t_cell_state[i,j];
          t_cell_state[i,j]<<-0;
        }
      }
    }
  }else if(trm_move_method=='levy'){
    #(Levy flight)
    #Using implementation on stack exchange and assuming min jump=1 (see levy_tests.R, implementation 2b)
    #http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870
    move_tcell<<-function(i,j,settings,parameters){
      theta<-runif(1)*2*pi; #pick a random direction
      l<-(runif(1)^(-1/parameters$trm_levy_alpha)); #Levy displacement, assumes min=1
      l[l>parameters$trm_levy_max]<-parameters$trm_levy_max; #truncate at max
      newi<-min(max(1,round(i+l*sin(theta))),settings$L);
      newj<-min(max(1,round(j+l*cos(theta))),settings$L);
      if(t_cell_state[newi,newj]%in%c(0,3)){ #move if site is empty
        t_cell_state[newi,newj]<<-t_cell_state[i,j]; 
        t_cell_state[i,j]<<-0;
      }
    }
  }
}






#CELL DIVISION of T-cells
tcell_div<-function(i,j,settings){
	cell_nbrs<-3;
	imin1<-max(1,i-(cell_nbrs-1)/2);imax1<-min(i+(cell_nbrs-1)/2,settings$L);
	jmin1<-max(1,j-(cell_nbrs-1)/2);jmax1<-min(j+(cell_nbrs-1)/2,settings$L);
	vacant_div<-which(t_cell_state[imin1:imax1,jmin1:jmax1]==0|
					  t_cell_state[imin1:imax1,jmin1:jmax1]==3,arr.ind=TRUE);
	if(nrow(vacant_div)>0){
		ind<-vacant_div[sample(seq(1,nrow(vacant_div)),1),]; 
		m=imin1+ind[1]-1; n=jmin1+ind[2]-1;
		t_cell_state[m,n]<<-t_cell_state[i,j];
	}
}


#DECAY of free virus
viral_decay<-function(i,j,parameters,dice){
  if(free_virus[i,j]==1){
		if(dice<parameters$freevirus_decay/parameters$max_rate){ #stochastic loss
			free_virus[i,j]<<-0;
		}
	}else{
		free_virus[i,j]<<-round(free_virus[i,j]*exp(-parameters$freevirus_decay/parameters$max_rate));
	}
}

#DECAY of cytokine
cyt_decay<-function(i,j,parameters){
  cytokine[i,j]<<-max(0,cytokine[i,j]*exp(-parameters$cyt_decay_rate/parameters$max_rate));
}


# OUTPUT 
output_results<-function(settings,results,header){
    output_file<-paste(settings$dirname,'/',settings$out_fname,sep='');
    
    if (!nchar(output_file) || output_file == "") {
		stop('No input file');
    }
    if (header) {
		cat(file=output_file,paste(names(dynamic),collapse=','),'\n',append=F);
    } else {
		cat(file=output_file,paste(dynamic,collapse=','),'\n',append=T);
    }
}

# Compute totals of cells and virus
compute_totals<-function(dynamic,parameters){
	dynamic$time<<-dynamic$num_updates/parameters$max_rate;
	dynamic$virions<<-sum(sum(free_virus));
	dynamic$susceptible_cells<<-sum(cell_state==1);
	dynamic$infected_nonprodcells<<-sum(cell_state==2);
	dynamic$infected_prodcells<<-sum(cell_state==3);
	dynamic$dead_cells<<-sum(cell_state==4);
	dynamic$trm<<-sum(t_cell_state==1);
	dynamic$tem<<-sum(t_cell_state==2);
	dynamic$dead_tcells<<-sum(t_cell_state==3);
	dynamic$cytokine<<-sum(cytokine);
	dynamic$viral_cells<<-sum(free_virus>0);
	dynamic$cyto_cells<<-sum(cytokine>0);
	if (dynamic$virions > dynamic$max_vl) {
	    dynamic$max_vl <<-  dynamic$virions;
	    dynamic$max_vl_time <<-  dynamic$time;
	}
}

#END of simulation - save grids
end_simulation<-function(settings,parameters){
	if(settings$save_grids){
	  save(free_virus,file=paste(settings$dirname,'/free_virus',sep=''));
	  save(cell_state,file=paste(settings$dirname,'/cell_state',sep=''));
	  save(t_cell_state,file=paste(settings$dirname,'/t_cell_state',sep=''));
	  save(cytokine,file=paste(settings$dirname,'/cytokine',sep=''));
	}
  # save(parameters,file=paste(settings$dirname,'/sim_parameters',sep=''));
	if(settings$visualize=='pdf'){
	    dev.off();
	}
	print(paste("Highest log VL",ifelse(dynamic$max_vl>0,log10(dynamic$max_vl),0),"at t=",dynamic$max_vl_time));
	print(paste("End of simulation at t=",dynamic$time));
	#print('End of simulation. So long and thanks for all the fish...')
}
