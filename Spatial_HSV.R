# Spatial model of HSV in skin (2D)
# If running for the first time on mac, may need to install X11 for visualization
# Pavitra Roychoudhury
# Fred Hutchinson Cancer Research Center
# Latest update: 27-Jan-16

# INITIALIZATION: lines to run when this file is sourced
rm(list=ls()); 

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
	settings$sim_length<-5*24;											#Length of simulation in hours
	
	#Visualization
	settings$visualize<-'off';                         #Options: 'screen','pdf','off'
	settings$col_cells<-c('grey','papayawhip','orange','turquoise','peachpuff4',        #0=Vacant, 1=uninf, 2=inf_np, 3=inf_p, 4=dead,
	                      'peachpuff2','orangered','turquoise4');                     #5=uninf-prot, 6=inf_np-prot, 7=inf_p-prot 
	settings$col_tcells<-c('grey','green4','orange','deeppink','black');#Vacant, patrolling TRM, activated TRM, activated TEM, dead cell
	settings$refr_freq<-24;												#Refresh frequency (hrs)
	
	#Place virus
	#settings$place_virus<-list(type='single point',num=1);				#Start with single virus at the center of grid
	#settings$place_virus<-list(type='random',num=10);					#Start with num viruses, randomly scattered
	#settings$place_virus<-list(type='none');							#No free virus at start of simulation
	settings$place_virus<-list(type='infected producer',num=1);
	
	return(settings);
}

# 1b. Get default parameters 
get_default_parameters<-function(){
	parameters<-list(); 
	
	#Rates (per site or cell per day)
	parameters$max_rate<-24;						  #Sets time resolution, should be greater than or equal to sum of all rates for a given state
	
	#Uninfected cells
	parameters$uninf_cell_death<-0;	#Rate of uninfected cell death 
	parameters$infectivity_free<-0.1;				#Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
	parameters$infectivity_bound<-0;			#Infectivity (per cell per day)
	parameters$cell_div<-0.077; 					#Rate of generation of new susceptible cells (repop of empty sites)
	
	#Uninfected, cytokine-protected cells
	parameters$uninf_prot_cell_death<-0;	#Rate of uninfected cell death
	parameters$infectivity_prot_free<-0.1;			#Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
	parameters$infectivity_prot_bound<-0;			#Infectivity (per cell per day)
	
	#Infected cells
	parameters$inf_cell_death<-1.25;			#Rate of infected cell death due to infection
	parameters$viral_lag<-8;						  #Rate at which inf non-producer becomes viral producing cell (per cell per day)
	parameters$vir_prod_rate_mean<-100000;	#Rate at which an infected cell produces free virus (per cell per day, mean and sd)
	parameters$vir_prod_rate_sd<-10000;
	
	#Infected, cytokine-protected cells
	parameters$inf_prot_cell_death<-2.5;			#Rate of infected cell death due to infection
	parameters$viral_lag_prot<-8;						  #Rate at which inf non-producer becomes viral producing cell (per cell per day)
	parameters$vir_prod_rate_prot_mean<-1000;	#Rate at which an infected cell produces free virus (per cell per day, mean and sd)
	parameters$vir_prod_rate_prot_sd<-50;
	
	#Virus
	parameters$diff_rate<-12; 						  #Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
	parameters$freevirus_decay<-8.8;		#Rate of decay of free virus (per day)
	parameters$neur_drip<-0;						  #Rate of arrival of virus from neurons (per site per day)

	#T-cells
	parameters$tem_trafficking<-0.01;				#Rate of arrival of activated TEM from LN (per site per day)
	parameters$patrolling_trm_decay<-0;			#Death rate of patrolling TRMs (per site per day, default = 0)
	parameters$activated_trm_decay<-0.001;	#Death rate of activated TRMs (per site per day, default = 1E-3)
	parameters$trm_conv_rate<-0.01;					#Rate at which activated T-cell becomes resting in the absence of infection (per site per day, default = 1E-4)
	parameters$activated_tem_decay<-0.01;		#Death rate of activated TEMs (per site per day, default = same as act TRM)
	parameters$activated_trm_dbl<-1;				#Proliferation (doubling) rate of activated TRMs (per cell per day)
	parameters$activated_tem_dbl<-1;				#Proliferation (doubling) rate of activated TEMs (per cell per day)
	parameters$tcell_killing_rate<-96;      #Rate of target killing by T-cells (both aTRM and TEM)
	parameters$trm_motility<-5;					  #T-cell motility (jumps/movements per day)
	parameters$trm_levy_alpha<-2.5;       #Levy exponent (should be between 1 and 3)
	parameters$trm_levy_max<-10;          #Maximum jump per timestep (sites)

	#Cytokine-related parameters
	parameters$cyt_diff_rate<-50;         #Cytokine diffusion/spread rate in sites/day
	parameters$cyt_decay_rate<-0.1;       #Rate of clearance of cytokine (/site/day)
	parameters$cyt_activ<-48;             #Rate of cytokine-based activation given presence of cytokine at site (/cell/day)
	
	#Probabilities, fractions, etc
	parameters$host_density<-1;           #Fraction of sites occupied by susceptible cells at start	
	parameters$diff_range<-3;						  #n x n site, default n = 3, i.e. immediate Moore neighbourhood. Change this only if you want very rapid diffusion (more than max_rate sites per day) 
	parameters$tcell_density<-0.1;				#Density (fraction of total sites) of T-cells at the start of the simulation
	parameters$tcell_detection_range<-3;	#Spatial range of detection of antigen from infected cell (n x n), default n = 3
	parameters$tcell_killing_range<-3;		#Spatial range of killing (n x n), default n = 3
	
	#Other
	parameters$pat_trm_activatedby<-'inf';		#what activates a patrolling trm: 'inf','act_trm', 'act_tem' or 'any'
	parameters$tcell_killing_target<-'infp';	#What are the T-cells killing? 'infp','inf','all'
	
	return(parameters);
}



# 2. Initialize grid and other arrays
init_grid<-function(settings,parameters){
	print('Starting simulation.. hold on to your hat!');

	#Compute simulation length in number of updates
	settings$tot_updates<-settings$sim_length*parameters$max_rate/24;
	
	#Create arrays
	cell_state<<-array(data=as.integer(0),dim=c(settings$L,settings$L));
	t_cell_state<<-cell_state;
	free_virus<<-cell_state;											
	cytokine<<-array(data=as.logical(0),dim=c(settings$L,settings$L));
	
	#Directories
	print(paste('Results directory =',settings$dirname));
	dir.create(settings$dirname,showWarnings=FALSE);
	save(parameters,file=paste(settings$dirname,'/simulationparameters',sep=''));
	save(settings,file=paste(settings$dirname,'/simulationsettings',sep=''));
	
	#Dynamic variables - this stores time elapsed and population counts
	dynamic<<-list(num_updates=0,time=0,
				   susceptible_cells=0,infected_nonprodcells=0,infected_prodcells=0,dead_cells=0,virions=0,
				   patrolling_trm=0,activ_trm=0,activ_tem=0,dead_tcells=0);
	output_results(settings,dynamic,T); #write header
	
	#Visualization setup
	if(settings$visualize!='off'){
	  if(settings$visualize=='pdf'){
	    pdf(file=paste(settings$dirname,'/Results.pdf',sep=''),
	        width=8,height=8); 
	  }else if(settings$visualize=='screen'){
	    quartz(width=6,height=6);   
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
	next_screen_grab<-dynamic$num_updates+settings$refr_freq*parameters$max_rate/24;
	
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
			  ##--------------
			  #State 1: cytokine present
			  if(cytokine[i,j]){
			    
			    #Diffusion of cytokine
			    if(dice_cyt<=parameters$cyt_diff_rate/parameters$max_rate){
			      nbr<-pick_nbr(settings,i,j);
			      if(!is.null(nbr)){
			        cytokine[nbr[1],nbr[2]]<<-as.logical(1);
			      }
			     
			    #Cytokine decay 
			    }else if(dice_cyt<=(parameters$cyt_diff_rate+parameters$cyt_decay_rate)/parameters$max_rate){
			      cytokine[i,j]<<-as.logical(0);
			    }
			  }
			  
			  
				##T-cells - these are antigen-specific cells
				##--------------------------------------------
				#State 0 or 4: no T-cell present or occupied by dead T-cell
				if(t_cell_state[i,j]==0||t_cell_state[i,j]==4){
				
					#TEM arrives from LN (trafficking) if either inf cell, act TEM or act TRM in cytokine nbhd
					if(detect_inf(settings,parameters$tcell_detection_range,i,j,what='inf')&&	
						dice_t<=parameters$tem_trafficking/parameters$max_rate){
						t_cell_state[i,j]<<-3;
					}
			
			
				#State 1: resting/patrolling TRM - able to recognize infected cells but can't proliferate
				}else if(t_cell_state[i,j]==1){
					#decay/die
					if(dice_t<=parameters$patrolling_trm_decay/parameters$max_rate){
						t_cell_state[i,j]<<-4;
				
					#check neighbourhood for infected cells, if infected cell present in nbhd, change to state 2
					#or become activated due to cytokine
					}else if(detect_inf(settings,parameters$tcell_detection_range,i,j,
                              what=parameters$pat_trm_activatedby)||cytokine[i,j]){
						t_cell_state[i,j]<<-2;
				
					#move
					}else if(dice_t<=(parameters$patrolling_trm_decay+parameters$trm_motility)/parameters$max_rate){
					  move_tcell(i,j,settings,parameters);
					}
			
				#State 2: activated TRM - can proliferate, don't patrol
				}else if(t_cell_state[i,j]==2){
				  
				  #Release cytokine
				  cytokine[i,j]<<-as.logical(1);
				  
          #Case 1: Antigen present
          if(detect_inf(settings,parameters$tcell_detection_range,i,j,what=parameters$pat_trm_activatedby)){
            
            #proliferate
            if(dice_t<=parameters$activated_trm_dbl/parameters$max_rate){
              tcell_div(i,j,settings);
            }
				  
          #Case 2: Antigen absent
          }else{
            
            #die/decay/return to blood
            if(dice_t<=parameters$activated_trm_decay/parameters$max_rate){
              t_cell_state[i,j]<<-4;
              
            #go back to resting state (pTRM)
            }else if(dice_t<=(parameters$activated_trm_decay+parameters$trm_conv_rate)/parameters$max_rate){
              t_cell_state[i,j]<<-1;
            }
            
          }
  			
				#State 3: activated TEM (need to have trafficked from LN) -- same functions as aTRM
				}else if(t_cell_state[i,j]==3){
					
				  #Release cytokine
				  cytokine[i,j]<<-as.logical(1);
				  
				  #Case 1: Antigen present
				  if(detect_inf(settings,parameters$tcell_detection_range,i,j,what=parameters$pat_trm_activatedby)){
				    
				    #proliferate
				    if(dice_t<=parameters$activated_tem_dbl/parameters$max_rate){
				      tcell_div(i,j,settings);
				    }
				    
				  #Case 2: Antigen absent
				  }else{
				    
				    #die/decay/return to blood
				    if(dice_t<=parameters$activated_tem_decay/parameters$max_rate){
				      t_cell_state[i,j]<<-4;
				      
				      #go back to resting state (pTRM)
				    }else if(dice_t<=(parameters$activated_tem_decay+parameters$trm_conv_rate)/parameters$max_rate){
				      t_cell_state[i,j]<<-1;
				    }				    
				  }									
				}
							
				##Host cells
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
					if(free_virus[i,j]>0&&dice_cells<=(parameters$neur_drip+parameters$infectivity_free*free_virus[i,j])/parameters$max_rate){
						cell_state[i,j]<<-2; #becomes an infected cell
						free_virus[i,j]<<-free_virus[i,j]-1;						
					
					#Event 1.2: cell death
					}else if(dice_cells<=(parameters$neur_drip+parameters$uninf_cell_death)/parameters$max_rate){
						cell_state[i,j]<<-4; #death of uninfected cell
					
					#Event 1.3: cytokine protection
					}else if(cytokine[i,j]&&
					         dice_cells<=(parameters$neur_drip+parameters$uninf_cell_death+parameters$cyt_activ)/parameters$max_rate){
					  cell_state[i,j]<<-5;
					}
			
				#State 2: infected cell (non-virus producing)
				}else if(cell_state[i,j]==2){					
								
					#Event 2.1: become a virus-producing cell
					if(dice_cells<=(parameters$neur_drip+parameters$viral_lag)/parameters$max_rate){
						cell_state[i,j]<<-3; #becomes virus-producing infected cell
						
					#Event 2.2: cell death
					}else if(dice_cells<=(parameters$neur_drip+parameters$viral_lag+parameters$uninf_cell_death)/parameters$max_rate){
						cell_state[i,j]<<-4; #death of infected cell			
					
					#Event 2.3: cytokine protection
					}else if(cytokine[i,j]&&
					         dice_cells<=(parameters$neur_drip+parameters$viral_lag+parameters$uninf_cell_death+
					                      parameters$cyt_activ)/parameters$max_rate){
					  cell_state[i,j]<<-6; 
					}
				  
				#State 3: infected cell (producing virus)
				}else if(cell_state[i,j]==3){
					
					#Event 3.0: produce free virus 
					free_virus[i,j]<<-free_virus[i,j]+max(0,round(rnorm(1,mean=parameters$vir_prod_rate_mean/parameters$max_rate,
																	sd=parameters$vir_prod_rate_sd/parameters$max_rate)));
				
					#Event 3.1: infect neighbouring cell
					nbr<-pick_nbr(settings,i,j);
					if(!is.null(nbr)&&cell_state[nbr[1],nbr[2]]==1&&
									  dice_cells<=(parameters$neur_drip+parameters$infectivity_bound)/parameters$max_rate){
						cell_state[nbr[1],nbr[2]]<<-2;	
				
					#Event 3.2: cell death due to virus
					}else if(dice_cells<=(parameters$neur_drip+parameters$inf_cell_death)/parameters$max_rate){
						cell_state[i,j]<<-4; #death of uninfected cell			
					
					#Event 3.3: cell death due to T-cell
					}else if(detect_inf(settings,parameters$tcell_detection_range,i,j,what='act_t_cell')&&
					         dice_cells<=(parameters$neur_drip+parameters$inf_cell_death+parameters$tcell_killing_rate)/parameters$max_rate){
					  cell_state[i,j]<<-4; #death of infected cell
					
					#Event 3.4: cytokine protection
					}else if(cytokine[i,j]&&
					         dice_cells<=(parameters$neur_drip+parameters$inf_cell_death+parameters$tcell_killing_rate+
					                      parameters$cyt_activ)/parameters$max_rate){
					  cell_state[i,j]<<-7;
					}
				
				#State 5: uninfected, protected cell 
				}else if(cell_state[i,j]==5){
				  
				  #Event 5.1: infection due to free virus at the site
				  if(free_virus[i,j]>0&&dice_cells<=(parameters$neur_drip+parameters$infectivity_prot_free*free_virus[i,j])/parameters$max_rate){
				    cell_state[i,j]<<-6; #becomes an infected cell
				    free_virus[i,j]<<-free_virus[i,j]-1;						
				    
				    #Event 5.2: cell death
				  }else if(dice_cells<=(parameters$neur_drip+parameters$uninf_prot_cell_death)/parameters$max_rate){
				    cell_state[i,j]<<-4; #death of uninfected cell
				  }
				
				#State 6: infected cell (non-virus producing), protected
				}else if(cell_state[i,j]==6){					
				  
				  #Event 6.1: become a virus-producing (cytokine-protected) cell
				  if(dice_cells<=(parameters$neur_drip+parameters$viral_lag_prot)/parameters$max_rate){
				    cell_state[i,j]<<-7; #becomes virus-producing infected cell
				    
				    #Event 6.2: cell death
				  }else if(dice_cells<=(parameters$neur_drip+parameters$viral_lag_prot+parameters$uninf_prot_cell_death)/parameters$max_rate){
				    cell_state[i,j]<<-4; #death of infected cell			
				    
				  }
				  
				#State 7: infected cell (producing virus), protected
				}else if(cell_state[i,j]==7){
				  
				  #Event 7.0: produce free virus 
				  free_virus[i,j]<<-free_virus[i,j]+max(0,round(rnorm(1,mean=parameters$vir_prod_rate_prot_mean/parameters$max_rate,
				                                                      sd=parameters$vir_prod_rate_prot_sd/parameters$max_rate)));
				  
				  #Event 7.1: infect neighbouring cell
				  nbr<-pick_nbr(settings,i,j);
				  if(!is.null(nbr)&&cell_state[nbr[1],nbr[2]]==1&&
				     dice_cells<=(parameters$neur_drip+parameters$infectivity_prot_bound)/parameters$max_rate){
				    cell_state[nbr[1],nbr[2]]<<-2;
				  }else if(!is.null(nbr)&&cell_state[nbr[1],nbr[2]]==5&&
				      dice_cells<=(parameters$neur_drip+parameters$infectivity_prot_bound)/parameters$max_rate){
				    cell_state[nbr[1],nbr[2]]<<-6;
				    
				  #Event 7.2: cell death due to virus
				  }else if(dice_cells<=(parameters$neur_drip+parameters$inf_prot_cell_death)/parameters$max_rate){
				    cell_state[i,j]<<-4; #death of uninfected cell			
				    
				  #Event 7.3: cell death due to T-cell
				  }else if(detect_inf(settings,parameters$tcell_detection_range,i,j,what='act_t_cell')&&
				           dice_cells<=(parameters$neur_drip+parameters$inf_prot_cell_death+parameters$tcell_killing_rate)/parameters$max_rate){
				    cell_state[i,j]<<-4; #death of infected cell
				  }
				}
			
				##Cell-free virus
				##--------------------------------------------
				if(free_virus[i,j]>0){

					#Event 0.1: Decay
					viral_decay(i,j,parameters,dice_v);
				
					#Event 0.2: Diffusion of remaining free virus at site
					if(dice_v<=parameters$diff_rate/parameters$max_rate){
						diffusion(settings,i,j,gaussian_mask_diffusion,parameters);}			
				
				}
			}			
			
		}
		
		dynamic$num_updates<<-dynamic$num_updates+1;
		compute_totals(dynamic,parameters);
		
		#print(dynamic$num_updates);
		if(settings$visualize!='off'&&dynamic$num_updates==next_screen_grab){
			visualize(settings);
			next_screen_grab<-dynamic$num_updates+settings$refr_freq*parameters$max_rate/24;
		}
		output_results(settings,results,F);
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

#PLACE T-CELLS
place_tcells<-function(settings,parameters){
	if(parameters$tcell_density>0){
		locs<-sample(seq(1,length(t_cell_state)),length(t_cell_state),replace=FALSE);
		locs<-locs[runif(length(locs),min=0,max=1)<parameters$tcell_density];
		t_cell_state[locs]<<-1; rm(locs);
	}
}

#NEURONAL DRIP
neur_drip<-function(i,j){
	free_virus[i,j]<<-free_virus[i,j]+1;
}

#VISUALIZE
visualize<-function(settings){
  cols<-settings$col_cells;
  types<-as.numeric(as.data.frame(table(cell_state),stringsAsFactors=FALSE)[,1]);
  image(cell_state,col=cols,breaks=c(-1:7),xaxt='n',yaxt='n',main='Cells')
  leg_text<-c('0: vacant','1: uninf','2: inf_np','3: inf_p','4: dead',
              '5: uninf_prot','6: inf_np-prot','7: inf_p-prot')
  legend('topright',leg_text[types+1],fill=cols[types+1],bty='n',ncol=2)
  image(free_virus>0,col=c('grey','blue'),xaxt='n',yaxt='n',main='Free virus');
  legend('topright',c('virus absent','virus present'),fill=c('grey','blue'), bty='n')
  cols<-settings$col_tcells;
  types<-as.numeric(as.data.frame(table(t_cell_state),stringsAsFactors=FALSE)[,1]);
  image(t_cell_state,col=cols,breaks=c(-1:4),xaxt='n',yaxt='n',main='T-cells')
  leg_text<-c('0 = vacant','1 = patrolling TRM','2 = activated TRM',
              '3 = activated TEM', '4 = dead T-cell');
  legend('topright',leg_text[types+1],fill=cols[types+1],bty='n',ncol=2)
  image(cytokine==1,col=c('grey','red'),xaxt='n',yaxt='n',main='Cytokine');
  legend('topright',c('cytokine absent','cytokine present'),fill=c('grey','red'), bty='n')
  mtext(paste('Time =',dynamic$time,'h'),side=1,outer=T)
}

#COMPUTE GAUSSIAN kernel - only run once, at the beginning of simulation
compute_gaussian<-function(mask_range){
	sig<-1.2; temp<-(mask_range-1)/2;
	y<-matrix(rep(seq(-temp,temp),mask_range),mask_range); x<-t(y);
	gaussian_mask<-(1/(2*pi*(sig^2)))*exp(-(x^2+y^2)/(2*sig^2));
	gaussian_mask<-gaussian_mask/sum(gaussian_mask); #this makes everything sum to 1
	return(gaussian_mask);
	rm(sig,y,x,temp); gc();
}

#DIFFUSION of free virus
diffusion<-function(settings,i,j,gaussian_mask,parameters){
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

#DETECT INFECTED CELL or activated T-cell within detection range
#If return_pos is true, function will return the position of the detected object
#If there are multiple objects detected, one will be picked at random and its position 
#returned
detect_inf<-function(settings,det_range,i,j,what,return_pos=F){
	if(what=='none'){
		return(TRUE);
	}else{
		imin<-max(1,i-(det_range-1)/2);imax<-min(i+(det_range-1)/2,settings$L);
		jmin<-max(1,j-(det_range-1)/2);jmax<-min(j+(det_range-1)/2,settings$L);
		
		if(return_pos==F){		
			if(what=='any'){
				det_inf<-any(t_cell_state[imin:imax,jmin:jmax]%in%c(2,3)|
				               cell_state[imin:imax,jmin:jmax]%in%c(2,3));
			}else if(what=='inf'){
				det_inf<-any(cell_state[imin:imax,jmin:jmax]==2|
							 cell_state[imin:imax,jmin:jmax]==3);
			}else if(what=='infp'){
				det_inf<-any(cell_state[imin:imax,jmin:jmax]==3);
			}else if(what=='infn'){
				det_inf<-any(cell_state[imin:imax,jmin:jmax]==2);
			}else if(what=='all_cells'){
			  det_inf<-any(cell_state[imin:imax,jmin:jmax]!=0);
			}else if(what=='act_trm'){
				det_inf<-any(t_cell_state[imin:imax,jmin:jmax]==2);
			}else if(what=='act_tem'){
				det_inf<-any(t_cell_state[imin:imax,jmin:jmax]==3);
			}else if(what=='act_t_cell'){
			  det_inf<-any(t_cell_state[imin:imax,jmin:jmax]%in%c(2,3));
			}else if(what=='uninf'){
			  det_inf<-any(cell_state[imin:imax,jmin:jmax]==1);
			}
			return(det_inf);
		
		}else if(return_pos==T){
			if(what=='any'){
				det_inf<-which(cell_state[imin:imax,jmin:jmax]==2|
							   cell_state[imin:imax,jmin:jmax]==3|
							   t_cell_state[imin:imax,jmin:jmax]==2|
							   t_cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
			}else if(what=='inf'){
				det_inf<-which(cell_state[imin:imax,jmin:jmax]==2|
							   cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
			}else if(what=='infp'){
				det_inf<-which(cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
			}else if(what=='infn'){
				det_inf<-which(cell_state[imin:imax,jmin:jmax]==2,arr.ind=T);
			}else if(what=='all cells'){
			  det_inf<-which(cell_state[imin:imax,jmin:jmax]!=0,arr.ind=T);
			}else if(what=='act_trm'){
				det_inf<-which(t_cell_state[imin:imax,jmin:jmax]==2,arr.ind=T);
			}else if(what=='act_tem'){
				det_inf<-which(t_cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
			}else if(what=='act_t_cell'){
			  det_inf<-which(t_cell_state[imin:imax,jmin:jmax]==2|
			                 t_cell_state[imin:imax,jmin:jmax]==3,arr.ind=T);
			}else if(what=='uninf'){
			  det_inf<-which(cell_state[imin:imax,jmin:jmax]==1,arr.ind=T);
			}
			
			if(length(det_inf)==0){ 
				return(NULL);
			}else{
				if(nrow(det_inf)>1) det_inf<-det_inf[sample(seq(1,nrow(det_inf)),1),];
				det_inf<-det_inf+c(imin,jmin)-1;
				return(det_inf);
			}		
		}else{
			stop('Not a valid option for return_pos in detect_inf().')
		}	
	}
}

#T-CELL MOVEMENT (Levy flight)
#Using implementation on stack exchange and assuming min jump=1 (see levy_tests.R, implementation 2b)
#http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870
move_tcell<-function(i,j,settings,parameters){
  theta<-runif(1)*2*pi; #pick a random direction
  l<-(runif(1)^(-1/parameters$trm_levy_alpha)); #Levy displacement, assumes min=1
  l[l>parameters$trm_levy_max]<-parameters$trm_levy_max; #truncate at max
  newi<-min(max(1,round(i+l*sin(theta))),settings$L);
  newj<-min(max(1,round(j+l*cos(theta))),settings$L);
  if(t_cell_state[newi,newj]%in%c(0,4)){ #move if site is empty
    t_cell_state[newi,newj]<<-t_cell_state[i,j]; 
    t_cell_state[i,j]<<-0;
  }
}


#CELL DIVISION of T-cells
tcell_div<-function(i,j,settings){
	cell_nbrs<-3;
	imin1<-max(1,i-(cell_nbrs-1)/2);imax1<-min(i+(cell_nbrs-1)/2,settings$L);
	jmin1<-max(1,j-(cell_nbrs-1)/2);jmax1<-min(j+(cell_nbrs-1)/2,settings$L);
	vacant_div<-which(t_cell_state[imin1:imax1,jmin1:jmax1]==0|
					  t_cell_state[imin1:imax1,jmin1:jmax1]==4,arr.ind=TRUE);
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
	dynamic$time<<-dynamic$num_updates*24/parameters$max_rate;
	dynamic$virions<<-sum(sum(free_virus));
	dynamic$susceptible_cells<<-sum(cell_state==1);
	dynamic$infected_nonprodcells<<-sum(cell_state==2);
	dynamic$infected_prodcells<<-sum(cell_state==3);
	dynamic$dead_cells<<-sum(cell_state==4);
	dynamic$patrolling_trm<<-sum(t_cell_state==1);
	dynamic$activ_trm<<-sum(t_cell_state==2);
	dynamic$activ_tem<<-sum(t_cell_state==3);
	dynamic$dead_tcells<<-sum(t_cell_state==4);
}

#END of simulation - save grids
end_simulation<-function(settings,parameters){
	save(free_virus,file=paste(settings$dirname,'/free_virus',sep=''));
	save(cell_state,file=paste(settings$dirname,'/cell_state',sep=''));
	save(t_cell_state,file=paste(settings$dirname,'/t_cell_state',sep=''));
	write.csv(t(as.data.frame(parameters)),file=paste(settings$dirname,'/simulationparameters.csv',sep=''));
	
	if(settings$visualize=='pdf'){
	    dev.off();
	}
	print('End of simulation. So long and thanks for all the fish...')
}


