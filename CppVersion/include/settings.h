#ifndef SETTINGS_H
#define SETTINGS_H

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

struct global_settings {
	gsl_rng * ur;		// Ramdom number generator

	FILE *dataF1;

	const char *dirname;   								//Directory to store results
	const char *out_fname;							//File to store cell and virus counts
	const char *et_file;							//File to input E/T ratios for T-cell init

	int verbose;	// extra messages

	int L;							//size of matrix (LxL)
	int Runs;						//number of patients or episodes to simulate
	int sim_length;						//Length of simulation in hours

	int panels;     //Number of panels to draw when plotting w/ R (passed arg)
	int save_plot_data;     //Save state data for R-style plots
	int save_ratio_data;    //Save running e/t data for R-style plots
	int save_state_files;   //Save unning tate files during R-style plotting
	int stop_early;		    // Stop sim when no more virus or infected cells
	int cyto_act;	   // If 1 then cytokines activate patrolling Trms

	int cyto_act_cyto;   // If 1 then cytokines activate cyto production 
	int cyto_act_prolif;   // If 1 then cytokines activate Trm prolif

	const char *process_cmd;	// Rscript to run for state plots
	const char *et_process_cmd;	// Rscript to run for E/T ratio plots
	
	double site_scale; //how many microns is one site
	double refr_freq; //how often to dump data for screen prints, etc.
	double data_freq; //how often to dump data to results file (sampling freq)
	
	//Place virus
	//settings->place_virus=list(type='single point',num=1);				//Start with single virus at the center of grid
	//settings->place_virus=list(type='random',num=10);					//Start with num viruses, randomly scattered
	//settings->place_virus=list(type='none');							//No free virus at start of simulation
	int place_virus_type;
	int place_virus_num;
	int start_plaques;		// number of cells either infected or with virus at the start
	
	//TRM 
	int trm_init; //random or gaussian
	int trm_move_method; //random', 'directed' or 'levy'
	int trm_random_gaussian;	// place center of gaussian dist randomly on grid
	int add_cd4s;  // are we doing just CD8s or adding in CD4s (BYST)
	double cd4_to_cd8_ratio; // this is 1:3 from histologies

	double et_ratio_p_max;
	double et_ratio_p_min;
	double et_ratio_p_mean;
	double et_ratio_p_std;
	double et_ratio_p_k;
	double et_ratio_p_a;
	double et_ratio_p_b;
	double et_ratio_p_c;
	double et_ratio_p_d;
	double et_ratio_c_max;
	double et_ratio_c_min;
	double et_ratio_c_mean;
	double et_ratio_c_std;
	double et_ratio_c_k;
	double et_ratio_c_a;
	double et_ratio_c_b;
	int et_ratio_function;
	int et_cluster_size;
	int et_num_clusters;
	int et_num_samples;

	int plot_cytokines;
};

struct global_parameters {
	//Rates (per site or cell per day)
	int max_rate;						  //Sets time resolution, should be greater than or equal to sum of all rates for a given state
	
	//Uninfected cells
	double uninf_cell_death;	//Rate of uninfected cell death 
	double infectivity_free;				//Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
	int max_tcell_div; 					//Max number of t cell divisions
	double cell_div; 					//Rate of generation of new susceptible cells (repop of empty sites)
	
	//Infected cells
	double inf_cell_death;			//Rate of infected cell death due to infection
	double viral_lag;						  //Rate at which inf non-producer becomes viral producing cell (per cell per day)
	int vir_prod_rate_mean;	//Rate at which an infected cell produces free virus (per cell per day, mean)
	
	//Virus
	double diff_rate_virus; 						  //Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
	double freevirus_decay;		//Rate of decay of free virus (per day)
	double neur_drip;						  //Rate of arrival of virus from neurons (per site per day)
	double psi;	//fraction of virus picked up by swab

	//T-cells
	double hsv_fract;	//fraction of Trms that are hsv specific (vs byst)
	double tem_trafficking;	//Rate of arrival of TEM from LN (per site per day)
	double tem_exiting;	//Rate of exit of TEM from tissue (per site per day in absense of cytos or antigen)
	double trm_decay;	//Death rate of TRMs and TEMs in absence of Ag (per site per day, default = 1E-3)
	double trm_conv_rate;					//Rate at which TEM becomes TRM in the absence of infection (per site per day, default = 1E-4)
	double trm_dbl;				//Proliferation (doubling) rate of TRMs/TEMs (per cell per day)
	double tem_delay;      //Delay in the arrival of TEM (per day, default = 1 = 24h delay)
	double tcell_killing_rate;      //Rate of target killing by T-cells (both aTRM and TEM)
	double trm_motility;					  //T-cell motility (jumps/movements per day)
	double trm_levy_alpha;       //Levy exponent (should be between 1 and 3)
	double trm_levy_max;          //Maximum jump per timestep (sites)

	double cyt_absorption_rate_trm;  //Rate of cytokine absorption by T-cells (pg/cell/day)
	double cyt_absorption_rate_cell;  //Rate of cytokine absorption by epithelial cells (pg/cell/day)

	//Cytokine-related parameters
	double max_cyt_days;     //days a t-cell or epithelial can secrete cytokine
	double cyt_secretion_rate_tcell_mean;     //Rate of cytokine secretion by T-cells, mean (pg/cell/day)
	double cyt_secretion_rate_infcell_mean;     //Rate of cytokine secretion by infected cells, mean (pg/cell/day)

	double cyt_diff_rate;      //Cytokine diffusion/spread rate in sites/day
	double cyt_uptake;      // Rate of cellular cytokine uptake 
	double cyt_decay_rate;       //Rate of clearance of cytokine (/day)

	double cyt_trm_ic50;
	double beta_ic50;
	double inf_ic50;
	double prod_ic50;

	double cyt_hill;
	double cyt_prot_expiration;

	int cyt_effect_by_type;
	bool cyt_effect_infectivity; //Do cytokines affect infectivity: yes (1) or no (0)
	bool cyt_effect_virprod; //Do cytokines affect virus production rate: yes (1) or no (0)
	double cyt_effect_lag; //Max effect of cytokines on lag (multiplier, remember that lag here is a rate)
	double cyt_effect_uninfcelldeath; //Max death rate of uninfected cells in the presence of cytokine, 0 means no effect
	double cyt_effect_infncelldeath; //Max death rate of pre-productively infected cells cells in the presence of cytokine, 0 means no effect
	double cyt_effect_infpcelldeath; //Max death rate of infected cells in the presence of cytokine, 0 means no effect
	
	//Probabilities, fractions, etc
	double host_density;           //Fraction of sites occupied by susceptible cells at start	
	int diff_range;						  //n x n site, default n = 3, i.e. immediate Moore neighbourhood. 
	double tcell_density;				//Density (fraction of total sites) of T-cells at the start of the simulation
	int tcell_detection_range;	//Spatial range of detection of antigen from infected cell (n x n), default n = 3

	//Other
	int pat_trm_activatedby;		//what activates a T-cell: 'inf','trm', 'act_tem' or 'any'
	int tcell_killing_target;	//What are the T-cells killing? 'infp','inf','all'. 'inf' means both early and late
	int limit_cyt_spread;	// limit cytokine diffusion to this dist from plaque
	
};

struct global_dynamics {
	int runnum;	//number of current run (patient)
	double time;
	int virions;
	int susceptible_cells;
	int infected_nonprodcells;
	int infected_prodcells;
	int dead_cells;
	int pat_hsv;
	int act_hsv;
	int pat_byst;
	int act_byst;
	int dead_tcells;
	double cytokine;

	double max_vl_time;
	int max_vl;
	int viral_cells;
	int cyto_cells;

	double first_samp_vl_time;
	int first_samp_vl;

	double max_samp_vl_time;
	int max_samp_vl;

	double gauss_offset;
	double avg_et_ratio;
	double et_ratio_3;
	double et_ratio_5;
	double et_ratio_7;
	double et_ratio_9;
	double et_ratio_11;
	double perc_zeros;
        double closestHSV;
        double closestTcell;

	double trm_act_time;
	double trm_cyt_time;

	int num_updates;
	int tot_updates;

	int tot_infs;
	unsigned long int tot_virons;

	int hsv_tems;
	int byst_tems;
	int tem_exits;

	int cyto_kills;
	int trm_kills;

	int *visited;

	int **cell_state;
	int **t_cell_state;

	int **t_cell_potential; // how many times can this t-cell further divide
	double **t_cell_cyto_dur;  // how many more days can this cell secrete cytokines

	int **viral_matrix;
	double **cytokine_matrix;
	
	int **delta_viral_matrix;
	double **delta_cytokine_matrix;
};
#endif
