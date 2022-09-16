//NOTE: This read-only file is generated from the file ../include/full_settings.h.
//Any editing should be done in that file.

#ifndef SETTINGS_H
#define SETTINGS_H

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define NUM_STRAINS 3 

struct global_settings {
	gsl_rng * ur;		// Ramdom number generator

	FILE *dataF1;

	const char *dirname;   								//Directory to store results
	const char *out_fname;							//File to store cell and virus counts

	int verbose;	// extra messages

	int L;							//size of matrix (LxL)
	int Runs;						//number of patients or episodes to simulate
	int sim_length;						//Length of simulation in hours

	int save_plot_data;     //Save state data for R-style plots
	int save_ratio_data;    //Save running e/t data for R-style plots
	int save_state_files;   //Save unning tate files during R-style plotting
	int stop_early;		    // Stop sim when no more virus or infected cells
	const char *process_cmd;	// Rscript to run for state plots
	
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
	int max_start_strains;		// number of strains to seed (either virions or infected cells)
	int start_neur_strain;		// strains to begin seeding from neuron (either virions or infected cells)
	int max_neur_drips;		// number of drips to allow (virus or infected cells)
	

};

struct global_parameters {
	//Rates (per site or cell per day)
	int max_rate;						  //Sets time resolution, should be greater than or equal to sum of all rates for a given state
	double host_density;           //Fraction of sites occupied by susceptible cells at start	
	
	//Uninfected cells
	double uninf_cell_death;	//Rate of uninfected cell death 
	double infectivity_free[NUM_STRAINS];				//Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
	double cell_div; 					//Rate of generation of new susceptible cells (repop of empty sites)
	
	//Infected cells
	double inf_cell_death[NUM_STRAINS];			//Rate of infected cell death due to infection
	double viral_lag;						  //Rate at which inf non-producer becomes viral producing cell (per cell per day)
	int vir_prod_rate_mean;	//Rate at which an infected cell produces free virus (per cell per day, mean)
	
	//Virus
	double diff_rate_virus; 						  //Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
	int diff_range;						  //n x n site, default n = 3, i.e. immediate Moore neighbourhood. 
	double freevirus_decay;		//Rate of decay of free virus (per day)
	double neur_drip;						  //Rate of arrival of virus from neurons (per site per day)
	double neur_infect;						  //Rate of cellular infection from neurons (per site per day)
	double psi;	//fraction of virus picked up by swab

	double prob_recomb;	// probability of co-infection resulting in recombinant virus
	double prob_coinfect;	// probability of co-infection


	
};

struct global_dynamics {
	int runnum;	//number of current run (patient)
	double time;
	int virions;
	int virionsA;
	int virionsB;
	int virionsAB;
	int susceptible_cells;
	int infected_nonprodcells;
	int infn_A;
	int infn_B;
	int infn_ABn;
	int infn_ABr;
	int infected_prodcells;
	int infp_A;
	int infp_B;
	int infp_ABn;
	int infp_ABr;
	int dead_cells;
	int pat_hsv;
	int act_hsv;
	int pat_byst;
	int act_byst;
	int dead_tcells;

	double max_vl_time;
	int max_vl;
	int viral_cells;

	double first_samp_vl_time;
	int first_samp_vl;

	double max_samp_vl_time;
	int max_samp_vl;

	int num_updates;
	int tot_updates;

	int tot_infs;
	unsigned long int tot_virons;

	int *visited;

	int **cell_state;
	int ***viral_matrix;
	int ***delta_viral_matrix;


};
#endif
