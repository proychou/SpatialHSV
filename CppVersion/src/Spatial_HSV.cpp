//NOTE: This read-only file is generated from the file ../src/full_spatial.cpp.
//Any editing should be done in that file.

/************
// Spatial model of HSV in skin (2D)
// If running for the first time on mac, may need to install X11 for visualization
// Pavitra Roychoudhury (translated to c++ by D. Swan)
// Fred Hutchinson Cancer Research Center
// Translated : 6-March-18
***************/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<math.h>
#include<signal.h>
#include<stdarg.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <map>

#ifdef __sun
#include <strings.h>
#else
#include <string.h>
#endif
#include <stdlib.h>
// Some STL Headers
using namespace std;
 
// Some GSL Headers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "Spatial_HSV.h"


#define FRACT_W_ZEROS 0.333
#define DETECTION_THRESHOLD 100

// Types of distributions employed 
#define NONE 0
#define WEIBULL 1
#define NORMAL 2
#define LOGNORMAL 3
#define POLYNOMIAL 4
#define POWER_LAW 5
#define EXPONENTIAL 6

#define NO_VIRUS 0
#define SINGLE_POINT 1
#define RANDOM_VIRUS 2
#define INFECTED_PRODUCER 3

// cell states
#define EMPTY 0
#define SUSCEPTIBLE 1
#define UNPROTECTABLE 2
#define INFN_A 3
#define INFN_B 4
#define INFN_ABn 5
#define INFN_ABr 6
#define INFP_A 7
#define INFP_B 8
#define INFP_ABn 9
#define INFP_ABr 10
#define DEAD 11
#define CYTO_KILLED 12
#define TRM_KILLED 13

// categories for checking neighbors
#define ANY 1
#define INF 2
#define INFP 3
#define INFN 4
#define ALL_CELLS 5
#define PAT_HSV 6
#define ACT_HSV 7
#define PAT_BYST 8
#define ACT_BYST 9
#define T_CELL 10
#define T_KILL 11
#define UNINF 12
#define CYTOKINE 13
#define VIRUS 14

// methods of canvasing "neighborhood"
#define CHECK 0
#define COORDS 1
#define SUM 2

#define MAX_LINE 200
#define MAX_ET_ENTRIES 2000
#define MAX_CLUSTERS 100

// Macros for looking for various input file variables
#define CHECK_FOR_SETTING(name,param) \
    if (inputs.find(name) != inputs.end()) { settings->param = inputs[name]; }

#define CHECK_FOR_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { parameters->param = inputs[name]; }

#define CHECK_FOR_INT_SETTING(name,param) \
    if (inputs.find(name) != inputs.end()) { settings->param = static_cast<int>(inputs[name]); }

#define CHECK_FOR_INT_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { parameters->param = static_cast<int>(inputs[name]); } 

// special malloc to inform user of failures
void *my_malloc(int size)
{
    void *ret_val;
    ret_val=malloc(size);
    if (ret_val == NULL)
    {
	fprintf(stderr,"Encountered a problem allocating memory of size=%d. Exiting.\n",size);
	exit(1);
    }
    return ret_val;
}

// Routine to read and process the input file
void read_input_file(char *inp_file, global_settings *settings, global_parameters *parameters)
{
    ////////////////////////////////////////////////////////////////////////
    ///// read input parameters through the input file 
    char tmpline[MAX_LINE];
    char *valuep;
    int i=0;
    FILE *inf;

    map<string,double> inputs;
    string par;
    double parv;

    if ((inf = fopen (inp_file,"r")) == NULL) {
	cerr << "Could not open input file "<<inp_file<<"\n";
	exit(1);
    }
    while (fgets(tmpline, MAX_LINE-1,inf) != NULL) {
	i++;
	tmpline[MAX_LINE-1] = '\0';
	if (tmpline[0]=='#' || tmpline[0]=='\r' || tmpline[0]=='\n')
	{
	    cout << tmpline;
	    cout<<endl;
	    continue;
	}

	valuep = index(tmpline,' ');
	if (valuep == NULL) {
	    cerr << "Error while reading parameter name from "<<inp_file<<" at parameter #"<<i<<"\n";
	    exit(1);
	}
	*valuep = '\0';
	par = tmpline;
	
	if (sscanf(valuep+1,"%lf",&parv) != 1) {
	    cerr << "Error while reading value for "<<par<<" in "<<inp_file<<" (parameter #"<<i<<")\n";
	    exit(1);
	}
	inputs[par] = parv;
	cout <<"#"<< i <<" "<<par <<" = "<< inputs[par];
	cout<<endl;
    }
    cout<<endl;
    cout <<"Finished reading input file "<<inp_file<<".\n";
    fclose (inf);

    ////////////////////////////////////////////////////////////////////////
    ///// read parameters through piping a file that has them //////////////

    //check for settings in the file:
    CHECK_FOR_INT_SETTING("L",L);
    CHECK_FOR_INT_SETTING("Runs",Runs);
    CHECK_FOR_INT_SETTING("sim_length",sim_length);
    CHECK_FOR_PARAMETER("host_density",host_density);
    CHECK_FOR_INT_SETTING("place_virus_type",place_virus_type);
    CHECK_FOR_INT_SETTING("place_virus_num",place_virus_num);
    CHECK_FOR_INT_SETTING("start_plaques",start_plaques);
    CHECK_FOR_INT_SETTING("max_start_strains",max_start_strains);
    CHECK_FOR_INT_SETTING("start_neur_strain",start_neur_strain);
    CHECK_FOR_INT_SETTING("max_neur_drips",max_neur_drips);
    CHECK_FOR_INT_PARAMETER("max_rate",max_rate);

    CHECK_FOR_PARAMETER("uninf_cell_death",uninf_cell_death);
    CHECK_FOR_PARAMETER("infectivity_free_A",infectivity_free[0]);
    CHECK_FOR_PARAMETER("infectivity_free_B",infectivity_free[1]);
    CHECK_FOR_PARAMETER("infectivity_free_AB",infectivity_free[2]);
    CHECK_FOR_PARAMETER("cell_div",cell_div);
    CHECK_FOR_PARAMETER("inf_cell_death_A",inf_cell_death[0]);
    CHECK_FOR_PARAMETER("inf_cell_death_B",inf_cell_death[1]);
    CHECK_FOR_PARAMETER("inf_cell_death_AB",inf_cell_death[2]);
    CHECK_FOR_PARAMETER("viral_lag",viral_lag);
    CHECK_FOR_PARAMETER("vir_prod_rate_mean",vir_prod_rate_mean);
    CHECK_FOR_PARAMETER("diff_rate_virus",diff_rate_virus);
    CHECK_FOR_PARAMETER("freevirus_decay",freevirus_decay);
    CHECK_FOR_PARAMETER("neur_drip",neur_drip);
    CHECK_FOR_PARAMETER("neur_infect",neur_infect);
    CHECK_FOR_PARAMETER("prob_coinfect",prob_coinfect);
    CHECK_FOR_PARAMETER("prob_recomb",prob_recomb);

    CHECK_FOR_PARAMETER("psi",psi);

    CHECK_FOR_SETTING("site_scale",site_scale);
    CHECK_FOR_SETTING("refr_freq",refr_freq);
    CHECK_FOR_SETTING("data_freq",data_freq);
    CHECK_FOR_INT_SETTING("save_plot_data",save_plot_data);
    CHECK_FOR_INT_SETTING("save_state_files",save_state_files);
    CHECK_FOR_INT_SETTING("stop_early",stop_early);



}


// Usage printout (for -h option)
void usage(char *prog_name)
{

    fprintf(stderr,
	"Usage: %s [-h][-f <input_file>][-d][-o <dir>][-p <prog>][-r][-s <seed>][-v]\n",prog_name);
    fprintf(stderr,"\t-h = this help\n");
    fprintf(stderr,"\t-f = optional input file\n");
    fprintf(stderr,"\t-o = output directory for results.csv and cell state files\n");
    fprintf(stderr,"\t-p = Rscript program for cell state files (plotting)\n");
    fprintf(stderr,"\t-r = random number seed truely random\n");
    fprintf(stderr,"\t-s <seed> = random number seed for repeating run (unsigned)\n");
}

// 0.  MAIN FUNCTION-only used if there is no wrapper file
int main(int argc, char **argv){
    
    long seed;
    int verbose = 0;
    char *inp_file=NULL;

    global_settings settings;
    global_parameters parameters;
    global_dynamics dynamics;

    set_default_settings(&settings);

    set_default_parameters(&parameters);
    
    for(int i = 0; i < argc; i++) {
	if ((!strcmp(argv[i],"-f") || !strcmp(argv[i],"-F")) && i +1 < argc) {
	    inp_file = argv[i+1];
	    i++;
	}
	if ((!strcmp(argv[i],"-o") || !strcmp(argv[i],"-O")) && i +1 < argc) {
	    settings.dirname = argv[i+1];
	    DIR* directory = opendir(settings.dirname);
	    if(directory == NULL){
		fprintf(stdout,"Failed to set Output directory to %s\n",settings.dirname);
		exit(1);
	    } else {
		fprintf(stdout,"Output directory set to %s\n",settings.dirname);
		closedir(directory);
	    }
	    i++;
	}
	if ((!strcmp(argv[i],"-p") || !strcmp(argv[i],"-P")) && i +1 < argc) {
	    settings.process_cmd = argv[i+1];
	    i++;
	}
	if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"-R")) {
	    seed = time (NULL) * getpid();    
	    gsl_rng_set (settings.ur, seed);
	    fprintf(stdout,"Random generator seed is %ld\n",seed);
	}
	if (i +1 < argc && !strcmp(argv[i],"-s") && sscanf(argv[i+1],"%ld",&seed)==1) {
	    gsl_rng_set (settings.ur, seed);
	    fprintf(stdout,"Random generator seed is %ld\n",seed);
	    i++;
	}
	if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"-V")) {
	    verbose = 1;
	}
	if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"-H")) {
	    usage(argv[0]);
	    exit(1);
	}
    }
    if (inp_file != NULL)
	read_input_file(inp_file, &settings, &parameters);


    string data1;
    if (settings.dirname != NULL) {
	data1=settings.dirname;
	data1+="/";
    }
    else
	data1="";

    data1+=settings.out_fname;

    settings.verbose = verbose;

    if((settings.dataF1 = fopen(data1.c_str(),"wt")) == NULL){
	fprintf(stderr,"Could not open data file %s\n",data1.c_str());
	exit(1);
    }
    //modify settings and parameters here or use a wrapper, e.g. Spatial_wrapper.R
    //initialize and run sim
    dynamics.runnum=0;

    alloc_grid(&settings,&dynamics);

    output_results(dynamics.runnum,&settings,&parameters,&dynamics,true);
    while (dynamics.runnum < settings.Runs) {
	init_grid(dynamics.runnum+1,&settings,&parameters,&dynamics);

	simulation(dynamics.runnum+1,&settings,&parameters,&dynamics);
	dynamics.runnum++;
    }

    if(settings.dataF1 != NULL) fclose(settings.dataF1);
}

// FUNCTIONS: all the functions that are needed for the simulation are below this line

// 1a.  Get default settings for the simulation (can be overriden in input file)
void set_default_settings(global_settings *settings){

	const gsl_rng_type * T;//pointer to gsl rand generator type

	/* create a genrator chosen by the environment variable GSL_RNG_TYPE */
	gsl_rng_env_setup();

	T = gsl_rng_default;

	settings->Runs = 1;

	settings->ur = gsl_rng_alloc (T);

	settings->out_fname="results.csv";							//File to store cell and virus counts
	settings->dirname=NULL;

	settings->process_cmd=NULL;

	settings->L=200;											//Length of simulation in hours	
	settings->sim_length=5*24;											//Length of simulation in hours
	settings->refr_freq=0.1;  // frequency for screen capture
	settings->data_freq=0.1;  // frequency for data capture

	settings->save_plot_data=1;  // save cell/t-cell/viral data for R-style plots
	settings->save_ratio_data=0;  // save e/t ratio data for R-style plots
	settings->save_state_files=0;  // save cell&cyto state files during plotting
	settings->stop_early=1;  // stop when virus is gone
	
	settings->site_scale=50; //how many microns is one site
	
	//Place virus
	settings->place_virus_type=INFECTED_PRODUCER;
	settings->place_virus_num=1000;
	settings->start_plaques=NUM_STRAINS;
	settings->max_start_strains=NUM_STRAINS;
	settings->start_neur_strain=0;
	settings->max_neur_drips=1;
	
}

// 1b. Get default values for parameters  (can be overriden in input file)
void set_default_parameters(global_parameters *parameters){
	
	//Rates (per site or cell per day)
	parameters->max_rate=24;						  //Sets time resolution, should be greater than or equal to sum of all rates for a given state
	
	//Uninfected cells
	parameters->uninf_cell_death=0;	//Rate of uninfected cell death 
	for (int s=0; s < NUM_STRAINS; s++)
	    parameters->infectivity_free[s]=0.1;				//Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site

	parameters->cell_div=0.077; 					//Rate of generation of new susceptible cells (repop of empty sites)
	parameters->host_density=1;           //Fraction of sites occupied by susceptible cells at start	
	
	//Infected cells
	parameters->inf_cell_death[0]=1.25;			//Rate of infected cell death due to infection (strain A)
	parameters->inf_cell_death[1]=1.25;			//Rate of infected cell death due to infection (strain B)
	parameters->inf_cell_death[2]=0.25;			//Rate of infected cell death due to infection (strain AB)

	parameters->viral_lag=8;						  //Rate at which inf non-producer becomes viral producing cell (per cell per day)
	parameters->vir_prod_rate_mean=100000;	//Rate at which an infected cell produces free virus (per cell per day, mean and sd)
	
	//Virus
	parameters->diff_rate_virus=12; 						  //Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
	parameters->freevirus_decay=8.8;		//Rate of decay of free virus (per day)
	parameters->neur_drip=0;			//Rate of arrival of virus from neurons (per site per day)
	parameters->neur_infect=0;			//Rate of cellular infection from neurons (per site per day)
	parameters->psi=1;
	parameters->diff_range=3;						  //n x n site, default n = 3, i.e. immediate Moore neighbourhood. 


}



// 2a. Initialize grid and other arrays (memory allocation)
void alloc_grid(global_settings *settings,global_dynamics *dynamics){

	
	int *visited = (int  *)my_malloc(settings->L*settings->L * sizeof (int));
	if (visited == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating cell Matrix. Exiting.");
	    exit(1);
	}
	dynamics->visited = visited;

	//Create arrays
	int **cell_state = (int  **)my_malloc(settings->L * sizeof (int *));
	if (cell_state == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating cell Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    cell_state[i] = (int *)my_malloc(settings->L * sizeof (int));
	    if (cell_state[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating cell Matrix. Exiting.");
		exit(1);
	    }
	}
	dynamics->cell_state = cell_state;

	int ***viral_matrix = (int ***)my_malloc(NUM_STRAINS * sizeof (int**));
	if (viral_matrix == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating viral_matrix Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < NUM_STRAINS; i++)
	{
	    viral_matrix[i] = (int **)my_malloc(settings->L * sizeof (int *));
	    if (viral_matrix[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating viral_matrix Matrix. Exiting.");
		exit(1);
	    }
	    for (int j=0; j < settings->L; j++)
	    {
		viral_matrix[i][j] = (int *)my_malloc(settings->L * sizeof (int));
		if (viral_matrix[i][j] == NULL)
		{
		    fprintf(stderr,"Encountered a problem creating viral_matrix Matrix. Exiting.");
		    exit(1);
		}
	    }
	}
	dynamics->viral_matrix = viral_matrix;
	
	int ***delta_viral_matrix = (int ***)my_malloc(NUM_STRAINS * sizeof (int **));
	if (delta_viral_matrix == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating delta_viral_matrix Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < NUM_STRAINS; i++)
	{
	    delta_viral_matrix[i] = (int **)my_malloc(settings->L * sizeof (int *));
	    if (delta_viral_matrix[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating delta_viral_matrix Matrix. Exiting.");
		exit(1);
	    }
	    for (int j=0; j < settings->L; j++)
	    {
		delta_viral_matrix[i][j] = (int *)my_malloc(settings->L * sizeof (int));
		if (delta_viral_matrix[i][j] == NULL)
		{
		    fprintf(stderr,"Encountered a problem creating delta_viral_matrix Matrix. Exiting.");
		    exit(1);
		}
	    }
	}
	dynamics->delta_viral_matrix = delta_viral_matrix;
	
	
}

// 2b. Initialize grid and other arrays (initial values for fresh runs)
void init_grid(int runnum, global_settings *settings,global_parameters *parameters,global_dynamics *dynamics){
	//Compute simulation length in number of updates
	dynamics->tot_updates=settings->sim_length*parameters->max_rate;
	fprintf(stderr,"Starting run %d (%d updates expected).. hold on to your hat!\n",
		runnum, dynamics->tot_updates);

	//Dynamic variables - this stores time elapsed and population counts
	dynamics->num_updates=0;
	dynamics->time=0;
	dynamics->susceptible_cells=0;
	dynamics->infected_nonprodcells=0;
	dynamics->infn_A=0;
	dynamics->infn_B=0;
	dynamics->infn_ABn=0;
	dynamics->infn_ABr=0;
	dynamics->infected_prodcells=0;
	dynamics->infp_A=0;
	dynamics->infp_B=0;
	dynamics->infp_ABn=0;
	dynamics->infp_ABr=0;
	dynamics->dead_cells=0;
	dynamics->virions=0;
	dynamics->virionsA=0;
	dynamics->virionsB=0;
	dynamics->virionsAB=0;
	dynamics->max_vl_time=0;
	dynamics->max_vl=0;
	dynamics->viral_cells=0;
	dynamics->first_samp_vl_time=0;
	dynamics->first_samp_vl=0;
	dynamics->max_samp_vl_time=0;
	dynamics->max_samp_vl=0;
	
	clear_counts(settings,parameters,dynamics);
}

// reset dynamic state array counts for each new run...
void clear_counts(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics) {

	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++) {
		for (int s=0; s < NUM_STRAINS; s++)
		    dynamics->viral_matrix[s][i][j]=0;
		dynamics->cell_state[i][j]=EMPTY;
	    }

	dynamics->tot_infs=0;
	dynamics->tot_virons=0;



}

// clear deltas associated with a new time step (applied together at end of the step)
void clear_deltas(global_settings *settings,global_dynamics *dynamics) {
	for (int s=0; s < NUM_STRAINS; s++)
	    for (int i=0; i < settings->L; i++)
		for (int j=0; j < settings->L; j++)
		    dynamics->delta_viral_matrix[s][i][j]=0;
}

// apply deltas associated with a new time step 
void propagate_deltas(global_settings *settings,global_dynamics *dynamics) {
	for (int s=0; s < NUM_STRAINS; s++)
	    for (int i=0; i < settings->L; i++)
		for (int j=0; j < settings->L; j++)
		    dynamics->viral_matrix[s][i][j]=MAX(0,dynamics->viral_matrix[s][i][j] + dynamics->delta_viral_matrix[s][i][j]);


}



//PLACE HOST cells (epithelials)
void place_host(global_settings *settings,global_parameters *parameters, global_dynamics *dynamics){
  if(parameters->host_density==1){
	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<settings->L;j++)
		dynamics->cell_state[i][j]=SUSCEPTIBLE;
  }else{
    for (int i=0; i<settings->L;i++)
	for (int j=0; j<settings->L;j++)
	    if (gsl_rng_uniform(settings->ur)<parameters->host_density)
		dynamics->cell_state[i][j]=SUSCEPTIBLE;
  }
}

//PLACE VIRUS (initial innoculum - if any)
void place_virus(global_settings *settings, global_dynamics *dynamics){
    if(settings->place_virus_num==0 && settings->start_plaques == 0)
      return; //do nothing and exit

    if(settings->place_virus_type==SINGLE_POINT){ //place virus at the center of the grid
      dynamics->viral_matrix[0][settings->L/2][settings->L/2]=settings->place_virus_num; 
      
    }else if(settings->place_virus_type==RANDOM_VIRUS){ //virus spread randomly on grid

      int s=0;
      for (int plaqs=0; plaqs<settings->start_plaques;plaqs++) {
	  int i = gsl_rng_uniform_int(settings->ur,settings->L);
	  int j = gsl_rng_uniform_int(settings->ur,settings->L);
	  dynamics->viral_matrix[s][i][j]=settings->place_virus_num;
	  s++;
	  if (s >= settings->max_start_strains) s=0;
      }
      
    }else if(settings->place_virus_type==INFECTED_PRODUCER){ //place infected cell(s) instead of virus
      if(settings->start_plaques==1){ //at the center of the grid
	dynamics->cell_state[settings->L/2][settings->L/2]=INFP_A; 
	dynamics->tot_infs++;
      }else{ //randomly on the grid
	int s=0;
	for (int plaqs=0; plaqs<settings->start_plaques;plaqs++) {
	    int i = gsl_rng_uniform_int(settings->ur,settings->L);
	    int j = gsl_rng_uniform_int(settings->ur,settings->L);
	    if (s==0)
		dynamics->cell_state[i][j]=INFP_A;
	    else if (s==1)
		dynamics->cell_state[i][j]=INFP_B;
	    else
		dynamics->cell_state[i][j]=INFP_ABr;

	    dynamics->tot_infs++;
	    s++;
	    if (s >= settings->max_start_strains) s=0;
	}
      }
    }
}


//NEURONAL DRIP (future feature - untested in this single episode version)
void neur_drip(int s, int i,int j, global_dynamics *dynamics){
	dynamics->delta_viral_matrix[s][i][j]+=1;
}

//NEURONAL INFECT (future feature - untested in this single episode version)
void neur_infect(int s, int i,int j, global_dynamics *dynamics){
	fprintf(stdout,"Infected cell at site: %d,%d with type %d at time = %g\n,",i,j,s,dynamics->time);
	if (s==0)
	    dynamics->cell_state[i][j] = INFP_A;
	else if (s==1)
	    dynamics->cell_state[i][j] = INFP_B;
	else 
	    dynamics->cell_state[i][j] = INFP_ABr;
}

//Calculates plaque diameter in mm
double calc_plaque_size(global_settings *settings, global_dynamics *dynamics){
   return(sqrt(dynamics->dead_cells)*settings->site_scale/1000);
}


//COMPUTE GAUSSIAN kernel - only run once, at the beginning of simulation
gsl_matrix *compute_gaussian(double mask_range,double sig,double amp){
	gsl_matrix *gaussian_mask;
	gaussian_mask=gsl_matrix_alloc(mask_range,mask_range);
	double temp=(mask_range-1)/2;
	double x=-temp;
	double y=-temp;
	double gaus_sum = 0;

        for (int row=0; row < mask_range; row++) {
	    for (int col=0; col < mask_range; col++) {
		double value =(amp/(2*M_PI*(sig*sig)))*exp(-(x*x+y*y)/(2*sig*sig)); 
		gsl_matrix_set(gaussian_mask,row,col,value);
		gaus_sum += value;
		x=x+1;
	    }
	    x=-temp;
	    y=y+1;
	}

	//this makes everything sum to 1
        for (int row=0; row < mask_range; row++)
	    for (int col=0; col < mask_range; col++)
	{
	    double value = gsl_matrix_get(gaussian_mask,row,col);
	    gsl_matrix_set(gaussian_mask,row,col,value/gaus_sum);
	}

	return(gaussian_mask);
}

//DIFFUSION of free virus
void diffusion_virus(global_settings *settings,int s,int i, int j,gsl_matrix *gaussian_mask,global_parameters *parameters, global_dynamics *dynamics){
  bool got_neighbor;
  int nbr[2];

  dynamics->delta_viral_matrix[s][i][j]-=dynamics->viral_matrix[s][i][j]; 
  if(dynamics->viral_matrix[s][i][j]==1){
    got_neighbor=pick_nbr(settings,i,j,nbr);
    if(got_neighbor){
	    dynamics->delta_viral_matrix[s][nbr[0]][nbr[1]]+=1;
    }
  }else{
	double burst = dynamics->viral_matrix[s][i][j];
	//virus at focal site [i][j] will be spread across diffusion range
	// note: some will be placed at [i][j] as well
	distrib_progeny(settings,parameters,dynamics,VIRUS,s, i,j,burst,gaussian_mask); 
  }
}

//DISTRIBUTE PROGENY in free virus or cytokine  array for burst or diffusion
//set into delat arrays for distributing at the end of the processing loop
//(avoids double diffusing) 
void distrib_progeny(global_settings *settings,global_parameters *parameters,
	global_dynamics *dynamics, int what, int s, int i, int j,double burst,
	gsl_matrix *gaussian_mask){
    int dist_range = parameters->diff_range;
    int temp2=(dist_range-1)/2;
    int imin=MAX(0,i-temp2);
    int imax=MIN(i+temp2,settings->L-1); //deal with edges
    int jmin=MAX(0,j-temp2);
    int jmax=MIN(j+temp2,settings->L-1);
    int distributed=0;
    for (int row=imin; row <= imax; row++)
	for (int col=jmin; col <= jmax; col++)
	{
	    int i1=MAX(0,row-imin);
	    int j1=MAX(0,col-jmin);
	    if (what == VIRUS) {
		int this_batch = 
		    round(burst * gsl_matrix_get(gaussian_mask,i1,j1));
		distributed+=this_batch;
		dynamics->delta_viral_matrix[s][row][col]+=this_batch;
	    }
	}
    if (distributed != burst && what == VIRUS) {
	int extra=distributed-(int)burst;
	while(extra != 0)
	{
	    int row=imin+gsl_rng_uniform_int(settings->ur,imax-imin+1);
	    int col=jmin+gsl_rng_uniform_int(settings->ur,jmax-jmin+1);

	    if (distributed>burst) {
		dynamics->delta_viral_matrix[s][row][col]--;
		extra--;
	    } else {
		dynamics->delta_viral_matrix[s][row][col]++;
		extra++;
	    }
	}
    }
}


//PICK NEIGHBOUR from immediate moore neighbourhood
bool pick_nbr(global_settings *settings,int i,int j, int *coords){
	int delta_i=gsl_rng_uniform_int(settings->ur,3);
	int delta_j=gsl_rng_uniform_int(settings->ur,3);

	// repick if no move selected
	while (delta_i == 1 && delta_j == 1) {
	    delta_i=gsl_rng_uniform_int(settings->ur,3);
	    delta_j=gsl_rng_uniform_int(settings->ur,3);
	}

	coords[0]=i+delta_i-1;
	coords[1]=j+delta_j-1;

	if (coords[0] < 0 || coords[0] >= settings->L ||
	    coords[1] < 0 || coords[1] >= settings->L) {
		return(false); //outside grid
	}else{
		return(true);
	}
}

//DETECT stuff within detection range
//If return_pos == COORDS, function will return the position of the detected object (i,j coords)
//If there are multiple objects detected, one will be picked at random and its position 
//returned
//If return_pos == SUM, function will return the total number of "what"s in the nbhd
//If return_pos == CHECK, function will return one if any number of "what"s are in the nbhd
int detect_inf(global_settings *settings, global_dynamics *dynamics, 
	int det_range,int i, int j,int what,int return_pos,int *coords){

  int det_inf = 0;

  if(what==NONE){
    return(1);
  }else{
    int imin=MAX(0,i-(det_range-1)/2);
    int imax=MIN(i+(det_range-1)/2,settings->L-1);
    int jmin=MAX(0,j-(det_range-1)/2);
    int jmax=MIN(j+(det_range-1)/2,settings->L-1);
    
    if(return_pos==CHECK){		
      if(what==ANY){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD)
		    return 1;
	    }
      }else if(what==INF){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr ||
		    dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
		    return 1;
	    }
      }else if(what==INFP){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr)
		    return 1;
	    }
      }else if(what==INFN){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
		    return 1;
	    }
      }else if(what==ALL_CELLS){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD)
		    return 1;
	    }
      }else if(what==UNINF){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == SUSCEPTIBLE ||
		    dynamics->cell_state[row][col] == UNPROTECTABLE)
		    return 1;
	    }
      }
      return(det_inf);
      
    }else if(return_pos==COORDS || return_pos==SUM){
      if(what==ANY){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==INF){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr ||
		    dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr ||
		    dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==INFP){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->cell_state[row][col] == INFP_A ||
		    dynamics->cell_state[row][col] == INFP_B ||
		    dynamics->cell_state[row][col] == INFP_ABn ||
		    dynamics->cell_state[row][col] == INFP_ABr)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==INFN){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->cell_state[row][col] == INFN_A ||
		    dynamics->cell_state[row][col] == INFN_B ||
		    dynamics->cell_state[row][col] == INFN_ABn ||
		    dynamics->cell_state[row][col] == INFN_ABr)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==ALL_CELLS){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->cell_state[row][col] > EMPTY &&
		    dynamics->cell_state[row][col] < DEAD)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==UNINF){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == SUSCEPTIBLE ||
		    dynamics->cell_state[row][col] == UNPROTECTABLE)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->cell_state[row][col] == SUSCEPTIBLE ||
		    dynamics->cell_state[row][col] == UNPROTECTABLE)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else{
	    fprintf(stderr,"Not a valid option for what in detect_inf().\n");
	    exit(1);
      }	
      return(det_inf);
    }else{
	fprintf(stderr,"Not a valid option for return_pos in detect_inf().\n");
	exit(1);
    }	
  }
}


//DECAY of free virus
void viral_decay(int s,int i, int j,global_parameters *parameters, global_dynamics *dynamics){
    dynamics->viral_matrix[s][i][j]=MAX(0,dynamics->viral_matrix[s][i][j]*exp(-parameters->freevirus_decay/parameters->max_rate));
}


void set_infection_type(int s, int i, int j, global_settings *settings, global_parameters *parameters, global_dynamics *dynamics, bool *new_infection){

      // Let A infect susc & change B into ABn or ABr
      if (s==0 && dynamics->cell_state[i][j]!=INFN_A && dynamics->cell_state[i][j]!=INFN_ABn && dynamics->cell_state[i][j]!=INFN_ABr)
      {
	  if (dynamics->cell_state[i][j]!=INFN_B)
	  {
	      dynamics->cell_state[i][j]=INFN_A; //becomes an infected cell (strain A)
	      *new_infection=true;
	  }
	  else if (gsl_rng_uniform(settings->ur) < parameters->prob_coinfect)
	  {
	      if (gsl_rng_uniform(settings->ur) < parameters->prob_recomb)
		  dynamics->cell_state[i][j]=INFN_ABr; //becomes an infected cell (strain ABr)
	      else
		  dynamics->cell_state[i][j]=INFN_ABn; //becomes an infected cell (strain ABn)
	      *new_infection=true;
	  }
      }
      // Let B infect susc & change A into ABn or ABr
      else if (s==1 && dynamics->cell_state[i][j]!=INFN_B && dynamics->cell_state[i][j]!=INFN_ABn && dynamics->cell_state[i][j]!=INFN_ABr)
      {
	  if (dynamics->cell_state[i][j]!=INFN_A)
	  {
	      dynamics->cell_state[i][j]=INFN_B; //becomes an infected cell (strain B)
	      *new_infection=true;
	  }
	  else if (gsl_rng_uniform(settings->ur) < parameters->prob_coinfect)
	  {
	      if (gsl_rng_uniform(settings->ur) < parameters->prob_recomb)
		  dynamics->cell_state[i][j]=INFN_ABr; //becomes an infected cell (strain ABr)
	      else
		  dynamics->cell_state[i][j]=INFN_ABn; //becomes an infected cell (strain ABn)
	      *new_infection=true;
	  }
      }
      // Let AB recomb virus coinfect other types & if so switch to recomb (or just susc as ABr)
      else if (s==2) // AB recomb virus (check for coinfection w/ other types & if so switch to recomb
      {
	  if (dynamics->cell_state[i][j]==INFN_A || dynamics->cell_state[i][j]==INFN_B || dynamics->cell_state[i][j]==INFN_ABn)
	  {
		if (gsl_rng_uniform(settings->ur) < parameters->prob_coinfect)
		{
		  dynamics->cell_state[i][j]=INFN_ABr; //becomes an infected cell (strain ABr)
		  *new_infection=true;
		}
	  }
	  else
	  {
		dynamics->cell_state[i][j]=INFN_ABr; //becomes an infected cell (strain ABr)
		*new_infection=true;
	  }
      }
}
// 2. SIMULATION
void simulation(int runnum, global_settings *settings,global_parameters *parameters,global_dynamics *dynamics){
    
    int coords[2];
    double next_data_capture = 0;
    double next_screen_capture = 0;

    int infections_avoided=0;
    unsigned long int virions_avoided=0;

    int vdist_range = MAX(parameters->diff_range,1+2*parameters->diff_rate_virus/parameters->max_rate);
    gsl_matrix *gaussian_mask_vir_diff=compute_gaussian(vdist_range,1.2,1); 

     
    place_host(settings,parameters,dynamics); //Place host cells

    if(settings->place_virus_type!=NO_VIRUS) 
	    place_virus(settings,dynamics); //Place virus	

    
    dynamics->num_updates=0;
    compute_totals(settings,parameters,dynamics);
    output_results(runnum,settings,parameters,dynamics,false);
    if (settings->save_plot_data || settings->save_state_files)
	dump_data(settings,parameters,dynamics);
    

    next_screen_capture += settings->refr_freq;
    next_data_capture += settings->data_freq;

    int last_neuro_type=settings->start_neur_strain;  // allows any strain to interact with starting plaques
    int num_drips=0;

    while(dynamics->num_updates<dynamics->tot_updates){
	
	clear_deltas(settings,dynamics);

	////Any site on grid
	////------------------
	////Viral drip from neuron (DOES NOT PRECLUDE OTHER ACTIONS)
	if (parameters->neur_drip!=0 ) {
	    double dice_drip=gsl_rng_uniform(settings->ur); // neuro drip
	    int drip_x=gsl_rng_uniform_int(settings->ur,settings->L);
	    int drip_y=gsl_rng_uniform_int(settings->ur,settings->L);
	    if(num_drips < settings->max_neur_drips && dice_drip<=parameters->neur_drip/parameters->max_rate){
		  neur_drip(last_neuro_type,drip_x,drip_y,dynamics);
		  last_neuro_type++;
		  if (last_neuro_type >= settings->max_start_strains) last_neuro_type=0;
		  num_drips++;
	    }
	}
      
	////Cellular infection from neuron (DOES NOT PRECLUDE OTHER ACTIONS)
	if (parameters->neur_infect!=0 ) {
	    double dice_infect=gsl_rng_uniform(settings->ur); // neuro infect
	    int infect_x=gsl_rng_uniform_int(settings->ur,settings->L);
	    int infect_y=gsl_rng_uniform_int(settings->ur,settings->L);
	    if(num_drips < settings->max_neur_drips && dice_infect<=parameters->neur_infect/parameters->max_rate){
		  neur_infect(last_neuro_type,infect_x,infect_y,dynamics);
		  last_neuro_type++;
		  if (last_neuro_type >= settings->max_start_strains) last_neuro_type=0;
		  num_drips++;
	    }
	}
      
	// Like the R version, the C++ version evaluates cells
	// in a random order.  In addition, it calculates all deltas first 
	// before applying them across the grid (to avoid double moves, etc)
	for(int sites=0; sites <settings->L*settings->L;sites++)
	    dynamics->visited[sites]=0;

	for(int sites=0; sites <settings->L*settings->L;sites++) {
	     int this_site;
	     do {
		this_site=gsl_rng_uniform_int(settings->ur,settings->L*settings->L);
	     } while(dynamics->visited[this_site]==1);
	     dynamics->visited[this_site]=1;
	     int i = this_site/settings->L;
	     int j = this_site - (i*settings->L);

		double dice=gsl_rng_uniform(settings->ur); // generic dice
		double dice_v=gsl_rng_uniform(settings->ur); // virus related
		double dice_t=gsl_rng_uniform(settings->ur); // t-cell related
		double dice_cyt=gsl_rng_uniform(settings->ur); // cytokine related
		double dice_cell_death=gsl_rng_uniform(settings->ur); // cell death
		double dice_cell_lag=gsl_rng_uniform(settings->ur); // lag on p

		double cytokine_factor=1;
		    	
					      
		////Host cells
		//// 0: vacant site
		//// 1: susceptible keratinocyte
		//// 2: infected, not producing virus
		//// 3: infected, producing virus
		//// 4: dead keratinocyte
		////--------------------------------------------
		//State 0 or 4: vacant site or occupied by dead cell
		if(dynamics->cell_state[i][j]==EMPTY||
		      dynamics->cell_state[i][j]>=DEAD){
	    
		    //Event 0.1: cell renewal, occupation by susceptible
		    if(dice<=(parameters->cell_div)/parameters->max_rate){
			    dynamics->cell_state[i][j]=SUSCEPTIBLE;
		    }
    
	    //State 1: susceptible or infected cell (cytokines can reduce
	    //infectivity while in this state buit only for a period set using
	    // parameters->cyt_prot_expiration (1 day?)
	    // infected cells can be coinfected

	    //State 2: infection/production
		}else if(dynamics->cell_state[i][j]==SUSCEPTIBLE ||
			dynamics->cell_state[i][j]==INFN_A ||
			dynamics->cell_state[i][j]==INFN_B ||
			dynamics->cell_state[i][j]==INFN_ABn ||
			dynamics->cell_state[i][j]==INFN_ABr ||
			dynamics->cell_state[i][j]==UNPROTECTABLE){

		
		    //Event 1.1: infection due to free virus at the site
		    bool virus_present=false;
		    bool new_infection=false;
		    for (int s=0; s < NUM_STRAINS; s++)
		    {
		      if(dynamics->viral_matrix[s][i][j]>0){
			virus_present=true;

			    if(dice<=(
				parameters->infectivity_free[s]*
				dynamics->viral_matrix[s][i][j])/
				parameters->max_rate){
			      set_infection_type(s, i, j, settings, parameters, dynamics, &new_infection);
			}
			if (new_infection)
			{
			  dynamics->tot_infs++;
			  dynamics->delta_viral_matrix[s][i][j]--;
			}
			new_infection=false;
		      }
		    }

		    //Event 1.2: apoptosis
		      //without cytokine effects
		    if (!virus_present)
		    {
			if(
			    
			     dice<=(parameters->uninf_cell_death)/parameters->max_rate){
			  dynamics->cell_state[i][j]=DEAD;
			}

		    } 
	
		//State 2b: infected cell (non-virus producing, not presenting Ag)
		    if (dynamics->cell_state[i][j] == INFN_A ||
			dynamics->cell_state[i][j] == INFN_B ||
			dynamics->cell_state[i][j] == INFN_ABn ||
			dynamics->cell_state[i][j] == INFN_ABr) {

			double inf_cell_death=0;
			if (dynamics->cell_state[i][j] == INFN_ABr)
			    inf_cell_death = parameters->inf_cell_death[2];
			else if (dynamics->cell_state[i][j] == INFN_B)
			    inf_cell_death = parameters->inf_cell_death[1];
			else 
			    inf_cell_death = parameters->inf_cell_death[0];
			//Event 2.0: produce cytokine

		      if(dice_cell_lag<=(parameters->viral_lag)/parameters->max_rate){
			if (dynamics->cell_state[i][j] == INFN_A) dynamics->cell_state[i][j]=INFP_A; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_B) dynamics->cell_state[i][j]=INFP_B; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABn) dynamics->cell_state[i][j]=INFP_ABn; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABr) dynamics->cell_state[i][j]=INFP_ABr; //becomes virus-producing infected cell
		      
			  //with cytokine effects
				    
			    //Event 2.2: dies 
		      
	    
		      }else {
			if(dice_cell_death<=(inf_cell_death)/parameters->max_rate){
			    dynamics->cell_state[i][j]=DEAD; //death of infected cell
			    
			    
			  }
		      }
		  }
		  
		//State 3: infected cell (producing virus, presenting Ag)
		} else if (dynamics->cell_state[i][j] == INFP_A ||
		    dynamics->cell_state[i][j] == INFP_B ||
		    dynamics->cell_state[i][j] == INFP_ABn ||
		    dynamics->cell_state[i][j] == INFP_ABr) {
	
			double inf_cell_death=0;
			if (dynamics->cell_state[i][j] == INFP_ABr)
			    inf_cell_death = parameters->inf_cell_death[2];
			else if (dynamics->cell_state[i][j] == INFP_B)
			    inf_cell_death = parameters->inf_cell_death[1];
			else 
			    inf_cell_death = parameters->inf_cell_death[0];

			int s = gsl_rng_uniform_int(settings->ur,2);

			unsigned long int new_virions =
			    MAX(0, gsl_ran_weibull(settings->ur,parameters->vir_prod_rate_mean/parameters->max_rate,5.0));
		      //without cytokine effects
			    if (dynamics->cell_state[i][j] == INFP_A || dynamics->cell_state[i][j] == INFP_ABn)
				dynamics->delta_viral_matrix[0][i][j]+=new_virions;
			    if (dynamics->cell_state[i][j] == INFP_B || dynamics->cell_state[i][j] == INFP_ABn)
				dynamics->delta_viral_matrix[1][i][j]+=new_virions;
			    if (dynamics->cell_state[i][j] == INFP_ABr)
				dynamics->delta_viral_matrix[2][i][j]+=new_virions;
			    
			    dynamics->tot_virons+= new_virions;
			    if (settings->verbose)
			      fprintf(stderr,"Infected cell %d,%d made virus at t=%3.2lf\n",i,j,dynamics->time);
		  
			
			//Event 3.1: die
			if (dice_cell_death<=(inf_cell_death)/parameters->max_rate){
			  dynamics->cell_state[i][j]=DEAD; //death of infected cell
			  if (settings->verbose)
			    fprintf(stderr,"Infected cell %d,%d died at t=%3.2lf\n",i,j,dynamics->time);
			  
			}
		}
	

		////Cell-free virus
		////--------------------------------------------
		for (int s=0; s < NUM_STRAINS; s++)
		{
		    if(dynamics->viral_matrix[s][i][j]>0){

		      //Event 0a: Decay
			    viral_decay(s, i,j,parameters,dynamics);
		    
			    //Event 0b: Diffusion of remaining free virus at site
			    if(dice_v<=parameters->diff_rate_virus/parameters->max_rate){
				    diffusion_virus(settings,s,i,j,gaussian_mask_vir_diff,parameters,dynamics);}			
		    
		    }			
		}
	}
	propagate_deltas(settings,dynamics);
	
	dynamics->num_updates=dynamics->num_updates+1;
	compute_totals(settings,parameters,dynamics);
	if (dynamics->time >= next_data_capture) {
	    output_results(runnum,settings,parameters,dynamics,false);
	    next_data_capture += settings->data_freq;
	}
	if (dynamics->time >= next_screen_capture) {
	    if (settings->save_plot_data || settings->save_state_files)
		dump_data(settings,parameters,dynamics);
	    next_screen_capture += settings->refr_freq;
	}
	if(settings->stop_early && dynamics->virions==0 && //dynamics->virions<100 && 
	   parameters->neur_drip==0 &&
	   parameters->neur_infect==0 &&
	   dynamics->infected_prodcells==0){
	  break;
	}
    }
    

    
    end_simulation(dynamics);
}

// Compute totals of cells and virus
void compute_totals(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics){
	dynamics->time=(double)dynamics->num_updates/parameters->max_rate; // time in days!
	dynamics->virions=0;
	dynamics->virionsA=0;
	dynamics->virionsB=0;
	dynamics->virionsAB=0;
	dynamics->viral_cells=0;
	dynamics->susceptible_cells=0;
	dynamics->infected_nonprodcells=0;
	dynamics->infn_A=0;
	dynamics->infn_B=0;
	dynamics->infn_ABn=0;
	dynamics->infn_ABr=0;
	dynamics->infected_prodcells=0;
	dynamics->infp_A=0;
	dynamics->infp_B=0;
	dynamics->infp_ABn=0;
	dynamics->infp_ABr=0;
	dynamics->dead_cells=0;

	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<settings->L;j++)
	{
	    if (dynamics->cell_state[i][j]==SUSCEPTIBLE ||
		dynamics->cell_state[i][j]==UNPROTECTABLE)
		dynamics->susceptible_cells++;
	    if (dynamics->cell_state[i][j] == INFN_A ||
		dynamics->cell_state[i][j] == INFN_B ||
		dynamics->cell_state[i][j] == INFN_ABn ||
		dynamics->cell_state[i][j] == INFN_ABr)
		dynamics->infected_nonprodcells++;
	    if (dynamics->cell_state[i][j] == INFN_A)
		dynamics->infn_A++;
	    if (dynamics->cell_state[i][j] == INFN_B)
		dynamics->infn_B++;
	    if (dynamics->cell_state[i][j] == INFN_ABn)
		dynamics->infn_ABn++;
	    if (dynamics->cell_state[i][j] == INFN_ABr)
		dynamics->infn_ABr++;
	    if (dynamics->cell_state[i][j] == INFP_A ||
		dynamics->cell_state[i][j] == INFP_B ||
		dynamics->cell_state[i][j] == INFP_ABn ||
		dynamics->cell_state[i][j] == INFP_ABr)
		dynamics->infected_prodcells++;
	    if (dynamics->cell_state[i][j] == INFP_A)
		dynamics->infp_A++;
	    if (dynamics->cell_state[i][j] == INFP_B)
		dynamics->infp_B++;
	    if (dynamics->cell_state[i][j] == INFP_ABn)
		dynamics->infp_ABn++;
	    if (dynamics->cell_state[i][j] == INFP_ABr)
		dynamics->infp_ABr++;
	    if (dynamics->cell_state[i][j]>=DEAD)
		dynamics->dead_cells++;

	    bool has_virus = false;
	    for (int s=0; s < NUM_STRAINS; s++)
	    {
		dynamics->virions+=dynamics->viral_matrix[s][i][j];
		if (dynamics->viral_matrix[s][i][j] > 0)
		    has_virus=true;
		if (s==0)
		    dynamics->virionsA+=dynamics->viral_matrix[s][i][j];
		else if (s==1)
		    dynamics->virionsB+=dynamics->viral_matrix[s][i][j];
		else 
		    dynamics->virionsAB+=dynamics->viral_matrix[s][i][j];
	    }
	    if (has_virus)
		dynamics->viral_cells++;
	}
	if (dynamics->virions*parameters->psi > dynamics->max_vl){
	    dynamics->max_vl=dynamics->virions*parameters->psi;
	    dynamics->max_vl_time=dynamics->time;
	}

}

// This routine generates state data so that it can be used to make
// plots either while the sim runs (using -p <script> option) or
// at a later time on this or other platforms (see scripts directory)
void dump_data(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics){
    string output_file1;
    string output_file2;

    string output_dir;
    FILE * outF= NULL;

    char temp_str[100];

    sprintf(temp_str,"cell_state_%d_%4.3lf.csv",
	dynamics->runnum+1,dynamics->time);

    if (settings->dirname != NULL) {
	output_file1=settings->dirname;
	output_file1+="/";
	output_file1+=temp_str;
    }
    else
	output_file1=temp_str;

    if((outF = fopen(output_file1.c_str(),"wt")) == NULL){
	fprintf(stderr,"Could not open data file %s\n",output_file1.c_str());
	exit(1);
    }
    for (int i=0; i<settings->L-1;i++) {
	    fprintf(outF,"%d,",i+1);
    }
    fprintf(outF,"%d\n",settings->L);
    for (int i=0; i<settings->L;i++) {
	for (int j=0; j<settings->L-1;j++) {
	    fprintf(outF,"%d,",dynamics->cell_state[i][j]);
	}
	fprintf(outF,"%d\n",dynamics->cell_state[i][settings->L-1]);
    }
    fclose(outF);

    for (int s=0; s < NUM_STRAINS; s++) {
	sprintf(temp_str,"virus_state%s_%d_%4.3lf.csv",
	    (s==0)?"A":((s==1)?"B":"AB"),dynamics->runnum+1,dynamics->time);

	if (settings->dirname != NULL) {
	    output_file2=settings->dirname;
	    output_file2+="/";
	    output_file2+=temp_str;
	}
	else
	    output_file2=temp_str;

	if((outF = fopen(output_file2.c_str(),"wt")) == NULL){
	    fprintf(stderr,"Could not open data file %s\n",output_file2.c_str());
	    exit(1);
	}
	for (int i=0; i<settings->L-1;i++) {
		fprintf(outF,"%d,",i+1);
	}
	fprintf(outF,"%d\n",settings->L);
	for (int i=0; i<settings->L;i++) {
	    for (int j=0; j<settings->L-1;j++) {
		fprintf(outF,"%d,",dynamics->viral_matrix[s][i][j]);
	    }
	    fprintf(outF,"%d\n",dynamics->viral_matrix[s][i][settings->L-1]);
	}
	fclose(outF);
    }



    if (settings->process_cmd != NULL) {
	string command;
	command = "Rscript ";
	command += settings->process_cmd;
	command += " ";
	if (settings->dirname != NULL) {
	    command+=settings->dirname;
	}
	else
	    command+=".";

	command+= "/ ";

	sprintf(temp_str, "%d %4.3lf %4.3lf %4.3lf",
	    dynamics->runnum+1,dynamics->time, (dynamics->max_vl>0)?log10(dynamics->max_vl):0, dynamics->max_vl_time);

	command+= temp_str;

	//if (settings->verbose)
	if (settings->save_state_files)
	    fprintf(stdout,"%s\n",command.c_str());

	if (settings->save_plot_data) {
	    fprintf(stderr,"Launching command: %s\n",command.c_str());

	    int ret = system(command.c_str());
	    if (WIFSIGNALED(ret) &&
		(WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT)) {
		fprintf(stderr,"Rscript aborted.  Exiting!\n");
		exit(1);
	    }
	    if (!settings->save_state_files) {
		remove(output_file1.c_str());
		remove(output_file2.c_str());
	    }
	}
    }
}

// This routine writes to the results.csv file (time history data)
void output_results(int runnum, global_settings *settings,
	global_parameters *parameters,global_dynamics *dynamics,bool header){

    if (settings->dataF1==NULL) {
	fprintf(stderr,"No output file\n");
	exit(1);
    }
    // update sampled peak and its timestamp
    if (dynamics->virions*parameters->psi > dynamics->max_samp_vl) {
	dynamics->max_samp_vl=dynamics->virions*parameters->psi;
	dynamics->max_samp_vl_time=dynamics->time;
    }
    if (dynamics->virions*parameters->psi > DETECTION_THRESHOLD &&
	dynamics->first_samp_vl == 0) {
	dynamics->first_samp_vl=dynamics->virions*parameters->psi;
	dynamics->first_samp_vl_time=dynamics->time;
    }
    if (header) {
	    fprintf(settings->dataF1,"run");
	    fprintf(settings->dataF1,",time");
	    fprintf(settings->dataF1,",susceptible");
	    fprintf(settings->dataF1,",inf_nonprod");
	    fprintf(settings->dataF1,",infn_A");
	    fprintf(settings->dataF1,",infn_B");
	    fprintf(settings->dataF1,",infn_ABn");
	    fprintf(settings->dataF1,",infn_ABr");
	    fprintf(settings->dataF1,",inf_prod");
	    fprintf(settings->dataF1,",infp_A");
	    fprintf(settings->dataF1,",infp_B");
	    fprintf(settings->dataF1,",infp_ABn");
	    fprintf(settings->dataF1,",infp_ABr");
	    fprintf(settings->dataF1,",dead_cells");
	    fprintf(settings->dataF1,",virions");
	    fprintf(settings->dataF1,",A virions");
	    fprintf(settings->dataF1,",B virions");
	    fprintf(settings->dataF1,",AB virions");
	    fprintf(settings->dataF1,",log virions");
	    fprintf(settings->dataF1,",max log VL");
	    fprintf(settings->dataF1,",max time");
	    fprintf(settings->dataF1,",viral_cells");
	    fprintf(settings->dataF1,",first sampled log VL");
	    fprintf(settings->dataF1,",first sampled time");
	    fprintf(settings->dataF1,",max sampled log VL");
	    fprintf(settings->dataF1,",max sampled time");
	    fprintf(settings->dataF1,",plaque size");
	    fprintf(settings->dataF1,",total infected");
	    fprintf(settings->dataF1,",log total virus");
	    fprintf(settings->dataF1,"\n");
	    fflush(settings->dataF1);
    } else {
	    fprintf(settings->dataF1,"%d", runnum);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->time);
	    fprintf(settings->dataF1,",%d", dynamics->susceptible_cells);
	    fprintf(settings->dataF1,",%d", dynamics->infected_nonprodcells);
	    fprintf(settings->dataF1,",%d", dynamics->infn_A);
	    fprintf(settings->dataF1,",%d", dynamics->infn_B);
	    fprintf(settings->dataF1,",%d", dynamics->infn_ABn);
	    fprintf(settings->dataF1,",%d", dynamics->infn_ABr);
	    fprintf(settings->dataF1,",%d", dynamics->infected_prodcells);
	    fprintf(settings->dataF1,",%d", dynamics->infp_A);
	    fprintf(settings->dataF1,",%d", dynamics->infp_B);
	    fprintf(settings->dataF1,",%d", dynamics->infp_ABn);
	    fprintf(settings->dataF1,",%d", dynamics->infp_ABr);
	    fprintf(settings->dataF1,",%d", dynamics->dead_cells);
	    fprintf(settings->dataF1,",%d", (int)(dynamics->virions*parameters->psi));
	    fprintf(settings->dataF1,",%d", (int)(dynamics->virionsA*parameters->psi));
	    fprintf(settings->dataF1,",%d", (int)(dynamics->virionsB*parameters->psi));
	    fprintf(settings->dataF1,",%d", (int)(dynamics->virionsAB*parameters->psi));
	    fprintf(settings->dataF1,",%4.3lf", (dynamics->virions*parameters->psi>0)?log10(dynamics->virions*parameters->psi):0);
	    fprintf(settings->dataF1,",%4.3lf", (dynamics->max_vl>0)?log10(dynamics->max_vl):0);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->max_vl_time);
	    fprintf(settings->dataF1,",%d", dynamics->viral_cells);
	    fprintf(settings->dataF1,",%4.3lf", (dynamics->first_samp_vl>0)?log10(dynamics->first_samp_vl):0);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->first_samp_vl_time);
	    fprintf(settings->dataF1,",%4.3lf", (dynamics->max_samp_vl>0)?log10(dynamics->max_samp_vl):0);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->max_samp_vl_time);
	    fprintf(settings->dataF1,",%4.3lf", calc_plaque_size(settings, dynamics));
	    fprintf(settings->dataF1,",%d", dynamics->tot_infs);
	    fprintf(settings->dataF1,",%4.3lf", (dynamics->tot_virons)?log10(dynamics->tot_virons):0);
	    fprintf(settings->dataF1,"\n");
	    fflush(settings->dataF1);
    }
}

// OUTPUT 
//END of simulation 
void end_simulation(global_dynamics *dynamics){
	fprintf(stdout,"Highest log VL %3.2lf at t=%3.2lf\n",(dynamics->max_vl>0)?log10(dynamics->max_vl):0,dynamics->max_vl_time);
	fprintf(stdout,"End of simulation at t=%3.2lf\n",dynamics->time);
}
