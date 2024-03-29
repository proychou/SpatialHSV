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

#ifdef TCELL_CODE
// Constants for Trm setup and behavior
#define RANDOM 0
#define GAUSSIAN 1
#define DIRECTED 2
#define LEVY 3
#define DISTRIB 4
#define DRAWN 5
#endif

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

#ifdef TCELL_CODE
    CHECK_FOR_INT_SETTING("plot_trms",plot_trms);
    CHECK_FOR_INT_SETTING("save_ratio_data",save_ratio_data);
    CHECK_FOR_INT_SETTING("trm_init",trm_init);
    CHECK_FOR_PARAMETER("trm_decay",trm_decay);
    CHECK_FOR_PARAMETER("trm_conv_rate",trm_conv_rate);
    CHECK_FOR_PARAMETER("trm_dbl",trm_dbl);
    CHECK_FOR_PARAMETER("tem_trafficking",tem_trafficking);
    CHECK_FOR_PARAMETER("tem_exiting",tem_exiting);
    CHECK_FOR_PARAMETER("tem_delay",tem_delay);

    CHECK_FOR_INT_SETTING("trm_move_method",trm_move_method);
    CHECK_FOR_INT_SETTING("trm_random_gaussian",trm_random_gaussian);

    CHECK_FOR_INT_SETTING("add_cd4s",add_cd4s);
    CHECK_FOR_SETTING("cd4_to_cd8_ratio",cd4_to_cd8_ratio);
    CHECK_FOR_PARAMETER("hsv_fract",hsv_fract);

    CHECK_FOR_INT_SETTING("et_ratio_function",et_ratio_function);
    CHECK_FOR_SETTING("et_ratio_p_a",et_ratio_p_a);
    CHECK_FOR_SETTING("et_ratio_p_b",et_ratio_p_b);
    CHECK_FOR_SETTING("et_ratio_p_c",et_ratio_p_c);
    CHECK_FOR_SETTING("et_ratio_p_d",et_ratio_p_d);

    CHECK_FOR_SETTING("et_ratio_p_max",et_ratio_p_max);
    CHECK_FOR_SETTING("et_ratio_p_min",et_ratio_p_min);
    CHECK_FOR_SETTING("et_ratio_p_mean",et_ratio_p_mean);
    CHECK_FOR_SETTING("et_ratio_p_std",et_ratio_p_std);
    CHECK_FOR_SETTING("et_ratio_p_k",et_ratio_p_k);

    CHECK_FOR_SETTING("et_ratio_c_max",et_ratio_c_max);
    CHECK_FOR_SETTING("et_ratio_c_min",et_ratio_c_min);
    CHECK_FOR_SETTING("et_ratio_c_mean",et_ratio_c_mean);
    CHECK_FOR_SETTING("et_ratio_c_std",et_ratio_c_std);
    CHECK_FOR_SETTING("et_ratio_c_k",et_ratio_c_k);

    CHECK_FOR_SETTING("et_ratio_c_a",et_ratio_c_a);
    CHECK_FOR_SETTING("et_ratio_c_b",et_ratio_c_b);

    CHECK_FOR_INT_SETTING("et_cluster_size",et_cluster_size);
    CHECK_FOR_INT_SETTING("et_num_clusters",et_num_clusters);
    CHECK_FOR_INT_SETTING("et_num_samples",et_num_samples);

    CHECK_FOR_PARAMETER("tcell_killing_rate",tcell_killing_rate);
    CHECK_FOR_PARAMETER("trm_motility",trm_motility);
    CHECK_FOR_PARAMETER("trm_levy_alpha",trm_levy_alpha);
    CHECK_FOR_PARAMETER("trm_levy_max",trm_levy_max);
    CHECK_FOR_INT_PARAMETER("max_tcell_div",max_tcell_div);
    CHECK_FOR_PARAMETER("tcell_density",tcell_density);
    CHECK_FOR_INT_PARAMETER("tcell_detection_range",tcell_detection_range);
    CHECK_FOR_INT_PARAMETER("pat_trm_activatedby",pat_trm_activatedby);
    CHECK_FOR_INT_PARAMETER("tcell_killing_target",tcell_killing_target);
#endif

#ifdef CYTOKINE_CODE
    CHECK_FOR_INT_SETTING("plot_cytokines",plot_cytokines);
    CHECK_FOR_INT_SETTING("cyto_act",cyto_act);
    CHECK_FOR_INT_SETTING("cyto_act_cyto",cyto_act_cyto);
    CHECK_FOR_INT_SETTING("cyto_act_prolif",cyto_act_prolif);

    CHECK_FOR_PARAMETER("max_cyt_days",max_cyt_days);
    CHECK_FOR_INT_PARAMETER("limit_cyt_spread",limit_cyt_spread);

    CHECK_FOR_PARAMETER("cyt_absorption_rate_trm",cyt_absorption_rate_trm);
    CHECK_FOR_PARAMETER("cyt_absorption_rate_cell",cyt_absorption_rate_cell);

    CHECK_FOR_PARAMETER("cyt_secretion_rate_tcell_mean",cyt_secretion_rate_tcell_mean);
    CHECK_FOR_PARAMETER("cyt_secretion_rate_infcell_mean",cyt_secretion_rate_infcell_mean);
    CHECK_FOR_PARAMETER("cyt_diff_rate",cyt_diff_rate);
    CHECK_FOR_PARAMETER("cyt_decay_rate",cyt_decay_rate);
    CHECK_FOR_PARAMETER("cyt_uptake",cyt_uptake);
    CHECK_FOR_PARAMETER("cyt_trm_ic50",cyt_trm_ic50);
    CHECK_FOR_PARAMETER("beta_ic50",beta_ic50);
    CHECK_FOR_PARAMETER("inf_ic50",inf_ic50);
    CHECK_FOR_PARAMETER("prod_ic50",prod_ic50);
    CHECK_FOR_PARAMETER("cyt_hill",cyt_hill);
    CHECK_FOR_PARAMETER("cyt_prot_expiration",cyt_prot_expiration);

    CHECK_FOR_INT_PARAMETER("cyt_effect_infectivity",cyt_effect_infectivity);
    CHECK_FOR_INT_PARAMETER("cyt_effect_virprod",cyt_effect_virprod);

    CHECK_FOR_PARAMETER("cyt_effect_lag",cyt_effect_lag);
    CHECK_FOR_PARAMETER("cyt_effect_uninfcelldeath",cyt_effect_uninfcelldeath);
    CHECK_FOR_PARAMETER("cyt_effect_infncelldeath",cyt_effect_infncelldeath);
    CHECK_FOR_PARAMETER("cyt_effect_infpcelldeath",cyt_effect_infpcelldeath);

    CHECK_FOR_INT_PARAMETER("diff_range",diff_range);
#endif

}

#ifdef TCELL_CODE
// Routine to read and process an E/T ratio input file (special option)
void read_et_ratios(const char *fname, int *num_read, double **et_array)
{
    FILE *et_file=NULL;
    if((et_file = fopen(fname,"r")) == NULL){
	fprintf(stderr,"Could not open E/T ratio file %s\nExiting!\n",fname);
	exit(1);
    }
    if (((*et_array) = (double *)malloc(MAX_ET_ENTRIES * sizeof(double))) == NULL) {
	fprintf(stderr,"Could not allocate E/T ratio array for %s\nExiting!\n",fname);
	exit(1);
    }
    char tmpline[MAX_LINE];
    int num_entries=0;

    while (fgets(tmpline, MAX_LINE-1,et_file) != NULL && num_entries < MAX_ET_ENTRIES) {
	if (sscanf(tmpline,"%lf",(*et_array)+num_entries) != 1) {
	    fprintf(stderr,"E/T ratio %d read failure (file=%s)\nExiting!\n",
		num_entries+1,fname);
	    exit(1);
	}
	num_entries++;
    }
    fclose(et_file);
    *num_read = num_entries;
}
#endif

// Usage printout (for -h option)
void usage(char *prog_name)
{

#ifdef TCELL_CODE
    fprintf(stderr,
	"Usage: %s [-h][-f <input_file>][-d][-e <e/t ratio file>][-o <dir>][-p|-E <prog>][-r][-s <seed>][-v]\n",prog_name);
#else
    fprintf(stderr,
	"Usage: %s [-h][-f <input_file>][-d][-o <dir>][-p <prog>][-r][-s <seed>][-v]\n",prog_name);
#endif
    fprintf(stderr,"\t-h = this help\n");
    fprintf(stderr,"\t-f = optional input file\n");
    fprintf(stderr,"\t-o = output directory for results.csv and cell state files\n");
    fprintf(stderr,"\t-p = Rscript program for cell state files (plotting)\n");
    fprintf(stderr,"\t-r = random number seed truely random\n");
    fprintf(stderr,"\t-s <seed> = random number seed for repeating run (unsigned)\n");
#ifdef TCELL_CODE
    fprintf(stderr,"\t-d = calc e/t ratios as exp w/ draw p & c values\n");
    fprintf(stderr,"\t-e = optional e/t ratio file\n");
    fprintf(stderr,"\t-E = Rscript program for E/T ratio files (plotting)\n");
#endif
}

// 0.  MAIN FUNCTION-only used if there is no wrapper file
int main(int argc, char **argv){
    
    long seed;
    int verbose = 0;
    char *inp_file=NULL;
#ifdef TCELL_CODE
    char *et_file=NULL;
    bool draw_et_ratios=false;
#endif

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
#ifdef TCELL_CODE
	if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"-D")){
	    draw_et_ratios=true;
	    if (et_file != NULL) {
		fprintf(stderr,
		    "Choose either -d or -e but not both! Exiting.\n");
		exit(1);
	    }
	}
	if (!strcmp(argv[i],"-e") && i +1 < argc) {
	    if (draw_et_ratios==true) {
		fprintf(stderr,
		    "Choose either -d or -e but not both! Exiting.\n");
		exit(1);
	    }
	    et_file = argv[i+1];
	    i++;
	}
	if (!strcmp(argv[i],"-E") && i +1 < argc) {
	    settings.et_process_cmd = argv[i+1];
	    i++;
	}
#endif
    }
    if (inp_file != NULL)
	read_input_file(inp_file, &settings, &parameters);

#ifdef TCELL_CODE
    if (et_file != NULL) {
	settings.et_file = et_file;
        settings.trm_init=DISTRIB;
    }
    if (draw_et_ratios)
        settings.trm_init=DRAWN;
#endif

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
	
#ifdef TCELL_CODE
	settings->et_file=NULL;
	settings->et_process_cmd=NULL;

	//TRM_HSV & TRM_BYST
	settings->trm_init=GAUSSIAN; //random or gaussian
	settings->trm_move_method=RANDOM; //random', 'directed' or 'levy'
	settings->trm_random_gaussian=1; //random center for gaussian distribution

	settings->add_cd4s=1;  // adding in CD4s as BYST
	settings->cd4_to_cd8_ratio=0.333; // this is 1:3 from histologies

	settings->cyto_act=0;
	settings->cyto_act_cyto=1;
	settings->cyto_act_prolif=1;

	settings->et_ratio_function=POLYNOMIAL;
	settings->et_cluster_size = 100;
	settings->et_num_clusters = MAX_CLUSTERS;
	settings->et_num_samples = 18;

	settings->et_ratio_p_max = 2;
	settings->et_ratio_p_min = 0.1;
	settings->et_ratio_p_mean = 0.293835;
	settings->et_ratio_p_std = 0.2770515;
	settings->et_ratio_p_k = 1.9028157;

	if (settings->et_ratio_function == LOGNORMAL) {
	    settings->et_ratio_p_a = -0.32;
	    settings->et_ratio_p_b = 0.9889;

	} else if (settings->et_ratio_function == POLYNOMIAL) {
	    settings->et_ratio_p_a = -0.0005;
	    settings->et_ratio_p_b = 0.0214;
	    settings->et_ratio_p_c = -0.2751;
	    settings->et_ratio_p_d = 1.3136;
	} else if (settings->et_ratio_function == POWER_LAW) {
	    settings->et_ratio_p_a = -0.0005;
	    settings->et_ratio_p_b = 0.0214;
	} else if (settings->et_ratio_function == EXPONENTIAL) {
	    settings->et_ratio_p_a = 1.1542;
	    settings->et_ratio_p_b = -0.14;
	}

	settings->et_ratio_c_max = -0.01;
	settings->et_ratio_c_min = -0.2;
	settings->et_ratio_c_mean = -0.09565;
	settings->et_ratio_c_std = 0.0307969;
	settings->et_ratio_c_k = 0.5686187;

	// C must be pulled then negated if not from a normal dist!
	if (settings->et_ratio_function == POWER_LAW) {
	    settings->et_ratio_c_a = 0.0265;
	    settings->et_ratio_c_b = -0.1034;
	} else if (settings->et_ratio_function == EXPONENTIAL) {
	    settings->et_ratio_c_a = 0.206;
	    settings->et_ratio_c_b = -0.04;
	}
#endif
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

#ifdef TCELL_CODE
	//T-cells
	parameters->trm_decay=0;	//Death rate of TRMs in absence of Ag (per site per day, default = 1E-3)
	parameters->trm_dbl=1;				//Proliferation (doubling) rate of TRMs (per cell per day)
	parameters->tcell_killing_rate=96;      //Rate of target killing by T-cells (by HSV-specific TRMs)
	parameters->trm_motility=5;					  //T-cell motility (jumps/movements per day)
	parameters->trm_levy_alpha=2.5;       //Levy exponent (should be between 1 and 3)
	parameters->trm_levy_max=10;          //Maximum jump per timestep (sites)

	//Probabilities, fractions, etc
	parameters->tcell_density=0.1;				//Density (fraction of total sites) of T-cells at the start of the simulation
	parameters->tcell_detection_range=3;	//Spatial range of detection of antigen from infected cell (n x n), default n = 3

	//Other
	parameters->pat_trm_activatedby=INFP;		//what activates a T-cell: INF,'trm', 'act_tem' or ANY
	parameters->tcell_killing_target=INFP;	//What are the T-cells killing? 'infp',INF,'all'. INF means both early and late
	
	parameters->hsv_fract=0.4;

	parameters->tem_trafficking=0;	//Rate of arrival of TEM from LN (per site per day)
	parameters->tem_exiting=100;	//Rate of exit of TEM from tissue
	parameters->tem_delay=1;      //Delay in the arrival of TEM (per day, default = 1 = 24h delay)
#endif

#ifdef CYTOKINE_CODE
	//Cytokine-related parameters
	parameters->cyt_secretion_rate_tcell_mean=0.05;     //Rate of cytokine_matrix secretion by T-cells, mean for normal dist (pg/cell/day)
	parameters->cyt_secretion_rate_infcell_mean=0.005;     //Rate of cytokine_matrix secretion by infected cells, mean for normal dist (pg/cell/day)
	parameters->cyt_diff_rate=65.5;      //Cytokine diffusion/spread rate in sites/day
	parameters->cyt_decay_rate=2.16;       //Rate of clearance of cytokine (/day)
	parameters->cyt_uptake=0.1;       //Rate of uptake of cytokine (/day)
	parameters->cyt_trm_ic50=7.406*1.25e-4;
	parameters->beta_ic50 = parameters->cyt_trm_ic50;
	parameters->inf_ic50 = parameters->cyt_trm_ic50;
	parameters->prod_ic50 = parameters->cyt_trm_ic50;
	parameters->cyt_effect_infectivity=0; //Do cytokines affect infectivity: yes (1) or no (0)
	parameters->cyt_effect_virprod=0; //Do cytokines affect virus production rate: yes (1) or no (0)
	parameters->cyt_effect_lag=0; //Max effect of cytokines on lag (multiplier, remember that lag here is a rate)
	parameters->cyt_effect_uninfcelldeath=48; //Max death rate of uninfected cells in the presence of cytokine_matrix, 0 means no effect
	parameters->cyt_effect_infncelldeath=48; //Max death rate of pre-productively infected cells cells in the presence of cytokine_matrix, 0 means no effect
	parameters->cyt_effect_infpcelldeath=48; //Max death rate of infected cells in the presence of cytokine_matrix, 0 means no effect
	parameters->limit_cyt_spread=0;

	parameters->cyt_prot_expiration=1;

#endif
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
	
#ifdef TCELL_CODE
	int **t_cell_state = (int **)my_malloc(settings->L * sizeof (int*));
	if (t_cell_state == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating t_cell Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    t_cell_state[i] = (int *)my_malloc(settings->L * sizeof (int));
	    if (t_cell_state[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating t_cell Matrix. Exiting.");
		exit(1);
	    }
	}
	dynamics->t_cell_state = t_cell_state;
	
	int **t_cell_potential = (int **)my_malloc(settings->L * sizeof (int*));
	if (t_cell_potential == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating t_cell_potential Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    t_cell_potential[i] = (int *)my_malloc(settings->L * sizeof (int));
	    if (t_cell_potential[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating t_cell_potential Matrix. Exiting.");
		exit(1);
	    }
	}
	dynamics->t_cell_potential = t_cell_potential;
#endif
	
#ifdef CYTOKINE_CODE
	double **t_cell_cyto_dur = (double **)my_malloc(settings->L * sizeof (double*));
	if (t_cell_cyto_dur == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating t_cell_cyto_dur Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    t_cell_cyto_dur[i] = (double *)my_malloc(settings->L * sizeof (double));
	    if (t_cell_cyto_dur[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating t_cell_cyto_dur Matrix. Exiting.");
		exit(1);
	    }
	}
	dynamics->t_cell_cyto_dur = t_cell_cyto_dur;
	
	double **cytokine_matrix = (double **)my_malloc(settings->L * sizeof (double *));
	if (cytokine_matrix == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating cytokine_matrix Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    cytokine_matrix[i] = (double *)my_malloc(settings->L * sizeof (double));
	    if (cytokine_matrix[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating cytokine_matrix Matrix. Exiting.");
		exit(1);
	    }
	}

	dynamics->cytokine_matrix = cytokine_matrix;
	double **delta_cytokine_matrix = (double **)my_malloc(settings->L * sizeof (double *));
	if (delta_cytokine_matrix == NULL)
	{
	    fprintf(stderr,"Encountered a problem creating delta_cytokine_matrix Matrix. Exiting.");
	    exit(1);
	}
	for (int i=0; i < settings->L; i++)
	{
	    delta_cytokine_matrix[i] = (double *)my_malloc(settings->L * sizeof (double));
	    if (delta_cytokine_matrix[i] == NULL)
	    {
		fprintf(stderr,"Encountered a problem creating delta_cytokine_matrix Matrix. Exiting.");
		exit(1);
	    }
	}
	dynamics->delta_cytokine_matrix = delta_cytokine_matrix;
#endif
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
	
#ifdef TCELL_CODE
	dynamics->gauss_offset=0;
	dynamics->pat_hsv=0;
	dynamics->pat_byst=0;
	dynamics->act_hsv=0;
	dynamics->act_byst=0;
	dynamics->dead_tcells=0;
	dynamics->avg_et_ratio=0;
	dynamics->perc_zeros=0;
	dynamics->et_ratio_3=0;
	dynamics->et_ratio_5=0;
	dynamics->et_ratio_7=0;
	dynamics->et_ratio_9=0;
	dynamics->et_ratio_11=0;
	dynamics->closestHSV=0;
	dynamics->closestHSV=0;
	dynamics->trm_act_time=1000;
	dynamics->trm_cyt_time=1000;
#endif
#ifdef CYTOKINE_CODE
	dynamics->cytokine=0;
#endif
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

#ifdef TCELL_CODE
	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++) {
		dynamics->t_cell_state[i][j]=EMPTY;
		dynamics->t_cell_potential[i][j]=parameters->max_tcell_div;
		dynamics->t_cell_cyto_dur[i][j]=parameters->max_cyt_days;
	    }
	dynamics->trm_kills=0;
	dynamics->trm_act_time=1000;
	dynamics->hsv_tems=0;
	dynamics->byst_tems=0;
	dynamics->tem_exits=0;
#endif

#ifdef CYTOKINE_CODE
	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++) {
	        dynamics->cytokine_matrix[i][j]=0;
	    }
	dynamics->cyto_kills=0;
	dynamics->trm_cyt_time=1000;
#endif

}

// clear deltas associated with a new time step (applied together at end of the step)
void clear_deltas(global_settings *settings,global_dynamics *dynamics) {
	for (int s=0; s < NUM_STRAINS; s++)
	    for (int i=0; i < settings->L; i++)
		for (int j=0; j < settings->L; j++)
		    dynamics->delta_viral_matrix[s][i][j]=0;
#ifdef CYTOKINE_CODE
	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++) {
	        dynamics->delta_cytokine_matrix[i][j]=0;
	    }
#endif
}

// apply deltas associated with a new time step 
void propagate_deltas(global_settings *settings,global_dynamics *dynamics) {
	for (int s=0; s < NUM_STRAINS; s++)
	    for (int i=0; i < settings->L; i++)
		for (int j=0; j < settings->L; j++)
		    dynamics->viral_matrix[s][i][j]=MAX(0,dynamics->viral_matrix[s][i][j] + dynamics->delta_viral_matrix[s][i][j]);

#ifdef TCELL_CODE
	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++)
		if (dynamics->t_cell_state[i][j] > 10)
		    dynamics->t_cell_state[i][j] -= 10;
#endif
#ifdef CYTOKINE_CODE
	for (int i=0; i < settings->L; i++)
	    for (int j=0; j < settings->L; j++)
	        dynamics->cytokine_matrix[i][j] =MAX(0,dynamics->cytokine_matrix[i][j] + dynamics->delta_cytokine_matrix[i][j]);
#endif

}


#ifdef CYTOKINE_CODE
//Effect of cytokine_matrix-based activation given conc of cytokine_matrix at site (pg/ml). In previous version, xmid and scal were estimated from Jia's experiments (see folder IFNG_jia) and multiplied to convert from units/ml to pg/cell
// This was changed to use an untuned version that allowed for varying IC50s to be tried (rather than tuned one)

// NOTE: the cytokine effect is a reductive effect (1 w/ no cytokine, 0 at highest
double new_cyt_effect(double dose, double cyt_ic50){
    return 1.0 /(1.0 + dose/cyt_ic50);
}
#endif

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

#ifdef TCELL_CODE
//PLACE T-CELLS (may be uniform, gaussian, from a file or distribution)
void place_tcells(global_settings *settings, global_parameters *parameters, global_dynamics *dynamics){
    int tot_t_cells = 0;
    int tot_t_cells9 = 0;
    int tot_t_cells25 = 0;
    int tot_t_cells49 = 0;
    int tot_t_cells81 = 0;
    int tot_t_cells121 = 0;

    int tot_cells = 0;
    int tot_cells9 = 0;
    int tot_cells25 = 0;
    int tot_cells49 = 0;
    int tot_cells81 = 0;
    int tot_cells121 = 0;
    int zeros=0;
    double minHSVdist = 100;
    double minTcelldist = 100;

    // if density is set to 0, do not place any T cells!
    if(parameters->tcell_density>0){

      // Option 1: draw et_ratio p & c from distributions then calculate 
      // for each of the clusters based on rank order
      if(settings->trm_init==DRAWN){
	double et_array[MAX_CLUSTERS];
	double et_dist_p;
	double et_dist_c;

	if (settings->et_ratio_function == NONE) {
	    et_dist_p = settings->et_ratio_p_mean;
	} else if (settings->et_ratio_function == WEIBULL) {
	    et_dist_p = gsl_ran_weibull(settings->ur,settings->et_ratio_p_mean,settings->et_ratio_p_k);
	} else if (settings->et_ratio_function == NORMAL) {
	    et_dist_p = settings->et_ratio_p_mean + gsl_ran_gaussian(settings->ur,settings->et_ratio_p_std);
	} else if (settings->et_ratio_function == LOGNORMAL) {
	    // p = a * ln(x) + b;
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    et_dist_p = settings->et_ratio_p_a * log(picked) + settings->et_ratio_p_b;
	} else if (settings->et_ratio_function == POLYNOMIAL) {
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    // p = a*x^3 + b*x^2 + c*x + d;
	    et_dist_p = settings->et_ratio_p_a * picked *picked * picked +
			settings->et_ratio_p_b * picked * picked +
			settings->et_ratio_p_c * picked + settings->et_ratio_p_d;
	} else if (settings->et_ratio_function == POWER_LAW) {
	    // p = a * x^b;
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    et_dist_p = settings->et_ratio_p_a * pow(picked,settings->et_ratio_p_b);
	} else if (settings->et_ratio_function == EXPONENTIAL) {
	    // p = a * e^bx;
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    et_dist_p = settings->et_ratio_p_a * pow(M_E,picked*settings->et_ratio_p_b);
	}
	// always pull c from a dist, (but with a negative sign if power law)
	if (settings->et_ratio_function == POWER_LAW) {
	    // c = - (a * x^b);
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    et_dist_c = -settings->et_ratio_c_a * pow(picked,settings->et_ratio_c_b);
	} else if (settings->et_ratio_function == EXPONENTIAL) {
	    // c = -a * e^bx;
	    int picked = 1 + gsl_rng_uniform_int(settings->ur,18);
	    et_dist_c = -settings->et_ratio_c_a * pow(M_E,picked*settings->et_ratio_c_b);
	} else
	    et_dist_c = (settings->et_ratio_c_mean + gsl_ran_gaussian(settings->ur,settings->et_ratio_c_std));

	if (et_dist_p > settings->et_ratio_p_max)
	    et_dist_p = settings->et_ratio_p_max;
	if (et_dist_p < settings->et_ratio_p_min)
	    et_dist_p = settings->et_ratio_p_min;

	if (et_dist_c > settings->et_ratio_c_max)
	    et_dist_c = settings->et_ratio_c_max;
	if (et_dist_c < settings->et_ratio_c_min)
	    et_dist_c = settings->et_ratio_c_min;

    fprintf(stderr,
	      "Calculating %d E/T ratios where p=%lf & c=%lf\n",
	      settings->et_num_clusters,et_dist_p,et_dist_c);

	double min_val=1;
	double max_val=0;
	double avg_val=0;

	int zero_clusters = settings->et_num_clusters*FRACT_W_ZEROS;
	int et_num_clusters=settings->et_num_clusters;
	int early_zero=0;

	for (int cluster=0; cluster < et_num_clusters; cluster++) {
	  if (cluster < et_num_clusters - zero_clusters)
	  {
	      et_array[cluster] = et_dist_p * pow(M_E,et_dist_c*(cluster+1));
	      if (early_zero == 0 && 
		  et_array[cluster] * settings->et_cluster_size < 1)
	      {
		et_array[cluster] = 0;
		early_zero = cluster;
	      }
	      if (et_array[cluster] < min_val)
		  min_val=et_array[cluster];
	      if (et_array[cluster] > max_val)
		  max_val=et_array[cluster];
	  }
	  else
	      et_array[cluster] = 0;

	  avg_val+=et_array[cluster];
	  fprintf(stderr,
	      "Rank=%d & e/t ratio=%lf\n",
	      cluster+1,et_array[cluster]);
	}
	if (early_zero > 0 && early_zero * 1.333 < et_num_clusters)
	    et_num_clusters = (int)(early_zero * 1.333);

	fprintf(stderr,"E/T ratio stats: %d regs min=%lf, max=%lf, avg=%lf\n",
	      et_num_clusters, min_val, max_val, avg_val / et_num_clusters);
	
	zeros = 0;
	int total_clusters = 0;
	int total_cd4s = 0;
	int dropped_cd4s = 0;
	int cluster_dim = (int)sqrt(settings->et_cluster_size);

	for (int cluster_row=0; cluster_row<(int)((settings->L/cluster_dim)+0.5);cluster_row++)
	    for (int cluster_col=0; cluster_col<(int)((settings->L/cluster_dim)+0.5);cluster_col++) {
		//double prob_zeros=gsl_rng_uniform(settings->ur);
		int cluster_num=gsl_rng_uniform_int(settings->ur,et_num_clusters);
		bool has_tcells=false;
		total_clusters++;
		for (int i=0; i < cluster_dim; i++)
		    for (int j=0; j < cluster_dim; j++) {
			int row = (cluster_dim)*cluster_row+i;
			int col = (cluster_dim)*cluster_col+j;
			double et_dice = gsl_rng_uniform(settings->ur);
			double cd4_dice = gsl_rng_uniform(settings->ur);
			bool place_tcell = false;
			double dist = sqrt(pow((row - settings->L/2),2.0)+pow((col - settings->L/2),2.0));
			// E/T ratio functions are based on just CD8s.  If
			// CD4s are to be included, add them to the cluster 
			// according to the cd4:cd8 ratio from the histologies
			if (settings->add_cd4s && 
			    et_dice<et_array[cluster_num] &&
			    cd4_dice < settings->cd4_to_cd8_ratio)
			{
			    // place a cd4 (bystander) in any open slot
			    // in the immediate 5x5 neighborhood
			    bool placed_cd4=false;
			    for (int k=-2; k <= 2; k++) {
				for (int l=-2; l <= 2; l++) {
				    if (row+k >= cluster_dim*cluster_row &&
					col+l >= cluster_dim*cluster_col &&
					row+k < cluster_dim*cluster_row+cluster_dim &&
					col+l < cluster_dim*cluster_col+cluster_dim &&
					(k != 0 || l != 0) &&
					dynamics->t_cell_state[row+k][col+l]==EMPTY)
				    {
					dynamics->t_cell_state[row+k][col+l]=PAT_BYST;
					if (dist < minTcelldist)
					    minTcelldist=dist;
					has_tcells=true;
					placed_cd4=true;
					total_cd4s++;
					tot_t_cells++;
					break;
				    }
				}
				if (placed_cd4)
				    break;
			    }
			    if (!placed_cd4)
				dropped_cd4s++;
			}

			    
			if (et_dice<et_array[cluster_num]*parameters->hsv_fract) {
			    // if a bystander (CD4) beat us here, 
			    // toss it out in favor of the HSV-specific one!
			    if (dynamics->t_cell_state[row][col]!=EMPTY) {
				place_tcell=true;
				dropped_cd4s++;
				total_cd4s--;
				tot_t_cells--;
			    }
			    dynamics->t_cell_state[row][col]=PAT_HSV;
			    if (dist < minHSVdist)
				minHSVdist=dist;
			    has_tcells=true;
			} else if (dynamics->t_cell_state[row][col]==EMPTY &&
				et_dice<et_array[cluster_num]){
			    dynamics->t_cell_state[row][col]=PAT_BYST;
			    if (dist < minTcelldist)
				minTcelldist=dist;
			    place_tcell=true;
			    has_tcells=true;
			}
			if (place_tcell) {
			    tot_t_cells++;
			    // immediate neighborhood of center point?
			    if (abs(row - settings->L/2) <= 1 && 
				abs(col - settings->L/2) <= 1)
				tot_t_cells9++;
			    if (abs(row - settings->L/2) <= 2 && 
				abs(col - settings->L/2) <= 2)
				tot_t_cells25++;
			    if (abs(row - settings->L/2) <= 3 && 
				abs(col - settings->L/2) <= 3)
				tot_t_cells49++;
			    if (abs(row - settings->L/2) <= 4 && 
				abs(col - settings->L/2) <= 4)
				tot_t_cells81++;
			    if (abs(row - settings->L/2) <= 5 && 
				abs(col - settings->L/2) <= 5)
				tot_t_cells121++;
			}
			if (abs(row - settings->L/2) <= 1 && 
			    abs(col - settings->L/2) <= 1 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells9++;
			if (abs(row - settings->L/2) <= 2 && 
			    abs(col - settings->L/2) <= 2 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells25++;
			if (abs(row - settings->L/2) <= 3 && 
			    abs(col - settings->L/2) <= 3 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells49++;
			if (abs(row - settings->L/2) <= 4 && 
			    abs(col - settings->L/2) <= 4 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells81++;
			if (abs(row - settings->L/2) <= 5 && 
			    abs(col - settings->L/2) <= 5 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells121++;
			if (dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells++;
		    }
		if (has_tcells == false)
		    zeros++;
	}
	dynamics->perc_zeros= 100*((double)zeros)/total_clusters;
        fprintf(stderr,"%d clusters had no effector cells (%lf%%)\n",
      		zeros,dynamics->perc_zeros);
	if (settings->add_cd4s && tot_t_cells > total_cd4s)
	    fprintf(stderr,"Added %d CD4s as bystander T cells (%lf%% of T cells) - dropped %d potential Cd4s\n",
		    total_cd4s, 100.0*total_cd4s/tot_t_cells,dropped_cd4s);

      // Option 2: read E/T ratios for clusters of <et_cluster_size> cells from file 
      } else if(settings->trm_init==DISTRIB){
	int clusters_read=0;
	double *et_array=NULL;
	read_et_ratios(settings->et_file,&clusters_read,&et_array);
	if (clusters_read == 0 || et_array==NULL){
	    fprintf(stderr,"Failed to read E/T ratios from file %s! Exiting.\n",
		settings->et_file);
	    exit(1);
	}
	//if (settings->verbose) {
	  fprintf(stderr,"Read %d E/T ratios from file %s\n",
		clusters_read,settings->et_file);
	  double min_val=1;
	  double max_val=0;
	  double avg_val=0;
	  for (int cluster=0; cluster < clusters_read; cluster++) {
	    if (et_array[cluster] < min_val)
		min_val=et_array[cluster];
	    if (et_array[cluster] > max_val)
		max_val=et_array[cluster];
	    avg_val+=et_array[cluster];
	  }
	  fprintf(stderr,"E/T ratio stats: min=%lf, max=%lf, avg=%lf\n",
		min_val, max_val, avg_val / clusters_read);
	//}
	
	zeros = 0;
	int total_clusters = 0;
	int cluster_dim = (int)sqrt(settings->et_cluster_size);
	for (int cluster_row=0; cluster_row<settings->L/cluster_dim;cluster_row++)
	    for (int cluster_col=0; cluster_col<settings->L/cluster_dim;cluster_col++) {
		int cluster_num=gsl_rng_uniform_int(settings->ur,clusters_read);
		total_clusters++;
		bool has_tcells=false;
		for (int i=0; i < cluster_dim; i++)
		    for (int j=0; j < cluster_dim; j++) {
			int row = cluster_dim*cluster_row+i;
			int col = cluster_dim*cluster_col+j;
			double et_dice = gsl_rng_uniform(settings->ur);
			bool place_tcell = false;
			double dist = sqrt(pow((row - settings->L/2),2.0)+pow((col - settings->L/2),2.0));

			if (et_dice<et_array[cluster_num]*parameters->hsv_fract) {
			    dynamics->t_cell_state[row][col]=PAT_HSV;
			    if (dist < minHSVdist)
				minHSVdist=dist;
			    place_tcell=true;
			} else if (et_dice<et_array[cluster_num]){
			    dynamics->t_cell_state[row][col]=PAT_BYST;
			    if (dist < minTcelldist)
				minTcelldist=dist;
			    place_tcell=true;
			}
			if (place_tcell) {
			    has_tcells=true;
			    tot_t_cells++;
			    // immediate neighborhood of center point?
			    if (abs(row - settings->L/2) <= 1 && 
				abs(col - settings->L/2) <= 1)
				tot_t_cells9++;
			    if (abs(row - settings->L/2) <= 2 && 
				abs(col - settings->L/2) <= 2)
				tot_t_cells25++;
			    if (abs(row - settings->L/2) <= 3 && 
				abs(col - settings->L/2) <= 3)
				tot_t_cells49++;
			    if (abs(row - settings->L/2) <= 4 && 
				abs(col - settings->L/2) <= 4)
				tot_t_cells81++;
			    if (abs(row - settings->L/2) <= 5 && 
				abs(col - settings->L/2) <= 5)
				tot_t_cells121++;
			}
			if (abs(row - settings->L/2) <= 1 && 
			    abs(col - settings->L/2) <= 1 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells9++;
			if (abs(row - settings->L/2) <= 2 && 
			    abs(col - settings->L/2) <= 2 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells25++;
			if (abs(row - settings->L/2) <= 3 && 
			    abs(col - settings->L/2) <= 3 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells49++;
			if (abs(row - settings->L/2) <= 4 && 
			    abs(col - settings->L/2) <= 4 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells81++;
			if (abs(row - settings->L/2) <= 5 && 
			    abs(col - settings->L/2) <= 5 &&
			    dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells121++;
			if (dynamics->cell_state[row][col]!=EMPTY)
			    tot_cells++;
		    }
		if (has_tcells == false)
		    zeros++;
	}
	dynamics->perc_zeros= 100*((double)zeros)/total_clusters;
        fprintf(stderr,"%d clusters had no effector cells (%lf%%)\n",
      		zeros,dynamics->perc_zeros);

      // Option 3: place all Trms randomly using tcell_density parameter
      } else if(settings->trm_init==RANDOM){
	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<settings->L;j++) {
		if (gsl_rng_uniform(settings->ur)<parameters->tcell_density) {
		    if (gsl_rng_uniform(settings->ur)<parameters->hsv_fract) {
			dynamics->t_cell_state[i][j]=PAT_HSV;
		    } else {
			dynamics->t_cell_state[i][j]=PAT_BYST;
		    }
		    // immediate neighborhood of center point?
		    if (abs(i - settings->L/2) <= 1 && 
			abs(j - settings->L/2) <= 1)
			tot_t_cells9++;
		    if (abs(i - settings->L/2) <= 2 && 
			abs(j - settings->L/2) <= 2)
			tot_t_cells25++;
		    if (abs(i - settings->L/2) <= 3 && 
			abs(j - settings->L/2) <= 3)
			tot_t_cells49++;
		    if (abs(i - settings->L/2) <= 4 && 
			abs(j - settings->L/2) <= 4)
			tot_t_cells81++;
		    if (abs(i - settings->L/2) <= 5 && 
			abs(j - settings->L/2) <= 5)
			tot_t_cells121++;
		    tot_t_cells++;
		}
		if (abs(i - settings->L/2) <= 1 && 
		    abs(j - settings->L/2) <= 1 &&
		    dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells9++;
		if (abs(i - settings->L/2) <= 2 && 
		    abs(j - settings->L/2) <= 2 &&
		    dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells25++;
		if (abs(i - settings->L/2) <= 3 && 
		    abs(j - settings->L/2) <= 3 &&
		    dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells49++;
		if (abs(i - settings->L/2) <= 4 && 
		    abs(j - settings->L/2) <= 4 &&
		    dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells81++;
		if (abs(i - settings->L/2) <= 5 && 
		    abs(j - settings->L/2) <= 5 &&
		    dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells121++;

		if (dynamics->cell_state[i][j]!=EMPTY)
		    tot_cells++;
	    }
      // Option 3: place all Trms using a gaussian with the tcell_density 
      }else if(settings->trm_init==GAUSSIAN){
	gsl_matrix *gausskern=compute_gaussian(3*settings->L/2,settings->L/3,parameters->tcell_density);
	double max_gauss=0;
	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<3*settings->L/2;j++)
	    {
		double value = gsl_matrix_get(gausskern,i,j);
		if (value >= max_gauss)
		    max_gauss = value;
	    }
	int xshift = 0;
	int yshift = 0;
	int source_x = 0;
	int source_y = 0;
	int target_x = 0;
	int target_y = 0;

	if (settings->trm_random_gaussian) {
	    // shift X in the range -L/4 to L/4
	    xshift = settings->L/4  - gsl_rng_uniform_int(settings->ur,settings->L/2);
	    yshift = settings->L/4  - gsl_rng_uniform_int(settings->ur,settings->L/2);
	    dynamics->gauss_offset=sqrt(xshift*xshift+yshift*yshift)/(settings->L/2);
	    fprintf(stderr,"T-cell gaussian centered at %d,%d dist = %3.2lf\n",
		settings->L/2+xshift,settings->L/2+yshift,dynamics->gauss_offset);
	}
	for (int i=0; i<3*settings->L/2;i++) {
	    source_x = i + xshift;
	    if (source_x >= 0 && source_x < 3*settings->L/2 ){
		target_y=0;
		for (int j=0; j<3*settings->L/2;j++) {
		    source_y = j + yshift;
		    if(source_y >= 0 && source_y < 3*settings->L/2 &&
			target_x >= 0 && target_x < settings->L &&
			target_y >= 0 && target_y < settings->L) {
			double value = gsl_matrix_get(gausskern,source_x,source_y);
			double new_value = value*parameters->tcell_density/max_gauss;
			if (gsl_rng_uniform(settings->ur)<new_value) {
			    if (gsl_rng_uniform(settings->ur)<parameters->hsv_fract) {
				dynamics->t_cell_state[target_x][target_y]=PAT_HSV;
			    } else {
				dynamics->t_cell_state[target_x][target_y]=PAT_BYST;
			    }
			    // immediate neighborhood of center point?
			    if (abs(target_x - settings->L/2) <= 1 && 
				abs(target_y - settings->L/2) <= 1)
				tot_t_cells9++;
			    if (abs(target_x - settings->L/2) <= 2 && 
				abs(target_y - settings->L/2) <= 2)
				tot_t_cells25++;
			    if (abs(target_x - settings->L/2) <= 3 && 
				abs(target_y - settings->L/2) <= 3)
				tot_t_cells49++;
			    if (abs(target_x - settings->L/2) <= 4 && 
				abs(target_y - settings->L/2) <= 4)
				tot_t_cells81++;
			    if (abs(target_x - settings->L/2) <= 5 && 
				abs(target_y - settings->L/2) <= 5)
				tot_t_cells121++;
			    tot_t_cells++;
			}
			if (abs(target_x - settings->L/2) <= 1 && 
			    abs(target_y - settings->L/2) <= 1 &&
			    dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells9++;
			if (abs(target_x - settings->L/2) <= 2 && 
			    abs(target_y - settings->L/2) <= 2 &&
			    dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells25++;
			if (abs(target_x - settings->L/2) <= 3 && 
			    abs(target_y - settings->L/2) <= 3 &&
			    dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells49++;
			if (abs(target_x - settings->L/2) <= 4 && 
			    abs(target_y - settings->L/2) <= 4 &&
			    dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells81++;
			if (abs(target_x - settings->L/2) <= 5 && 
			    abs(target_y - settings->L/2) <= 5 &&
			    dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells121++;

			if (dynamics->cell_state[target_x][target_y]!=EMPTY)
			    tot_cells++;
			target_y=target_y+1;
		    }
		}
	    }
	    if(source_x >= 0 && source_x < 3*settings->L/2)
		target_x=target_x+1;
        }
      }
      // regardless of placement strategy, spit out starting E/T ratios 
      // for the various NxN neighborhoods around the infection seed point
      dynamics->et_ratio_3= (double)tot_t_cells9/tot_cells9;
      dynamics->et_ratio_5= (double)tot_t_cells25/tot_cells25;
      dynamics->et_ratio_7= (double)tot_t_cells49/tot_cells49;
      dynamics->et_ratio_9= (double)tot_t_cells81/tot_cells81;
      dynamics->et_ratio_11= (double)tot_t_cells121/tot_cells121;
      dynamics->closestHSV=minHSVdist;
      dynamics->closestTcell=minTcelldist;
      fprintf(stdout,"3x3 E/T ratio is %e (%d/%d)\n",dynamics->et_ratio_3,
      		tot_t_cells9,tot_cells9);
      fprintf(stdout,"5x5 E/T ratio is %e (%d/%d)\n",dynamics->et_ratio_5,
      		tot_t_cells25,tot_cells25);
      fprintf(stdout,"7x7 E/T ratio is %e (%d/%d)\n",dynamics->et_ratio_7,
      		tot_t_cells49,tot_cells49);
      fprintf(stdout,"9x9 E/T ratio is %e (%d/%d)\n",dynamics->et_ratio_9,
      		tot_t_cells81,tot_cells81);
      fprintf(stdout,"11x11 E/T ratio is %e (%d/%d)\n",dynamics->et_ratio_11,
      		tot_t_cells121,tot_cells121);

      dynamics->avg_et_ratio=(double)tot_t_cells/tot_cells;
      fprintf(stdout,"Avg E/T ratio is %e (%d/%d)\n",dynamics->avg_et_ratio,
      		tot_t_cells,tot_cells);
   }
}

// supports option to write out E/T ratios (and use with special R script)
// This allows the plotting of the ratios at various time points as the sim runs
void save_et_ratios(global_settings *settings,global_dynamics *dynamics){
    int total_clusters = 0;
    double cluster_dim = 10;
    double min_val=1;
    double max_val=0;
    double avg_val=0;
    int zeros=0;

    string output_file;

    string output_dir;
    FILE * outF= NULL;

    char temp_str[100];

    sprintf(temp_str,"et_ratio_%d_%4.3lf.csv",
	dynamics->runnum+1,dynamics->time);

    if (settings->dirname != NULL) {
	output_file=settings->dirname;
	output_file+="/";
	output_file+=temp_str;
    }
    else
	output_file=temp_str;

    if((outF = fopen(output_file.c_str(),"wt")) == NULL){
	fprintf(stderr,"Could not open data file %s\n",output_file.c_str());
	exit(1);
    }
    fprintf(outF,"cluster,e/t ratio\n");

    for (int cluster_row=0; cluster_row<(int)(settings->L/cluster_dim);cluster_row++)
	for (int cluster_col=0; cluster_col<(int)(settings->L/cluster_dim);cluster_col++) {
	    total_clusters++;
	    int cluster_cells=0;
	    int cluster_tcells=0;
	    for (int i=0; i < (int)cluster_dim; i++)
		for (int j=0; j < (int)cluster_dim; j++) {
		    int row = ((int)cluster_dim)*cluster_row+i;
		    int col = ((int)cluster_dim)*cluster_col+j;

		    if (dynamics->cell_state[row][col] < DEAD &&
		        dynamics->cell_state[row][col] != EMPTY) {
			cluster_cells++;
		    }
		    if (dynamics->t_cell_state[row][col] != DEAD &&
		        dynamics->t_cell_state[row][col] != EMPTY) {
			cluster_tcells++;
		    }
		}
	
	    double this_et_ratio = ((double)cluster_tcells)/cluster_cells;
	    if (this_et_ratio == 0)
		zeros++;
	    if (this_et_ratio < min_val)
		min_val=this_et_ratio;
	    if (this_et_ratio > max_val)
		max_val=this_et_ratio;
	    avg_val+=this_et_ratio;
	    fprintf(outF,"%d,%4.3lf\n",
		total_clusters,this_et_ratio);
	  }
    //fprintf(stderr,"E/T ratio stats: zeros=%lf%%, min=%lf, max=%lf, avg=%lf\n",
//	100*((double)zeros)/total_clusters,min_val, max_val, avg_val / total_clusters);
    fclose(outF);
    if (settings->et_process_cmd != NULL) {
	string command;
	command = "Rscript ";
	command += settings->et_process_cmd;
	command += " ";
	if (settings->dirname != NULL) {
	    command+=settings->dirname;
	}
	else
	    command+=".";

	command+= "/ ";

	sprintf(temp_str,"%d %4.3lf",
	    dynamics->runnum+1,dynamics->time);

	command+= temp_str;

	//if (settings->verbose)
	  fprintf(stderr,"Launching command: %s\n",command.c_str());

	int ret = system(command.c_str());
	if (WIFSIGNALED(ret) &&
	    (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT)) {
	    fprintf(stderr,"Rscript aborted.  Exiting!\n");
	    exit(1);
	}
	remove(output_file.c_str());
    }
}

#endif

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

#ifdef CYTOKINE_CODE
//DIFFUSION of cytokines
void diffusion_cyt(global_settings *settings,int i, int j,gsl_matrix *gaussian_mask,global_parameters *parameters, global_dynamics *dynamics){
    double burst = dynamics->cytokine_matrix[i][j];
    //virus at focal site [i][j] will be spread across diffusion range
    // note: some will be placed at [i][j] as well
    dynamics->delta_cytokine_matrix[i][j]-=dynamics->cytokine_matrix[i][j]; 
    distrib_progeny(settings,parameters,dynamics,CYTOKINE,0,i,j,burst,gaussian_mask); 
}

bool close_to_plaq(int i,int j, global_settings *settings,global_parameters *parameters, global_dynamics *dynamics){
    int temp2=parameters->limit_cyt_spread;
    int imin=MAX(0,i-temp2);
    int imax=MIN(i+temp2,settings->L-1); //deal with edges
    int jmin=MAX(0,j-temp2);
    int jmax=MIN(j+temp2,settings->L-1);
    for (int row=imin; row <= imax; row++)
	for (int col=jmin; col <= jmax; col++) {
	    if (dynamics->cell_state[row][col] == INFP_A ||
		dynamics->cell_state[row][col] == INFP_B ||
		dynamics->cell_state[row][col] == INFP_ABn ||
		dynamics->cell_state[row][col] == INFP_ABr) {
		double dist_to_plaq=sqrt((double)((i-row)*(i-row)+(j-col)*(j-col)));
		if (dist_to_plaq <= parameters->limit_cyt_spread)
		    return true;
	    }
    }
    return false;
}
#endif
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
#ifdef CYTOKINE_CODE
	    } else if (parameters->limit_cyt_spread== 0 ||
			close_to_plaq(row,col,settings,parameters,dynamics)){
		double this_batch = 
		    burst * gsl_matrix_get(gaussian_mask,i1,j1); 
		dynamics->delta_cytokine_matrix[row][col]+=this_batch;
#endif
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
#ifdef TCELL_CODE
      }else if(what==PAT_HSV){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_HSV)
		    return 1;
	    }
      }else if(what==ACT_HSV){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == ACT_HSV)
		    return 1;
	    }
      }else if(what==PAT_BYST){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_BYST)
		    return 1;
	    }
      }else if(what==ACT_BYST){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == ACT_BYST)
		    return 1;
	    }
      }else if(what==T_CELL){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_HSV ||
		    dynamics->t_cell_state[row][col] == ACT_HSV ||
		    dynamics->t_cell_state[row][col] == PAT_BYST ||
		    dynamics->t_cell_state[row][col] == ACT_BYST)
		    return 1;
	    }
      }else if(what==T_KILL){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		// Only activated HSV+ Trms kill infected cells in this version!
		if (dynamics->t_cell_state[row][col] == ACT_HSV)
		    return 1;
	    }
#endif
      }else if(what==UNINF){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cell_state[row][col] == SUSCEPTIBLE ||
		    dynamics->cell_state[row][col] == UNPROTECTABLE)
		    return 1;
	    }
#ifdef CYTOKINE_CODE
      }else if(what==CYTOKINE){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cytokine_matrix[row][col] > 0)
		    return 1;
	    }
#endif
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
#ifdef TCELL_CODE
      }else if(what==PAT_HSV){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_HSV)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_HSV
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==PAT_BYST){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_BYST)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->t_cell_state[row][col] == PAT_BYST)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==ACT_HSV){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == ACT_HSV)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == ACT_HSV
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==ACT_BYST){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == ACT_BYST)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->t_cell_state[row][col] == ACT_BYST)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==T_CELL){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->t_cell_state[row][col] == PAT_HSV ||
		    dynamics->t_cell_state[row][col] == ACT_HSV ||
		    dynamics->t_cell_state[row][col] == PAT_BYST ||
		    dynamics->t_cell_state[row][col] == ACT_BYST)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if ((dynamics->t_cell_state[row][col] == PAT_HSV ||
		    dynamics->t_cell_state[row][col] == ACT_HSV ||
		    dynamics->t_cell_state[row][col] == PAT_BYST ||
		    dynamics->t_cell_state[row][col] == ACT_BYST)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
      }else if(what==T_KILL){
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		// Only activated HSV+ Trms kill infected cells in this version!
		if (dynamics->t_cell_state[row][col] == ACT_HSV)
		    det_inf++;
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	int instance = gsl_rng_uniform_int(settings->ur,det_inf);
	int current=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
/*
		if ((dynamics->t_cell_state[row][col] == PAT_HSV ||
*/
		if ((dynamics->t_cell_state[row][col] == ACT_HSV)
			&& current++ == instance)
		{
		    coords[0] = row;
		    coords[1] = col;
		    return(det_inf);
		}
	    }
#endif
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
#ifdef CYTOKINE_CODE
      }else if(what==CYTOKINE){
	int highest_row=imin;
	int highest_col=jmin;
	double highest_cyto=0;
	for (int row=imin; row <= imax; row++)
	    for (int col=jmin; col <= jmax; col++)
	    {
		if (dynamics->cytokine_matrix[row][col] > 0) {
		    det_inf++;
		    if (dynamics->cytokine_matrix[row][col] > highest_cyto) {
			highest_cyto=dynamics->cytokine_matrix[row][col];
			highest_row=row;
			highest_col=col;
		    }
		}
	    }
	if (return_pos==SUM || det_inf == 0) return det_inf;
	coords[0] = highest_row;
	coords[1] = highest_col;
	return(det_inf);
#endif
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

#ifdef TCELL_CODE
//T-CELL MOVEMENT 
//Based on settings->trm_move_method, perform appropriate TRM patrolling function
void move_tcell(int i, int j,global_settings *settings, global_parameters *parameters, global_dynamics *dynamics){
  bool det_neighbor = false;
  bool det_cyt= false;
  bool det_inf= false;
  int nbr_cyt[2];
  int nbr_inf[2];
  int newpos[2];

  if(settings->trm_move_method==RANDOM){
    //Random walk
      det_neighbor=pick_nbr(settings,i,j,newpos);
      if (!det_neighbor)
	return;
  }else if(settings->trm_move_method==DIRECTED){
      det_cyt=detect_inf(settings,dynamics,parameters->tcell_detection_range,i,j,CYTOKINE,COORDS,nbr_cyt);
      det_inf=detect_inf(settings,dynamics,parameters->tcell_detection_range,i,j,INFP,COORDS,nbr_inf);
      if(!det_inf && !det_cyt){ 
	//random if there's no cytokine_matrix or infected cell in nbhd
	det_neighbor=pick_nbr(settings,i,j,newpos);
	if (!det_neighbor)
	    return;
      } else {
      //move towards cytokine or antigen (preference to infected cell)
	if(det_inf){ 
	    if (nbr_inf[0] <= i)
		newpos[0]=MAX(0,i-1);
	    else
		newpos[0]=MIN(i+1,settings->L-1);

	    if (nbr_inf[1] <= j)
		newpos[1]=MAX(0,j-1);
	    else
		newpos[1]=MIN(j+1,settings->L-1);
        }else{
	    if (nbr_cyt[0] <= i)
		newpos[0]=MAX(0,i-1);
	    else
		newpos[0]=MIN(i+1,settings->L-1);

	    if (nbr_cyt[1] <= j)
		newpos[1]=MAX(0,j-1);
	    else
		newpos[1]=MIN(j+1,settings->L-1);
        }
      }
  }else if(settings->trm_move_method==LEVY){
    //(Levy flight)
    //Using implementation on stack exchange and assuming min jump=1 (see levy_tests.R, implementation 2b)
    //http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870
      double theta=gsl_rng_uniform(settings->ur)*2*M_PI; //pick a random direction
      double l=pow(gsl_rng_uniform(settings->ur),(-1/parameters->trm_levy_alpha)); //Levy displacement, assumes min=1
      //l[l>parameters->trm_levy_max]=parameters->trm_levy_max; //truncate at max
      newpos[0]=MIN(MAX(0,round(i+l*sin(theta))),settings->L-1);
      newpos[1]=MIN(MAX(0,round(j+l*cos(theta))),settings->L-1);
  }
  //Check if outside grid or if occupied. If so, no move
  if (dynamics->t_cell_state[newpos[0]][newpos[1]] == EMPTY ||
	dynamics->t_cell_state[newpos[0]][newpos[1]] == DEAD){
      // add 10 to code to ensure no double moves!
      dynamics->t_cell_state[newpos[0]][newpos[1]]=dynamics->t_cell_state[i][j]+10;
      dynamics->t_cell_state[i][j]=0;
      dynamics->t_cell_potential[newpos[0]][newpos[1]]=dynamics->t_cell_potential[i][j];
      dynamics->t_cell_potential[i][j]=0;
      dynamics->t_cell_cyto_dur[newpos[0]][newpos[1]]=dynamics->t_cell_cyto_dur[i][j];
      dynamics->t_cell_cyto_dur[i][j]=0;
  }
  
}

//CELL DIVISION of T-cells
void tcell_div(int i, int j,global_settings *settings, global_parameters *parameters, global_dynamics *dynamics){
	int cell_nbrs=3;
	int imin1=MAX(0,i-(cell_nbrs-1)/2);
	int imax1=MIN(i+(cell_nbrs-1)/2,settings->L-1);
	int jmin1=MAX(0,j-(cell_nbrs-1)/2);
	int jmax1=MIN(j+(cell_nbrs-1)/2,settings->L-1);

	// if there are any neighboring cells w/o a resident t-cell
	// then, divide this one and occupy one of those cells
	// 1st count empty spots, then randomly pick one
	int vacants=0;
        for (int row=imin1; row <= imax1; row++)
	    for (int col=jmin1; col <= jmax1; col++)
	      if (dynamics->t_cell_state[row][col]==0)
		vacants++;
	
	if(vacants>0){
	    int ind=gsl_rng_uniform_int(settings->ur,vacants);
	    int curr=0;
	    for (int row=imin1; row <= imax1; row++)
		for (int col=jmin1; col <= jmax1; col++)
		  if (dynamics->t_cell_state[row][col]==0)
		  {
		    if (curr == ind)
		    {
			dynamics->t_cell_state[row][col]=dynamics->t_cell_state[i][j];
			dynamics->t_cell_potential[i][j]--;
			dynamics->t_cell_potential[row][col]=dynamics->t_cell_potential[i][j];
			dynamics->t_cell_cyto_dur[row][col]=parameters->max_cyt_days;
			return;
		    }
		    curr++;
		  }
	}
}
#endif

//DECAY of free virus
void viral_decay(int s,int i, int j,global_parameters *parameters, global_dynamics *dynamics){
    dynamics->viral_matrix[s][i][j]=MAX(0,dynamics->viral_matrix[s][i][j]*exp(-parameters->freevirus_decay/parameters->max_rate));
}

#ifdef CYTOKINE_CODE
//DECAY of cytokine_matrix
void cyt_decay(int i, int j,global_parameters *parameters, global_dynamics *dynamics){
  dynamics->cytokine_matrix[i][j]=MAX(0,dynamics->cytokine_matrix[i][j]*exp(-parameters->cyt_decay_rate/parameters->max_rate));
}

void cyt_uptake(int i, int j,global_parameters *parameters, global_dynamics *dynamics){
  dynamics->cytokine_matrix[i][j]=MAX(0,dynamics->cytokine_matrix[i][j]*exp(-parameters->cyt_uptake/parameters->max_rate));
}
#endif

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

#ifdef CYTOKINE_CODE
    int cdist_range = MAX(parameters->diff_range,1+2*parameters->cyt_diff_rate/parameters->max_rate);
    gsl_matrix *gaussian_mask_cyt_diff=compute_gaussian(cdist_range,1.2,1); 
#endif
     
    place_host(settings,parameters,dynamics); //Place host cells

    if(settings->place_virus_type!=NO_VIRUS) 
	    place_virus(settings,dynamics); //Place virus	

#ifdef TCELL_CODE
    place_tcells(settings,parameters,dynamics); //Place T-cells
#endif
    
    dynamics->num_updates=0;
    compute_totals(settings,parameters,dynamics);
    output_results(runnum,settings,parameters,dynamics,false);
    if (settings->save_plot_data || settings->save_state_files)
	dump_data(settings,parameters,dynamics);
    
#ifdef TCELL_CODE
    if (settings->save_ratio_data)
	save_et_ratios(settings,dynamics);
#endif

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
		    	
#ifdef TCELL_CODE
	      ////T-cells - these are antigen-specific cells
		//// 0: no T-cell present
		//// PAT_HSV (patrolling HSV-specific)
		//// ACT_HSV (active HSV-specific)
		//// PAT_BYST (patrolling nonHSV-specific)
		//// ACT_BYST (active nonHSV-specific)
		//// 5: dead T-cell
		////--------------------------------------------
    
	      bool detected_antigen =
		    detect_inf(settings,dynamics,parameters->tcell_detection_range,i,j,
			    parameters->pat_trm_activatedby,CHECK,coords);

	      bool detected_cytokines = detect_inf(settings,dynamics,
		    parameters->tcell_detection_range,
		    i,j, CYTOKINE,CHECK,coords);

	      // in the "enhanced state" t-cells still move, but also 
	      // proliferate and produce their own cytokines!
		if(dynamics->t_cell_state[i][j]==EMPTY||
			dynamics->t_cell_state[i][j]==DEAD){

		//TEM arrives from LN (trafficking) if enough time has elapsed 
		// and probability met (based on rate)
		    if(dynamics->time>(parameters->tem_delay)&&
		       dice_t<=parameters->tem_trafficking/parameters->max_rate) {
			    if (gsl_rng_uniform(settings->ur) < parameters->hsv_fract)
			    {
				dynamics->t_cell_state[i][j]=PAT_HSV;
				dynamics->hsv_tems++;
			    }
			    else
			    {
				dynamics->t_cell_state[i][j]=PAT_BYST;
				dynamics->byst_tems++;
			    }
			}
	      }
	      // in the "enhanced state" t-cells still move, but also 
	      // proliferate and produce their own cytokines!
	      else if (dynamics->t_cell_state[i][j] == ACT_HSV ||
		  dynamics->t_cell_state[i][j] == ACT_BYST) {

		  //Event 1a.0: produce cytokine_matrix 
		  if(dynamics->t_cell_cyto_dur[i][j] > 0) {
		      dynamics->delta_cytokine_matrix[i][j]+=
			MAX(0, gsl_ran_weibull(settings->ur,parameters->cyt_secretion_rate_tcell_mean/parameters->max_rate,5.0));
		      dynamics->t_cell_cyto_dur[i][j] -= 1.0/parameters->max_rate;
		  }
		
		  //Event 1a.1: proliferate
		  if((dynamics->t_cell_state[i][j] == ACT_HSV ||
		      dynamics->t_cell_state[i][j] == ACT_BYST) &&
		      (parameters->max_tcell_div == 0 || 
			dynamics->t_cell_potential[i][j] > 0) &&
		      dice_t<=parameters->trm_dbl/parameters->max_rate){
			tcell_div(i,j,settings,parameters,dynamics);
		  }
	      }
	      //Case a: Antigen present
	      //  HSV + Trms activate!
	      else if (dynamics->t_cell_state[i][j] == PAT_HSV &&
		    detected_antigen) {
		    dynamics->t_cell_state[i][j] = ACT_HSV;
		    if (dynamics->time < dynamics->trm_act_time){
			dynamics->trm_act_time=dynamics->time;
			fprintf(stderr,"HSV+ activation at t=%lf\n",
			    dynamics->time);
		    }
	      } else if (detected_cytokines && settings->cyto_act) {
		    if (dynamics->t_cell_state[i][j] == PAT_HSV ||
			dynamics->t_cell_state[i][j] == PAT_BYST ) {

		      //Event 1a.0: produce cytokine_matrix 
		      cytokine_factor = new_cyt_effect(dynamics->cytokine_matrix[i][j],
				parameters->cyt_trm_ic50);

		      if(settings->cyto_act_cyto && 
			dynamics->t_cell_cyto_dur[i][j] > 0) {
			  dynamics->delta_cytokine_matrix[i][j]+=
			    MAX(0, gsl_ran_weibull(settings->ur,
				(1-cytokine_factor)*parameters->cyt_secretion_rate_tcell_mean/parameters->max_rate,5.0));
			  dynamics->t_cell_cyto_dur[i][j] -= 1.0/parameters->max_rate;
		      }
		      //Event 1b.0: replicate when cytokines are present

		      if( settings->cyto_act_prolif &&
			   (parameters->max_tcell_div == 0 || 
			    dynamics->t_cell_potential[i][j] > 0) &&
			  dice_t<=(1-cytokine_factor)*parameters->trm_dbl/parameters->max_rate){
			    tcell_div(i,j,settings,parameters,dynamics);
		      }

		    }
		      //Event 1c.0: continue to patroll
		    if (dynamics->t_cell_state[i][j] != ACT_HSV || !detected_antigen){
			double move_probability=gsl_rng_uniform(settings->ur);

			if (move_probability <= (parameters->trm_motility)/parameters->max_rate)
			    move_tcell(i,j,settings,parameters,dynamics);
		    }
	      //Case b: Antigen and cytokines absent
	      }else{
		
		//Event 1b.0: leave (if Tem)
		if(parameters->tem_trafficking > 0 &&
		   (dynamics->t_cell_state[i][j]==PAT_HSV ||
		   dynamics->t_cell_state[i][j]==PAT_BYST) &&
		   dice_t<=parameters->tem_exiting/parameters->max_rate){
		      dynamics->t_cell_state[i][j]=EMPTY;
		      dynamics->tem_exits++;
		}
		//Event 1b.1: decay/die
		else if(dice_t<=parameters->trm_decay/parameters->max_rate){
		  dynamics->t_cell_state[i][j]=DEAD;
		  
		//Event 1b.2: move -All patrolling phenotypes, activated bystanders
		// and HSV+ Trms lacking antigen to combat
		}else {
		    if (dynamics->t_cell_state[i][j] == PAT_HSV ||
			dynamics->t_cell_state[i][j] == PAT_BYST ||
			(dynamics->t_cell_state[i][j] == ACT_HSV && !detected_antigen) ||
			dynamics->t_cell_state[i][j] == ACT_BYST) {
			double move_probability=gsl_rng_uniform(settings->ur);

			if (move_probability <= (parameters->trm_motility)/parameters->max_rate)
			    move_tcell(i,j,settings,parameters,dynamics);
		    }
		  }
		}
#endif
					      
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

#ifdef CYTOKINE_CODE
		    if(dynamics->cell_state[i][j]!=UNPROTECTABLE &&
			parameters->cyt_effect_infectivity){
			double dice_prot=gsl_rng_uniform(settings->ur);
			if (dynamics->cytokine_matrix[i][j] > 0 &&
			    dice_prot <= parameters->cyt_prot_expiration/parameters->max_rate) {
				dynamics->cell_state[i][j]=UNPROTECTABLE;
				cytokine_factor = 1;
			} else {
			    if (dynamics->cytokine_matrix[i][j] > 0) {
				cytokine_factor = new_cyt_effect(dynamics->cytokine_matrix[i][j],
					parameters->beta_ic50);
			    } else {
				cytokine_factor = 1;
			    }
			}
		    } else {
			cytokine_factor = 1;
		    }
#endif
		
		    //Event 1.1: infection due to free virus at the site
		    bool virus_present=false;
		    bool new_infection=false;
		    for (int s=0; s < NUM_STRAINS; s++)
		    {
		      if(dynamics->viral_matrix[s][i][j]>0){
			virus_present=true;

#ifdef CYTOKINE_CODE
			//with cytokine effects
			if(parameters->cyt_effect_infectivity==1) {
			  if(dice<=(
				(parameters->infectivity_free[s]*cytokine_factor)*
				dynamics->viral_matrix[s][i][j])/parameters->max_rate){

			      set_infection_type(s, i, j, settings, parameters, dynamics, &new_infection);
			//without cytokine effects
			    } else if(dice<=(
				parameters->infectivity_free[s]*
				dynamics->viral_matrix[s][i][j])/
				parameters->max_rate){
			      infections_avoided++;
			    }
			}else 
#endif
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
#ifdef CYTOKINE_CODE
			    (parameters->cyt_effect_uninfcelldeath==0||
			    dynamics->cytokine_matrix[i][j] == 0)&&
#endif
			    
			     dice<=(parameters->uninf_cell_death)/parameters->max_rate){
			  dynamics->cell_state[i][j]=DEAD;
			}

#ifdef CYTOKINE_CODE
		      //with cytokine effects
		    }else if(parameters->cyt_effect_uninfcelldeath>0&&
			dice<=(((1-cytokine_factor)*
			(parameters->cyt_effect_uninfcelldeath-
			 parameters->uninf_cell_death) + 
			parameters->uninf_cell_death)/parameters->max_rate)){

		      dynamics->cell_state[i][j]=5; 
		      dynamics->cyto_kills++;
#endif
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

#ifdef CYTOKINE_CODE
			dynamics->delta_cytokine_matrix[i][j]+=
			    MAX(0, gsl_ran_weibull(settings->ur,parameters->cyt_secretion_rate_infcell_mean/parameters->max_rate,5.0));

		      //Event 2.1: become a virus-producing cell
			  //without cytokine effects
#ifdef CYTOKINE_CODE
		      if(parameters->cyt_effect_lag==0&&
			  dice_cell_lag<=(parameters->viral_lag)/parameters->max_rate){
#else
		      if(dice_cell_lag<=(parameters->viral_lag)/parameters->max_rate){
#endif
			if (dynamics->cell_state[i][j] == INFN_A) dynamics->cell_state[i][j]=INFP_A; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_B) dynamics->cell_state[i][j]=INFP_B; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABn) dynamics->cell_state[i][j]=INFP_ABn; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABr) dynamics->cell_state[i][j]=INFP_ABr; //becomes virus-producing infected cell
		      
			  //with cytokine effects
#ifdef CYTOKINE_CODE
		      }else if(parameters->cyt_effect_lag>0&&
			  dice_cell_lag<=(parameters->viral_lag*cytokine_factor)/parameters->max_rate){
			if (dynamics->cell_state[i][j] == INFN_A) dynamics->cell_state[i][j]=INFP_A; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_B) dynamics->cell_state[i][j]=INFP_B; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABn) dynamics->cell_state[i][j]=INFP_ABn; //becomes virus-producing infected cell
			else if (dynamics->cell_state[i][j] == INFN_ABr) dynamics->cell_state[i][j]=INFP_ABr; //becomes virus-producing infected cell
#endif
				    
			    //Event 2.2: dies 
		      
	    
		      }else {
#ifdef CYTOKINE_CODE
			if (dynamics->cytokine_matrix[i][j] > 0) {
			    cytokine_factor = new_cyt_effect(dynamics->cytokine_matrix[i][j],
			      parameters->inf_ic50);
			} else {
			    cytokine_factor = 1;
			}

			  //without cytokine effects
			if((parameters->cyt_effect_infncelldeath==0 ||
				dynamics->cytokine_matrix[i][j] == 0)&&
				dice_cell_death<=(inf_cell_death)/parameters->max_rate){
#else
			if(dice_cell_death<=(inf_cell_death)/parameters->max_rate){
#endif
			    dynamics->cell_state[i][j]=DEAD; //death of infected cell
			    
#ifdef CYTOKINE_CODE
			  //with cytokine effects
			  }else if(parameters->cyt_effect_infncelldeath>0&&
				dynamics->cytokine_matrix[i][j] > 0 &&
				dice_cell_death<=(((1-cytokine_factor)*
				(parameters->cyt_effect_infncelldeath-
				 inf_cell_death) + 
				inf_cell_death)/parameters->max_rate)){

			    dynamics->cell_state[i][j]=CYTO_KILLED; //death of infected cell		
			    dynamics->cyto_kills++;
#endif
			    
#ifdef TCELL_CODE
			  //Event 2.3: killed by a T-cell
			  }else if(parameters->tcell_killing_target==INF){
			    int num_t_cells=detect_inf(settings,dynamics,3,i,j,T_KILL,SUM,coords);
			    
			    if(num_t_cells>0&&dice<=( parameters->tcell_killing_rate*num_t_cells)/parameters->max_rate){
				dynamics->cell_state[i][j]=TRM_KILLED; 
				dynamics->trm_kills++;
			    }
#endif
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

#ifdef CYTOKINE_CODE
			//Event 3.0: secrete cytokines and virus
			dynamics->delta_cytokine_matrix[i][j]+=
			    MAX(0, gsl_ran_weibull(settings->ur,parameters->cyt_secretion_rate_infcell_mean/parameters->max_rate,5.0));
#endif
			int s = gsl_rng_uniform_int(settings->ur,2);

			unsigned long int new_virions =
			    MAX(0, gsl_ran_weibull(settings->ur,parameters->vir_prod_rate_mean/parameters->max_rate,5.0));
		      //without cytokine effects
#ifdef CYTOKINE_CODE
			if(parameters->cyt_effect_virprod==0 ||
			    dynamics->cytokine_matrix[i][j] == 0) {
#endif
			    if (dynamics->cell_state[i][j] == INFP_A || dynamics->cell_state[i][j] == INFP_ABn)
				dynamics->delta_viral_matrix[0][i][j]+=new_virions;
			    if (dynamics->cell_state[i][j] == INFP_B || dynamics->cell_state[i][j] == INFP_ABn)
				dynamics->delta_viral_matrix[1][i][j]+=new_virions;
			    if (dynamics->cell_state[i][j] == INFP_ABr)
				dynamics->delta_viral_matrix[2][i][j]+=new_virions;
			    
			    dynamics->tot_virons+= new_virions;
			    if (settings->verbose)
			      fprintf(stderr,"Infected cell %d,%d made virus at t=%3.2lf\n",i,j,dynamics->time);
#ifdef CYTOKINE_CODE
			}else{
		      //with cytokine effects
			    cytokine_factor = new_cyt_effect(dynamics->cytokine_matrix[i][j],
			      parameters->prod_ic50);
			    if (dynamics->cell_state[i][j] == INFP_A || (s==0 && dynamics->cell_state[i][j] == INFP_ABn))
				dynamics->delta_viral_matrix[0][i][j]+=new_virions*cytokine_factor;
			    if (dynamics->cell_state[i][j] == INFP_B || (s==1 && dynamics->cell_state[i][j] == INFP_ABn))
				dynamics->delta_viral_matrix[1][i][j]+=new_virions*cytokine_factor;
			    if (dynamics->cell_state[i][j] == INFP_ABr)
				dynamics->delta_viral_matrix[2][i][j]+=new_virions*cytokine_factor;

			    dynamics->tot_virons+= new_virions*cytokine_factor;
			    virions_avoided += new_virions*cytokine_factor;
			    if (settings->verbose)
			      fprintf(stderr,"Infected cell %d,%d made less virus at t=%3.2lf\n",i,j,dynamics->time);
			}
			if (dynamics->cytokine_matrix[i][j] > 0) {
			    cytokine_factor = new_cyt_effect(dynamics->cytokine_matrix[i][j],
			      parameters->inf_ic50);
			} else {
			    cytokine_factor = 1;
			}
#endif
		  
#ifdef TCELL_CODE
			int num_t_cells=detect_inf(settings,dynamics,3,i,j,T_KILL,SUM,coords);
#endif
			
			//Event 3.1: die
#ifdef CYTOKINE_CODE
		      //without cytokine effects
			if ((parameters->cyt_effect_infpcelldeath==0||
			    dynamics->cytokine_matrix[i][j] == 0)&&
			   dice_cell_death<=(inf_cell_death)/parameters->max_rate){
#else
			if (dice_cell_death<=(inf_cell_death)/parameters->max_rate){
#endif
			  dynamics->cell_state[i][j]=DEAD; //death of infected cell
			  if (settings->verbose)
			    fprintf(stderr,"Infected cell %d,%d died at t=%3.2lf\n",i,j,dynamics->time);
			  
#ifdef CYTOKINE_CODE
		      //with cytokine effects
			}else if(parameters->cyt_effect_infpcelldeath>0&&
			    dynamics->cytokine_matrix[i][j] > 0 &&
			    dice_cell_death<=(((1-cytokine_factor)*
			    (parameters->cyt_effect_infpcelldeath-
			     inf_cell_death) + 
			    inf_cell_death)/parameters->max_rate)){

				dynamics->cell_state[i][j]=CYTO_KILLED; //death of infected cell	
				dynamics->cyto_kills++;
			
#endif
#ifdef TCELL_CODE
			//Event 3.2: killed by a T-cell
			}else if(num_t_cells>0&&dice<=(parameters->tcell_killing_rate*num_t_cells)/parameters->max_rate){
			  dynamics->cell_state[i][j]=TRM_KILLED; //death of infected cell
			  dynamics->trm_kills++;
			  if (settings->verbose)
			    fprintf(stderr,"Infected cell %d,%d killed at t=%3.2lf\n",i,j,dynamics->time);
			
#endif
			}
		}
	
#ifdef CYTOKINE_CODE
		////Cytokine layer
		//// 0: no cytokine_matrix, 1: cytokine_matrix present
		////--------------
		//State 1: cytokine_matrix present
		if(dynamics->cytokine_matrix[i][j]>0){
		  //Event 0a: Decay
		  cyt_decay(i,j,parameters,dynamics);

		  if (dynamics->cell_state[i][j] < DEAD &&
		      dynamics->cell_state[i][j] != EMPTY) {
		      cyt_uptake(i,j,parameters,dynamics);
		  }

		  if (dynamics->t_cell_state[i][j] != DEAD &&
		      dynamics->t_cell_state[i][j] != EMPTY) {
		      cyt_uptake(i,j,parameters,dynamics);
		  }
			
		  
		  //Diffusion of cytokine_matrix
		  if(dice_cyt<=parameters->cyt_diff_rate/parameters->max_rate){
		    diffusion_cyt(settings,i,j,gaussian_mask_cyt_diff,parameters,dynamics);}	
		}
#endif

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
#ifdef TCELL_CODE
	    if (settings->save_ratio_data)
		save_et_ratios(settings,dynamics);
#endif
	    next_screen_capture += settings->refr_freq;
	}
	if(settings->stop_early && dynamics->virions==0 && //dynamics->virions<100 && 
	   parameters->neur_drip==0 &&
	   parameters->neur_infect==0 &&
	   dynamics->infected_prodcells==0){
	  break;
	}
    }
    
#ifdef TCELL_CODE
    fprintf(stderr,"In %g days, Added %d HSV+ Tems and %d bystander Tems (%d exited)\n",
	dynamics->time,dynamics->hsv_tems,dynamics->byst_tems,
	dynamics->tem_exits);
#endif

#ifdef CYTOKINE_CODE
    fprintf(stderr,"Cytokines saved %d infections and %lf log virions\n",
	infections_avoided,(virions_avoided > 0)?log10(virions_avoided):0);
#endif
    
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

#ifdef TCELL_CODE
	dynamics->pat_hsv=0;
	dynamics->pat_byst=0;
	dynamics->act_hsv=0;
	dynamics->act_byst=0;
	dynamics->dead_tcells=0;

	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<settings->L;j++)
	{
	    if (dynamics->t_cell_state[i][j]==PAT_HSV)
		dynamics->pat_hsv++;
	    if (dynamics->t_cell_state[i][j]==ACT_HSV)
		dynamics->act_hsv++;
	    if (dynamics->t_cell_state[i][j]==PAT_BYST)
		dynamics->pat_byst++;
	    if (dynamics->t_cell_state[i][j]==ACT_BYST)
		dynamics->act_byst++;
	    if (dynamics->t_cell_state[i][j]==DEAD)
		dynamics->dead_tcells++;
	}
#endif
#ifdef CYTOKINE_CODE
	dynamics->cytokine=0;
	dynamics->cyto_cells=0;

	for (int i=0; i<settings->L;i++)
	    for (int j=0; j<settings->L;j++)
	{
	    dynamics->cytokine+=dynamics->cytokine_matrix[i][j];
	    if (dynamics->cytokine_matrix[i][j] > 0)
		dynamics->cyto_cells++;
	}

#endif
}

// This routine generates state data so that it can be used to make
// plots either while the sim runs (using -p <script> option) or
// at a later time on this or other platforms (see scripts directory)
void dump_data(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics){
    string output_file1;
    string output_file2;
#ifdef TCELL_CODE
    string output_file3;
#endif
#ifdef CYTOKINE_CODE
    string output_file4;
#endif

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

#ifdef TCELL_CODE
    sprintf(temp_str,"tcell_state_%d_%4.3lf.csv",
	dynamics->runnum+1,dynamics->time);

    if (settings->dirname != NULL) {
	output_file3=settings->dirname;
	output_file3+="/";
	output_file3+=temp_str;
    }
    else
	output_file3=temp_str;

    if((outF = fopen(output_file3.c_str(),"wt")) == NULL){
	fprintf(stderr,"Could not open data file %s\n",output_file3.c_str());
	exit(1);
    }
    for (int i=0; i<settings->L-1;i++) {
	    fprintf(outF,"%d,",i+1);
    }
    fprintf(outF,"%d\n",settings->L);
    for (int i=0; i<settings->L;i++) {
	for (int j=0; j<settings->L-1;j++) {
	    fprintf(outF,"%d,",dynamics->t_cell_state[i][j]);
	}
	fprintf(outF,"%d\n",dynamics->t_cell_state[i][settings->L-1]);
    }
    fclose(outF);
#endif

#ifdef CYTOKINE_CODE
    sprintf(temp_str,"cyto_state_%d_%4.3lf.csv",
	dynamics->runnum+1,dynamics->time);

    if (settings->dirname != NULL) {
	output_file4=settings->dirname;
	output_file4+="/";
	output_file4+=temp_str;
    }
    else
	output_file4=temp_str;

    if((outF = fopen(output_file4.c_str(),"wt")) == NULL){
	fprintf(stderr,"Could not open data file %s\n",output_file4.c_str());
	exit(1);
    }
    for (int i=0; i<settings->L-1;i++) {
	    fprintf(outF,"%d,",i+1);
    }
    fprintf(outF,"%d\n",settings->L);
    for (int i=0; i<settings->L;i++) {
	for (int j=0; j<settings->L-1;j++) {
	    fprintf(outF,"%lf,",dynamics->cytokine_matrix[i][j]);
	}
	fprintf(outF,"%lf\n",dynamics->cytokine_matrix[i][settings->L-1]);
    }
    fclose(outF);
#endif

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

#ifdef CYTOKINE_CODE
	double cyt_inf_ic50;
	double cyt_sus_ic50;
	int cyt_inf;
	int cyt_sus;

	if (parameters->cyt_effect_infncelldeath > 0 ||
		parameters->cyt_effect_infpcelldeath > 0) {
	    cyt_inf_ic50=parameters->inf_ic50;
	    cyt_inf=1;
	} else {
	    cyt_inf_ic50=1000;
	    cyt_inf=0;
	}
	if (parameters->cyt_effect_infectivity) {
	    cyt_sus_ic50=parameters->beta_ic50;
	    cyt_sus = 1;
	} else if (parameters->cyt_effect_virprod) {
	    cyt_sus_ic50=parameters->prod_ic50;
	    cyt_sus = 1;
	} else if (settings->cyto_act) {
	    cyt_sus_ic50=parameters->cyt_trm_ic50;
	    cyt_sus = 1;
	} else {
	    cyt_sus_ic50=1000;
	    cyt_sus=0;
	}


	sprintf(temp_str,
	    "%d %4.3lf %d %4.3lf %e %e %4.3lf %4.3lf %d %d %d %d",
	    dynamics->runnum+1,dynamics->time, settings->plot_cytokines,
	    parameters->cyt_prot_expiration, cyt_inf_ic50, cyt_sus_ic50,
	    (dynamics->max_vl>0)?log10(dynamics->max_vl):0, 
	    dynamics->max_vl_time,settings->plot_trms,settings->cyto_act,
	    cyt_inf, cyt_sus);
#else
	sprintf(temp_str, "%d %4.3lf %4.3lf %4.3lf",
	    dynamics->runnum+1,dynamics->time, (dynamics->max_vl>0)?log10(dynamics->max_vl):0, dynamics->max_vl_time);
#endif

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
#ifdef TCELL_CODE
		remove(output_file3.c_str());
#endif
#ifdef CYTOKINE_CODE
		remove(output_file4.c_str());
#endif
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
#ifdef TCELL_CODE
	    fprintf(settings->dataF1,",pat hsv Trm");
	    fprintf(settings->dataF1,",act hsv Trm");
	    fprintf(settings->dataF1,",pat byst Trm");
	    fprintf(settings->dataF1,",act byst Trm");
	    fprintf(settings->dataF1,",gauss offset");
	    fprintf(settings->dataF1,",avg E/T ratio");
	    fprintf(settings->dataF1,",3x3 E/T ratio");
	    fprintf(settings->dataF1,",5x5 E/T ratio");
	    fprintf(settings->dataF1,",7x7 E/T ratio");
	    fprintf(settings->dataF1,",9x9 E/T ratio");
	    fprintf(settings->dataF1,",11x11 E/T ratio");
	    fprintf(settings->dataF1,",percent 0s");
	    fprintf(settings->dataF1,",closest HSV+");
	    fprintf(settings->dataF1,",closest Tcell");
	    fprintf(settings->dataF1,",Tcell kills");
	    fprintf(settings->dataF1,",first trm antigen");
	    fprintf(settings->dataF1,",first trm cycto");
	    fprintf(settings->dataF1,",hsv tems");
	    fprintf(settings->dataF1,",byst tems");
	    fprintf(settings->dataF1,",tem exits");
#endif
#ifdef CYTOKINE_CODE
	    fprintf(settings->dataF1,",cytokines");
	    fprintf(settings->dataF1,",cyto_cells");
	    fprintf(settings->dataF1,",cytokine kills");
#endif
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
#ifdef TCELL_CODE
	    fprintf(settings->dataF1,",%d", dynamics->pat_hsv);
	    fprintf(settings->dataF1,",%d", dynamics->act_hsv);
	    fprintf(settings->dataF1,",%d", dynamics->pat_byst);
	    fprintf(settings->dataF1,",%d", dynamics->act_byst);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->gauss_offset);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->avg_et_ratio);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->et_ratio_3);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->et_ratio_5);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->et_ratio_7);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->et_ratio_9);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->et_ratio_11);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->perc_zeros);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->closestHSV);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->closestTcell);
	    fprintf(settings->dataF1,",%d", dynamics->trm_kills);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->trm_act_time);
	    fprintf(settings->dataF1,",%4.3lf", dynamics->trm_cyt_time);
	    fprintf(settings->dataF1,",%d", dynamics->hsv_tems);
	    fprintf(settings->dataF1,",%d", dynamics->byst_tems);
	    fprintf(settings->dataF1,",%d", dynamics->tem_exits);
#endif
#ifdef CYTOKINE_CODE
	    fprintf(settings->dataF1,",%g", dynamics->cytokine);
	    fprintf(settings->dataF1,",%d", dynamics->cyto_cells);
	    fprintf(settings->dataF1,",%d", dynamics->cyto_kills);
#endif
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
