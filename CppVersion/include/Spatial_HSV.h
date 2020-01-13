#ifndef SPATIAL_HSV_H
#define SPATIAL_HSV_H
#include "settings.h"

void set_default_settings(global_settings *settings);
void set_default_parameters(global_parameters *parameters);
void alloc_grid(global_settings *settings,global_dynamics *dynamics);
void init_grid(int runnum, global_settings *settings,global_parameters *parameters,global_dynamics *dynamics);
void simulation(int runnum, global_settings *settings,global_parameters *parameters,global_dynamics *dynamics);
gsl_matrix *compute_gaussian(double mask_range,double sig,double amp);
void place_host(global_parameters *parameters,global_dynamics *dynamics);
void place_virus(global_settings *settings,global_dynamics *dynamics);
void place_tcells(global_settings *settings, global_parameters *parameters,global_dynamics *dynamics);
void neur_drip(int i,int j,global_dynamics *dynamics);
double calc_plaque_size(global_settings *settings);
void diffusion_cyt(global_settings *settings,int i, int j,gsl_matrix *gaussian_mask,global_parameters *parameters, global_dynamics *dynamics);
void diffusion_virus(global_settings *settings,int i, int j,gsl_matrix * gaussian_mask,global_parameters *parameters);
void gaussburst(double burst_size,gsl_matrix *gaussian_mask);
bool close_to_plaq(int i,int j, global_settings *settings,global_parameters *parameters, global_dynamics *dynamics);
void distrib_progeny(global_settings *settings,global_parameters *parameters,
	global_dynamics *dynamics, int what, int i, int j,double burst,
	gsl_matrix *gaussian_mask);
bool pick_nbr(global_settings *settings,int i,int j, int *coords);
int detect_inf(global_settings *settings, global_parameters *parameters, double dice_cells, 
	global_dynamics *dynamics, int det_range,int i, int j,int what,int return_pos,int *coords);
void set_trm_move_function(int trm_move_method);
void tcell_div(int i, int j,global_settings *settings, global_parameters *parameters, global_dynamics *dynamics);
void viral_decay(int i, int j,double dice,global_parameters *parameters, global_dynamics *dynamics);
void cyt_decay(int i, int j,global_parameters *parameters, global_dynamics *dynamics);
void compute_totals(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics);
void clear_counts(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics);
void clear_deltas(global_settings *settings,global_dynamics *dynamics);
void propagate_deltas(global_settings *settings,global_dynamics *dynamics);
double cyt_effect(double dose, global_parameters *parameters);
double new_cyt_effect(double dose, global_parameters *parameters);
void move_tcell(int i, int j,double dice_cells,global_settings *settings, global_parameters *parameters, global_dynamics *dynamics);
void output_results(int runnum, global_settings *settings,
	global_parameters *parameters,global_dynamics *dynamics,bool header);
void dump_data(global_settings *settings,global_parameters *parameters,global_dynamics *dynamics);
void end_simulation(global_dynamics *dynamics);
//void output_results(global_settings *settings,results,header);
#endif
