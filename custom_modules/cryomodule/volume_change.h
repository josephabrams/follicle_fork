#ifndef __VOLUME_CHANGE_H__
#define __VOLUME_CHANGE_H__
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <unordered_map>
#include <vector>
#include "./ABFM.h"

#include "../../core/PhysiCell.h"
#include "../../modules/PhysiCell_standard_modules.h"
using namespace BioFVM;
using namespace PhysiCell;


void sum_solute_volume(std::vector<double> &solute_specific_volume, std::vector<double> &solute_moles, std::vector<double> *solute_volume );

void sum_total_volume(std::vector<double> &solute_volume, double water_volume, double solid_volume, double* new_volume);

void dVw_Osmolality( double &Lp,  double &surface_area,  const double &gas_constant,  double &temperature,  double &exterior_osmolality,  double &interior_osmolality, double *previous_dVw, double *dVw);

void dN_molarity( std::vector<double> &Ps,  double &surface_area,  std::vector<double> &exterior_molarity,  std::vector<double> &interior_molarity, std::vector<double> *previous_dN, std::vector<double> *dN);

void calculate_solute_moles(std::vector<double>* next_solute_moles,  std::vector<double>& solute_moles,  std::vector<double>& dN,  std::vector<double> previous_dN,  double dt);

/*void calculate_water_volume(double* next_water_volume,  double &water_volume,  double &dVw,  double &previous_dVw,  double dt);*/

/*void calculate_solute_uptake(std::vector<double>* solute_uptake,  std::vector<double>& next_solute_moles,  std::vector<double>& solute_moles);*/

void calculate_water_uptake( double *water_uptake,  double& next_water_volume,  double& water_volume);

/*void two_p_forward_step( double dt);*/

/*void two_p_update_volume();*/

#endif
