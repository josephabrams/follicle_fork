

#ifndef __BASEMENT_MEMBRANE_H__
#define __BASEMENT_MEMBRANE_H__
#include <array>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include "../../core/PhysiCell_cell.h"
#include "../multivoxel/multivoxel_functions.h"
#include "../spring_class/spring_class.h"
#include "../cryomodule/cryocell.h"
using namespace PhysiCell;
using namespace BioFVM;

void spherical_bounding_region(std::vector <double> &center_point, double &radius,std::vector <int> *return_bounding_region);


void get_basement_membrane_voxels(std::vector <double> &center_point, double &inner_radius, double &outter_radius, std::vector<int> *basement_voxels );

void python_plot_BM(double inner_radius, double outter_radius);
void python_plot_BM( std::vector<double> center, double inner_radius, double outter_radius, std::vector<int>& bm_voxels);

void python_plot_BM_cells(double inner_radius, double outter_radius);
void find_multivoxel_neighbors_region(std::vector<int> &voxel_region, std::vector<Cell*> *return_neighbors, double outter_radius, double inner_radius);


void get_initial_BM_neighbors(double outter_radius, double inner_radius);


void create_BM_springs(double outter_radius, double inner_radius);

void update_BM_neighbors(double outter_radius, double inner_radius);

void advance_BM_springs(double outter_radius, double inner_radius);
void linker_test();

extern std::vector<int> basement_membrane_voxels;
extern std::vector<Cell*> basement_neighbors;
extern std::vector<Cell*> basement_initial_neighbors;
extern std::vector<Cell*> outter_neighbors;

extern double BASEMENT_K;
#endif
