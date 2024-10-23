

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

void python_plot_BM( std::vector<double> center, double inner_radius, double outter_radius, std::vector<int>& bm_voxels);

#endif
