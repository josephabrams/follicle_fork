#ifndef __MULTIVOXEL_UPTAKE_H__
#define __MULTIVOXEL_UPTAKE_H__
#include <array>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>
#include "../../core/PhysiCell_cell.h"
using namespace BioFVM;
using namespace PhysiCell;


void uptake_in_one_voxel(int &voxel, double& water_uptake_per_voxel, std::vector<double>& solute_uptake_per_voxel);

#endif
