

#ifndef __MULTIVOXEL_NEIGHBORHOOD_H__
#define __MULTIVOXEL_NEIGHBORHOOD_H__
#include <array>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>
#include "../../core/PhysiCell_cell.h"
#include "./multivoxel_functions.h"
using namespace PhysiCell;
using namespace BioFVM;
// void general_voxel_bounding_box(std::vector<int> *return_bounding_box,std::vector<double> &starting_position, std::vector <double>&ending_position, double &voxel_length, BioFVM::Cartesian_Mesh &a_mesh);

// void general_voxel_bounding_box(std::vector<int> *return_bounding_box, std::vector<double> center, std::vector<double> half_dimensions, double voxel_length, BioFVM::Cartesian_Mesh &a_mesh);

void find_multivoxel_neighbors(Cell* pCell, std::vector<Cell*> *return_neighbors);//currently Cryocell is the only multivoxel cell-- search a multivoxel cell and all its voxels for neighboors eventually replace with better method

void find_multivoxel_neighbors_direct_contact(Cell* pCell, std::vector<Cell*> *return_neighbors);
// void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index);
// void intersecting_neighbor_voxels(Cell* pCell, Cell* pNeighbor, std::vector<int> *my_bounding_voxels, std::vector<int> *neighbor_bounding_voxels, std::vector<int> *return_voxels);
// void cells_in_me(Cell *pCell, std::vector<Cell*> *return_cells_in_me); // uses mechanics vectors to search for cells that are within or equal to pCell radius can also be used for bounding boxes
void python_plot_cell_with_Neighbors(Cell* pCell, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name, double z_height, std::vector<Cell*> &neighbors);
#endif
