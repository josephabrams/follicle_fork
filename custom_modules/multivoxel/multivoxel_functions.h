
#ifndef __MULTIVOXEL_FUNCTIONS_H__
#define __MULTIVOXEL_FUNCTIONS_H__
#include <array>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>
#include "../../core/PhysiCell_cell.h"
using namespace PhysiCell;
using namespace BioFVM;
// void general_voxel_bounding_box(std::vector<int> *return_bounding_box,std::vector<double> &starting_position, std::vector <double>&ending_position, double &voxel_length, BioFVM::Cartesian_Mesh &a_mesh);

void general_voxel_bounding_box(std::vector<int> *return_bounding_box, std::vector<double> center, std::vector<double> half_dimensions, double voxel_length, BioFVM::Cartesian_Mesh &a_mesh);


// void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index);

void get_intersecting_voxels(Cell* pCell,std::vector<int>& bounding_voxels,std::vector<int>* return_intersecting_voxel_indicies);
void get_voxel_corners(std::vector<double> &voxel_center, std::vector<std::vector<double>> *return_corners );
// void get_intersecting_voxels(Cell* pCell,std::vector<int>* return_intersecting_voxel_indicies);
void get_exterior_voxels(Cell *pCell, std::vector<int>* return_exterior_voxel_indicies);
void get_interior_voxels(Cell *pCell, std::vector<int>* return_interior_voxel_indicies);


void intersecting_neighbor_voxels(Cell* pCell, Cell* pNeighbor, std::vector<int> my_bounding_voxels, std::vector<int> neighbor_bounding_voxels, std::vector<int> *return_voxels);
// void cells_in_me(Cell *pCell, std::vector<Cell*> *return_cells_in_me); // uses mechanics vectors to search for cells that are within or equal to pCell radius can also be used for bounding boxes
void python_plot_cell_and_voxels(Cell* pCell, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name);

void python_plot_two_cells_and_voxels(Cell* pCell, Cell* pNeighbor, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name);
#endif
