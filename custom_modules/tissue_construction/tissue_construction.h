
#ifndef __TISSUE_CONSTRUCTION_H__
#define __TISSUE_CONSTRUCTION_H__
#include <array>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include "./../../core/PhysiCell.h"
#include "./../../core/PhysiCell_cell.h"
#include "./../../core/PhysiCell_cell_container.h"
#include "./../../modules/PhysiCell_standard_modules.h"
#include <string>
#include <vector>
using namespace BioFVM;
using namespace PhysiCell;

std::vector<std::vector<double>> create_spheroid_2D(double cell_radius, double sphere_radius); 
std::vector<std::vector<double>> create_spheroid(double cell_radius, double sphere_radius); 

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius,double inner_radius); 

/*void setup_single_cell(void);*/
/*void setup_4_cell_test_case(void);*/
/*void setup_secondary_stage_follicle(void);*/
/*void setup_toy_granulosa_model(void);*/





#endif //__TISSUE_CONSTRUCTION_H__
