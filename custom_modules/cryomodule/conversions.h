
#ifndef __CONVERSIONS_h__
#define __CONVERSIONS_h__

#include "../../core/PhysiCell.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>


double molarity_to_molality(double molarity, std::string component_name);

double binary_virial(double molality, std::string component_name);

double ternary_virial(double molality_1, double molality_2, std::string component_1, std::string component_2);

double moles_to_volume(double moles, std::string component_name);
#endif
 
