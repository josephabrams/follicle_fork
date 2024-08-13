#ifndef __CUSTOM_COLORING_H__
#define __CUSTOM_COLORING_H__
#include "../core/PhysiCell.h"

/*#include "../modules/PhysiCell_SVG.h"*/
/*#include "../BioFVM/BioFVM_utilities.h"*/
#include "../modules/PhysiCell_pathology.h"
#include <vector>
#include <string>

using namespace PhysiCell;
using namespace BioFVM;


std::string formatted_seconds_to_HHMMSS( double seconds );

void Custom_SVG_plot( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::vector<std::string> (*substrate_coloring_function)(double, double, double) = NULL, void (*cell_counts_function) (char*) = NULL); // done

#endif
