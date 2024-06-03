#include "./volume_change.h"


void get_external_concentration(Cell *pCell, std::vector <double> *average_concentrations, std::vector<int> &solute_indecies, std::vector<int> &voxels );
void update_external_concentration(Cell* pCell);
void update_internal_concentration(Cell* pCell);
void two_p_model(Cell* pCell);
