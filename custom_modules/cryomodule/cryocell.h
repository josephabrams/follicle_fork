
#ifndef __CRYOCELL_h__
#define __CRYOCELL_h__

#include "../../core/PhysiCell.h"
#include <vector>
#include "../multivoxel/multivoxel_functions.h"
#include "./conversions.h"/*class Cryocell_State;*/
#include "./volume_change.h"
#include "ABFM.h"
#include <string>
#include <omp.h>
using namespace PhysiCell;
using namespace BioFVM;
/*class Cryo_Concentrations;*/
/*class Cryo_Parameters;*/
class Cryocell;
class Cryo_Parameters
{
  private:
  public:
    double osmotically_inactive_fraction;
    // water flux
    double Lp;
    double dVw;//cell water flux
    double previous_dVw;
    // solute flux
    std::vector<double> Ps;
    std::vector<double> dN;//cell mole flux of solutes
    std::vector<double> previous_dN;//previous mole flux of solutes

    Cryo_Parameters(Cryocell* cCell);
};
class Cryo_Concentrations
{
  private:
  public:
    
    bool use_virial;
    double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
    double interior_osmolality;// total internal osmolality salt+CPA
    std::vector<double> interior_molarity;
    std::vector<double> exterior_molarity;
    std::vector<double> interior_component_molality;
    std::vector<double> exterior_component_molality;

    Cryo_Concentrations(Cryocell* cCell);
};
class Cryocell_State 
{
  private:
  public:
    double previous_radius;
    double surface_area;
    double temperature;
    
    double solute_volume;
    
    std::vector<double> solute_moles;
    std::vector<double> next_solute_moles;
   
    double next_water_volume;//for ABM 2nd order
    double water_volume;
    
    //voxel uptakes
    double water_uptake;//um^3
    std::vector <double> solute_uptake;
    std::vector <double> solute_uptake_per_voxel;
    double water_uptake_per_voxel;
    std::vector<double> uptake;//molar uptake/secretion of solutes
    std::vector <int> uptake_voxels; //voxels changing from uptake

    Cryocell_State(Cryocell* cCell);
    
};

class Cryocell : public PhysiCell::Cell {

  private:
  public:
    Cryo_Concentrations cryo_concentrations{this};
    Cryo_Parameters cryo_parameters{this};
    Cryocell_State cryocell_state{this};
    double solid_volume;
    double toxicity;
    std::vector <int> cell_voxels;
    std::vector <int> neighbor_voxels;
    Cryocell();
    ~Cryocell(){};
    void update_cell_voxels();
    void update_neighbor_voxels();
  //note it might make sense to just make a threadsafe_write function for return values

};

Cell* instantiate_Cryocell();

Cell* create_Cryocell(Cell_Definition& cd);
void update_all_cells_voxels();
void get_concentration_at_boundary();//currently only 2D, 3D should probably use boundary voxels
void get_exterior_molalities();

void get_exterior_osmolalities();

void get_interior_molalities();

void get_interior_osmolalities();

void update_exterior_concentrations();
void update_interior_concentrations();
void calculate_derivatives();
void advance_osmosis(double dt);
void calculate_uptakes(double dt);
void update_next_step(double dt);
void calculate_per_voxel_uptake();

void uptake_in_one_voxel(int &voxel, double& water_uptake_per_voxel, std::vector<double>& solute_uptake_per_voxel);

void advance_uptake();
void uptake(double dt);
void two_p_forward_step(double dt);

void two_p_update_volume();
/*void multistep_loading(double dt); */

/*bool isMultivoxel(PhysiCell::Cell* pCell);*/

/*bool isMultivoxel(PhysiCell::Cell_Definition * cellDef);*/

/*std::vector<PhysiCell::Cell_Definition*>* getMultivoxelCellDefinitions();*/
extern std::vector<Cryocell*> all_cryocells;
#endif
