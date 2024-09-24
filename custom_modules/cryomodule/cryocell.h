
#ifndef __CRYOCELL_h__
#define __CRYOCELL_h__

#include "../../core/PhysiCell.h"
#include <vector>
#include "../multivoxel/multivoxel_functions.h"

#include "../multivoxel/multivoxel_neighborhood.h"
#include "./conversions.h"/*class Cryocell_State;*/
#include "./volume_change.h"
#include "ABFM.h"
#include <string>
#include <omp.h>

#include "../spring_class/spring_class.h"
using namespace PhysiCell;
using namespace BioFVM;
/*class Cryo_Concentrations;*/
/*class Cryo_Parameters;*/
// class Spring_Connections;
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

    Cryo_Parameters();
    void sync_to_cell_definition(Cell_Definition& cd); 
};
class Cryo_Concentrations
{
  private:
  public:
    
    bool use_virial; //not currently doing anything
    double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
    double interior_osmolality;// total internal osmolality salt+CPA
    std::vector<double> interior_molarity;
    std::vector<double> exterior_molarity;
    std::vector<double> interior_component_molality;
    std::vector<double> exterior_component_molality;

    Cryo_Concentrations();
    // void sync_to_cell_definition(Cell_Definition& cd); not needed atm 
};
class Cryocell_State 
{
  private:
  public:
    bool is_cryocell; //incase someone casts a cell they shouldnt
    double previous_radius;
    double surface_area;
    double temperature;
    
    double solid_volume;
    double toxicity;
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

    Cryocell_State();
    void sync_to_cell_definition(Cell_Definition& cd, Cryo_Parameters& cryo_p); 
    void sync_moles_and_volume(Cell_Definition& cd, Cryo_Parameters& cryo_p, Cryo_Concentrations& cc);
};

class Cryocell : public PhysiCell::Cell {

  private:
  public:
    Cryo_Concentrations cryo_concentrations;
    Cryo_Parameters cryo_parameters;
    Cryocell_State cryocell_state;
    // Spring_Connections spring_connections;
    std::vector <int> cell_voxels;
    std::vector <int> neighbor_voxels;
    std::vector<Cell*> initial_neighbors;
    std::vector<Cell*> all_neighbors;
    std::vector<double> net_force;
    std::vector<double> previous_net_force;
    double mass;
    std::vector<double> old_position;
    Cryocell();
    ~Cryocell(){};
    void update_cell_voxels();
    void update_neighbor_voxels();
    void sync_spring_connections();
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

void update_multivoxel_neighboorhood();

void update_initial_neighbors();

void update_springs();

void update_velocity();
void update_all_spring_forces();
void sum_spring_forces(Cryocell* cCell);

void sum_youngs_modulus(Cryocell* cCell);
void cell_to_cell_youngs_modulus( Cryocell* pMe, Cell* pOther, std::vector<double> *return_force);

void update_net_force();
void calculate_position_from_acceleration(std::vector<double> &old_position, std::vector<double>&current_position, std::vector<double> &net_acceleration, double dt, std::vector<double> *new_position);
void validate_cell_position(Cryocell* cCell);
void update_position_from_net_force(double dt);
void two_p_update_volume();
/*void multistep_loading(double dt); */

/*bool isMultivoxel(PhysiCell::Cell* pCell);*/

/*bool isMultivoxel(PhysiCell::Cell_Definition * cellDef);*/

/*std::vector<PhysiCell::Cell_Definition*>* getMultivoxelCellDefinitions();*/
extern std::vector<Cryocell*> all_cryocells;
#endif
