

#ifndef __CRYOCELL_h__
#define __CRYOCELL_h__

#include "../../core/PhysiCell.h"
#include <vector>
#include "../multivoxel/multivoxel_functions.h"
#include "../cryomodule/ABFM.h"
#include <string>
#include <omp.h>
using namespace PhysiCell;
using namespace BioFVM;

//for future use
/*
class Youngs_Modulus{
private:
public:
  Cell* neighbor;
  double rest_length;
  Matrix stress;
  Matrix strain;
  Vector normal;
  double cross_sectional_area;
  void force_function();
};
*/
class Spring{
private:
public:
  Cell* m_neighbor;
  double m_rest_length;
  double m_spring_constant;
  Spring( Cell* neighbor, double rest_length, double spring_constant );
};

class Spring_Connections{
private:
public:
  Cell* m_pCell;
  std::vector<Spring*> neighbor_springs;
  std::vector<double> previous_net_force;
  std::vector<double> current_net_force;
  std::vector<double> previous_velocity;
  std::vector<double> current_velocity;
  double previous_mass;//volume is a proxy for mass
  double current_mass;//if I shrink my mass goes down so force goes up from retreating membrane and acceleration goes up from reduced mass

  void spring_function();
  
  Spring_Connections(Cell* pCell);
	//format should be void (*update_velocity)( Cell* pCell, Phenotype& phenotype, double dt ); 
  void spring_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt);//calculate velocity using ABM 2nd order of force
  void spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt); //zero rest length
  void calculate_spring_force(Spring& spring);
  void update_springs(std::vector<Spring*> &springs);

  void add_spring(Cell* pCell);
  void remove_spring(Cell* pCell);

};

#endif //__CRYOCELL_h__
