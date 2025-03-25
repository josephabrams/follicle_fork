#ifndef __SPRING_H__
#define __SPRING_H__

#include "../../core/PhysiCell.h"
#include <vector>
#include "../multivoxel/multivoxel_functions.h"
#include "../cryomodule/ABFM.h"
#include <string>
#include <omp.h>
#include <cstddef>
#include <algorithm>
using namespace PhysiCell;
using namespace BioFVM;

//prototype for future version
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
  Cell* m_me;
  Cell* m_neighbor;
  double m_rest_length;
  double m_spring_constant;
  bool m_is_TZP;
  bool m_is_broken;
  std::vector<double> m_force;
  std::vector<double> m_previous_force;
  Spring( Cell* me, Cell* neighbor, double rest_length, double spring_constant, bool is_TZP );
  ~Spring();
  void calculate_spring_force();
  
  void hookes_law_membrane_pressure(double spring_length, double membrane_pressure);
  void hookes_law(double spring_length);
  void calculate_youngs_modulus(double spring_length);
  void update_force_vector(std::vector<double> *my_return_force, std::vector<double> *neighbor_return_force);
  void update_spring_velocity();
  void remove_spring();
  void test_TZPs();
};

int TZP_count();
// class Spring_Connections{
// private:
// public:
//   Cell* m_pCell;
//   std::vector<Spring*> neighbor_springs;
//   std::vector<double> previous_net_force;
//   std::vector<double> current_net_force;
//   std::vector<double> next_velocity;
//   std::vector<double> current_velocity;
//   double previous_mass;//volume is a proxy for mass
//   double current_mass;//if I shrink my mass goes down so force goes up from retreating membrane and acceleration goes up from reduced mass
//
//   void spring_function();
//
//   Spring_Connections();
//
//   void sync_cell(Cell* pCell);
// 	//format should be void (*update_velocity)( Cell* pCell, Phenotype& phenotype, double dt ); 
//   void spring_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt);//calculate velocity using ABM 2nd order of force
//   void spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt); //zero rest length
//
//   void calculate_spring_force(Spring& spring);
//   void update_springs(std::vector<Spring*> &springs);
//
//   void add_spring(Spring* spring);
//   void remove_spring(Spring* spring);
//
// };

Spring* create_spring( Cell* me, Cell* neighbor, double rest_length, double spring_constant, bool is_TZP );
void delete_spring();
void calculate_all_spring_forces();
void calculate_spring_velocity();
void TZPs();


class Point_Spring{
private:
public:
  Cell* m_me;
  std::vector<double> m_force_normal;
  double m_rest_length;
  double m_spring_constant;
  double m_spring_length;
  std::vector<double> m_force;
  std::vector<double> m_previous_force;
  
  Point_Spring( Cell* me, double rest_length, double spring_constant);
  ~Point_Spring();

  void calculate_spring_force();
  void hookes_law_membrane_pressure(double spring_length, double membrane_pressure);
  void hookes_law(double spring_length);
  // void calculate_youngs_modulus(double spring_length);
  void update_force_vector(std::vector<double> *my_return_force);
  void update_spring_velocity();

};

Point_Spring* create_point_spring( Cell* me, double rest_length, double spring_constant );

void delete_point_spring();
void calculate_all_point_spring_forces();
void calculate_point_spring_velocity();

Spring* find_spring( Cell* me, Cell* neighbor);


void output_TZP_csv(double k_oocyte, double k_granulosa, double k_basement, std::string sim_num);

void create_output_TZP_csv(std::string sim_num);

void outter_constraint(double outter_bound);

void inner_constraint(double inner_bound);


extern std::vector<Spring*> all_springs;
extern std::vector<Point_Spring*> all_point_springs;
//external variables for passing parameters into main() when running script on HPC
extern double GRANULOSA_K;
extern double OOCYTE_K;
extern int TZP_COUNT;
extern int INIT_TZP_COUNT;
#endif //__SPRING_H__
