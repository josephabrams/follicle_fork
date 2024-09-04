#include "./spring_class.h"
#include <cstddef>
Spring::Spring( Cell* neighbor, double rest_length, double spring_constant )
    : m_neighbor{neighbor}, m_rest_length{rest_length}, m_spring_constant{spring_constant} {}

Spring_Connections::Spring_Connections(Cell* pCell)
  : m_pCell{pCell}
{
  this->neighbor_springs.clear();
  this->previous_net_force.resize(3,0.0);
  this->current_net_force.resize(3,0.0);
  this->previous_velocity.resize(3,0.0);
  this->current_velocity.resize(3,0.0);
  this->previous_mass=0;
  this->previous_mass=0;
}

void Spring_Connections::spring_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt)
{
  return;
}

void Spring_Connections::spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt)
{
  //if not connected neighbor (sync_connected_neighbors)
  //adhesive and repulsive with 0 spring length

  return;
}
void Spring_Connections::calculate_spring_force(Spring& spring)
{
  //if connected neighbor
  double displacement=norm(spring.m_neighbor->position-this->m_pCell->position)-spring.m_neighbor->phenotype.geometry.radius-this->m_pCell->phenotype.geometry.radius;
  //else displacement==0
  double delta_x=spring.m_rest_length-displacement;
  //figure out signs
  return;
}


// class Spring_Connections{
// private:
// public:
//   Cell* m_pCell
//   std::vector<Spring*> neighbor_springs;
//   std::vector<double> previous_net_force;
//   std::vector<double> current_net_force;
//   std::vector<double> previous_velocity;
//   std::vector<double> current_velocity;
//   double previous_mass;//volume is a proxy for mass
//   double current_mass;//if I shrink my mass goes down so force goes up from retreating membrane and acceleration goes up from reduced mass
//
//   void spring_function();
//
//   Spring_Connections();
// 	//format should be void (*update_velocity)( Cell* pCell, Phenotype& phenotype, double dt ); 
//   void spring_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt);//calculate velocity using ABM 2nd order of force
//   void spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt); //zero rest length
//   void calculate_spring_force(Spring* spring);
//   void update_springs(std::vector<Spring*> &springs);
// };
