#include "./spring_class.h"
Spring::Spring( Cell* neighbor, double rest_length, double spring_constant, bool is_TZP )
    : m_neighbor{neighbor}, m_rest_length{rest_length}, m_spring_constant{spring_constant}, m_is_TZP{is_TZP} {
  m_force.resize(3,0.0);
}

Spring_Connections::Spring_Connections()
{
  this->neighbor_springs.clear();
  this->previous_net_force.resize(3,0.0);
  this->current_net_force.resize(3,0.0);
  this->next_velocity.resize(3,0.0);
  this->current_velocity.resize(3,0.0);
  this->previous_mass=0;
  this->previous_mass=0;
}
void Spring_Connections::sync_cell(Cell* pCell)
{
 //done in cryocell.cpp 
  return;
}
void Spring_Connections::spring_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt)
{
  std::vector<double> previous_acceleration=(1/this->previous_mass)*this->previous_net_force;
  std::vector<double> current_acceleration=(1/this->previous_mass)*this->previous_net_force;

  Adams_Bashforth_2_vec(&this->next_velocity, this->current_velocity, current_acceleration, previous_acceleration,dt);
  return;
}

void Spring_Connections::spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt)
{
  //double check how this works with pCell->velocity, might be better to just handle everything
  //if not connected neighbor (sync_connected_neighbors)
  //adhesive and repulsive with 0 spring length

  return;
}
void Spring_Connections::calculate_spring_force(Spring& spring)
{ 
  //if connected neighbor
  double spring_length=norm(spring.m_neighbor->position-this->m_pCell->position)-spring.m_neighbor->phenotype.geometry.radius-this->m_pCell->phenotype.geometry.radius;
  double delta_x=spring_length-spring.m_rest_length;
  std::vector<double> unit_vec=1/norm(spring.m_neighbor->position-this->m_pCell->position)*(spring.m_neighbor->position-this->m_pCell->position); 
  std::vector<double> force=(spring.m_spring_constant*(delta_x))*unit_vec;
  spring.m_force=force;
  //else displacement==0

  //figure out signs
  return;
}

void Spring_Connections::update_springs(std::vector<Spring*> &springs){
 return; 
}
void Spring_Connections::add_spring(Spring* spring)
{
  #pragma omp critical
  {
    this->neighbor_springs.push_back(spring);
  }
  return;
}

void Spring_Connections::remove_spring(Spring* spring)
{
  std::vector<Spring*>::iterator it;
  int vec_position=0;
  it = find (this->neighbor_springs.begin(), this->neighbor_springs.end(), spring);
  if (it != neighbor_springs.end())
  {
    vec_position=it-this->neighbor_springs.begin();
	  this->neighbor_springs[vec_position]= this->neighbor_springs[neighbor_springs.size()-1];//swap spring locations and pop off the list
	  this->neighbor_springs[neighbor_springs.size()-1] = spring;
	  this->neighbor_springs.pop_back();	
  }
  else{
    std::cout << "TRIED TO REMOVE UNKNOWN SPRING\n";
  }
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
