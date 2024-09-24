#include "./spring_class.h"

std::vector<Spring*> all_springs;
Spring::Spring(Cell* me, Cell* neighbor, double rest_length, double spring_constant, bool is_TZP )
    : m_me{me}, m_neighbor{neighbor}, m_rest_length{rest_length}, m_spring_constant{spring_constant}, m_is_TZP{is_TZP} {
  m_force.resize(3,0.0);
  m_previous_force.resize(3,0.0);
  m_is_broken=false;
}
void Spring::calculate_spring_force()
{
    // hooks law should be thread safe
    double spring_length=norm(m_neighbor->position-m_me->position)-m_neighbor->phenotype.geometry.radius-m_me->phenotype.geometry.radius;
    double delta_x=spring_length-m_rest_length;
    std::vector<double> unit_vec=(1/norm(m_neighbor->position-m_me->position))*(m_neighbor->position-m_me->position); 
    std::vector<double> force=(m_spring_constant*(delta_x))*unit_vec;
    m_force=force;

  // std::vector<double> force(3,0.0);
  //if connected neighbor
  // std::cout<<"Spring length: "<< spring_length<<"\n";
  // std::cout<<"delta_x: "<< delta_x<<"\n";
  // std::cout<<"unit_vec: "<< unit_vec<< "\n";
  // std::cout<<"spring constant: "<< spring.m_spring_constant<<"\n";
  // std::cout<<"force magnitude: "<< spring.m_spring_constant*delta_x<<"\n";
  // std::cout<<"FORCE: "<< spring.m_force<<"\n";
  // else displacement==0
}
void Spring::update_force_vector(std::vector<double> *my_return_force, std::vector<double> *neighbor_return_force)
{
  if(std::fabs(norm(m_force))<1e-12)
  {
    m_force={0.0, 0.0, 0.0};
  }
  axpy(my_return_force, 1.0, m_force);
  axpy(neighbor_return_force,-1.0,m_force);
}
void Spring::update_spring_velocity()//for when using only springs!!!
{
  //calculate acceleration for me and neighbor using volume as a proxy for mass
  std::vector<double> me_acceleration= (1/m_me->custom_data["initial_volume"])*m_force;
  std::vector<double> neighbor_acceleration= (-1/m_neighbor->custom_data["initial_volume"])*m_force;
  // get velocity and add it
  // first time step apply forward euler
  if(PhysiCell_globals.current_time<mechanics_dt)
  {
    axpy(&m_me->velocity,mechanics_dt, me_acceleration);
    axpy(&m_neighbor->velocity,mechanics_dt, neighbor_acceleration);
  }
  else {
    std::vector<double> prev_me_acceleration= (1/m_me->custom_data["initial_volume"])*m_previous_force;
    std::vector<double> prev_neighbor_acceleration= (-1/m_neighbor->custom_data["initial_volume"])*m_previous_force;
    std::vector<double> prev_me_velocity=m_me->get_previous_velocity();
    std::vector<double> prev_neighbor_velocity=m_neighbor->get_previous_velocity();
    Adams_Bashforth_2_vec(&m_me->velocity, prev_me_velocity, me_acceleration, prev_me_acceleration, mechanics_dt); 
    Adams_Bashforth_2_vec(&m_neighbor->velocity, prev_neighbor_velocity, neighbor_acceleration, prev_neighbor_acceleration, mechanics_dt); 
  }
  // late time steps ABM (best would be to update position from force but requires messing with core code at the moment)
  // set previous force to current for calculating previous acceleration
    m_previous_force=m_force;
    m_force={0.0, 0.0, 0.0};
  //
}

void Spring::test_TZPs()
{
  if(!m_is_TZP)
  {
    return;
  }
  else {
    double spring_length=norm(m_neighbor->position-m_me->position)-m_neighbor->phenotype.geometry.radius-m_me->phenotype.geometry.radius;
    double delta_x=spring_length-m_rest_length;
    if(delta_x>parameters.doubles("max_TZP_length"))
    {
      m_spring_constant=0.0;
      m_is_broken=true;
    }
    
  
  }
  
}
void calculate_all_spring_forces()
{
  for(int i=0; i<all_springs.size(); i++)
  {
    Spring* pSpring=all_springs[i];
    pSpring->calculate_spring_force();
  }
  return;
}
void calculate_spring_velocity()
{
  for(int i=0; i<all_springs.size(); i++)
  {
    Spring* pSpring=all_springs[i];
    pSpring->update_spring_velocity();
  }
}
void TZPs()
{
  for(int i=0; i<all_springs.size(); i++)
  {
    Spring* pSpring=all_springs[i];
    pSpring->test_TZPs();
  }
}
// Spring_Connections::Spring_Connections()
// {
//   this->neighbor_springs.clear();
//   this->previous_net_force.resize(3,0.0);
//   this->current_net_force.resize(3,0.0);
//   this->next_velocity.resize(3,0.0);
//   this->current_velocity.resize(3,0.0);
//   this->previous_mass=0;
// }
// void Spring_Connections::sync_cell(Cell* pCell)
// {
//  //done in cryocell.cpp 
//   return;
// }
// void Spring_Connections::spring_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt)
// {
//
//
//   return;
// }
//
// void Spring_Connections::spring_contact_function(Cell* pCell, Phenotype& phenotype, double dt)
// {
//   //double check how this works with pCell->velocity, might be better to just handle everything
//   //if not connected neighbor (sync_connected_neighbors)
//   //adhesive and repulsive with 0 spring length
//
//   return;
// }





/// spring connections
// void Spring_Connections::calculate_spring_force(Spring& spring)
// { 
//   std::cout<< "Spring: "<< &spring<<"\n";
//   // std::vector<double> force(3,0.0);
//   //if connected neighbor
//   double spring_length=norm(spring.m_neighbor->position-this->m_pCell->position)-spring.m_neighbor->phenotype.geometry.radius-this->m_pCell->phenotype.geometry.radius;
//   double delta_x=spring_length-spring.m_rest_length;
//   std::vector<double> unit_vec=(1/norm(spring.m_neighbor->position-this->m_pCell->position))*(spring.m_neighbor->position-this->m_pCell->position); 
//   std::vector<double> force=(spring.m_spring_constant*(delta_x))*unit_vec;
//   spring.m_force=force;
//   // std::cout<<"Spring length: "<< spring_length<<"\n";
//   // std::cout<<"delta_x: "<< delta_x<<"\n";
//   // std::cout<<"unit_vec: "<< unit_vec<< "\n";
//   // std::cout<<"spring constant: "<< spring.m_spring_constant<<"\n";
//   // std::cout<<"force magnitude: "<< spring.m_spring_constant*delta_x<<"\n";
//   // std::cout<<"FORCE: "<< spring.m_force<<"\n";
//   // else displacement==0
//
//   //figure out signs
//   return;
// }
//
// void Spring_Connections::update_springs(std::vector<Spring*> &springs){
//  return; 
// }
// void Spring_Connections::add_spring(Spring* spring)
// {
//   #pragma omp critical
//   {
//     this->neighbor_springs.push_back(spring);
//   }
//   return;
// }
//
// void Spring_Connections::remove_spring(Spring* spring)
// {
//   std::vector<Spring*>::iterator it;
//   int vec_position=0;
//   it = find (this->neighbor_springs.begin(), this->neighbor_springs.end(), spring);
//   if (it != neighbor_springs.end())
//   {
//     vec_position=it-this->neighbor_springs.begin();
// 	  this->neighbor_springs[vec_position]= this->neighbor_springs[neighbor_springs.size()-1];//swap spring locations and pop off the list
// 	  this->neighbor_springs[neighbor_springs.size()-1] = spring;
// 	  this->neighbor_springs.pop_back();
//     delete spring;
//   }
//   else{
//     std::cout << "TRIED TO REMOVE UNKNOWN SPRING\n";
//   }
//   return;
// }

Spring* create_spring( Cell* me, Cell* neighbor, double rest_length, double spring_constant, bool is_TZP ){
  Spring* nSpring=new Spring(me, neighbor, rest_length, spring_constant, is_TZP);
  all_springs.push_back(nSpring);
  return nSpring;
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
