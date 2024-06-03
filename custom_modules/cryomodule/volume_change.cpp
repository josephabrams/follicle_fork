#include "volume_change.h"
// #include "ABFM.h"


void sum_solute_volume(std::vector<double> &solute_specific_volume, std::vector<double> &solute_moles, std::vector<double> *solute_volume )
{
  // std::cout<<"exterior_osmolality: "<<this->exterior_osmolality<<"\n";
  // std::cout<<"interior_osmolality: "<<this->interior_osmolality<<"\n";
    double total_solute_volume=0.0;
    // #pragma omp reduction(+:total_solute_volume)
    for (size_t i = 0; i < solute_moles.size(); i++)
    {
        total_solute_volume+=(solute_specific_volume[i])*(solute_moles[i]);
    }
    
  return;
}
void sum_total_volume(std::vector<double> &solute_volume, double water_volume, double solid_volume, double* new_volume){
  double total_solute_volume=0.0;
  #pragma omp critical
  {
    // std::cout<<"Solute volume: "<<total_solute_volume<<"\n";
    for (size_t i = 0; i < solute_volume.size(); i++)
    {
        total_solute_volume+=solute_volume[i];
    }
    // std::cout<<"water_volume: "<<this->water_volume<<"\n";
    *new_volume=total_solute_volume+water_volume+solid_volume;
    // std::cout<<"Total volume: "<<*new_volume<<"\n";
  }
}
//Note: this function is using molality
void dVw_Osmolality( double &Lp,  double &surface_area,  double &gas_ant,  double &temperature,  double &exterior_osmolality,  double &interior_osmolality, double *previous_dVw, double *dVw)
{ //-LpART(m_e-m_i)
  *previous_dVw=*dVw;
  *dVw=-1*Lp*surface_area*temperature*gas_ant*(exterior_osmolality-interior_osmolality);//total exterior osmolality salt+CPA (mole/kg)
  // std::cout<<"exterior_osmolality: "<<this->exterior_osmolality<<"\n";
  // std::cout<<"interior_osmolality: "<<this->interior_osmolality<<"\n";
  // std::cout<<"difference: "<< this->exterior_osmolality-this->interior_osmolality<<"\n";
  // std::cout<<"dVw: "<<this->dVw<<"\n";
  // std::cout<<"Lp: "<< this->Lp<<"\n";
  // std::cout<<"temperature: "<< this->temperature<<"\n";
  // std::cout<<"gas_ant: "<< this->gas_ant<<"\n";
  return;
}
//Note: this function is using molarity
void dN_molarity( std::vector<double> &Ps,  double &surface_area,  std::vector<double> &exterior_molarity,  std::vector<double> &interior_molarity, std::vector<double> *previous_dN, std::vector<double> *dN)
{
    //multisolute
    //function=dN/dt=Ps*A(mol^ext-mol^int) for the ith solute, Ps=0 for non-permeating
    //internal molarity depends on Vw
    *previous_dN=*dN; 
    for (size_t i = 0; i < (*dN).size(); i++)
    {
      (*dN)[i]=Ps[i]*surface_area*(exterior_molarity[i]-interior_molarity[i]);
    }
  // std::cout<<"Ps[0]: "<< this->Ps[0]<<"\n";
  // std::cout<<"dN[0]: "<<this->dN[0]<<"\n";
  // std::cout<<"exterior_molarity[0]: "<<this->exterior_molarity[0]<<"\n";
  // std::cout<<"interior_molarity: "<<this->interior_molarity[0]<<"\n";
   // std::cout<<"exterior_molarity 1: "<<this->exterior_molarity[1]<<"\n";
   // std::cout<<"interior_molarity 1: "<<this->interior_molarity[1]<<"\n";
    return;
}

void calculate_solute_moles(std::vector<double>* next_solute_moles,  std::vector<double>& solute_moles,  std::vector<double>& dN,  std::vector<double> previous_dN,  double dt){
  if (PhysiCell_globals.current_time<dt) {
    Forward_Euler_vec(next_solute_moles, solute_moles, dN, dt); 
  }
  else{
    Adams_Bashforth_2_vec(next_solute_moles,solute_moles, dN, previous_dN, dt);
  }
  return;
}
void calculate_water_volume(double* next_water_volume,  double &water_volume,  double &dVw,  double &previous_dVw,  double dt){

  if (PhysiCell_globals.current_time<dt) {
    Forward_Euler(next_water_volume, water_volume, dVw, previous_dVw); 
  }
  else{
    Adams_Bashforth_2(next_water_volume,water_volume, dVw, previous_dVw, dt);
  }
  return;
}
void calculate_solute_uptake(std::vector<double>* solute_uptake,  std::vector<double>& next_solute_moles,  std::vector<double>& solute_moles){
  for(int i=0; i<(*solute_uptake).size(); i++){
    *solute_uptake=next_solute_moles-solute_moles;
  }
  return;
}
void calculate_water_uptake( double *water_uptake,  double& next_water_volume,  double& water_volume)
{
  *water_uptake=next_water_volume-water_volume;
  return;
}
void two_p_forward_step( double dt)
{
  #pragma omp critical
  {
  }
  return;
}

void two_p_update_volume(){
  //pcell.geometry.update();//update surface_area, radius and nuclear radius
  return;
}
