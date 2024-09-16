#include "./cryocell.h"
#include "ABFM.h"
#include "conversions.h"
#include "volume_change.h"
#include <string>
using namespace PhysiCell;
using namespace BioFVM;
#define PI 3.14159265
std::vector<Cryocell*> all_cryocells;
constexpr double GAS_CONSTANT{0.08205};
Cryo_Parameters::Cryo_Parameters()
{
    dVw=0.0;//cell water flux
    previous_dVw=0.0;
  // solute flux
    Ps.resize(parameters.ints("number_of_solutes"),0.0);
    dN.resize(parameters.ints("number_of_solutes"), 0.0);//cell mole flux of solutes
    previous_dN.resize(parameters.ints("number_of_solutes"), 0.0);//previous mole flux of solutes
}
void Cryo_Parameters::sync_to_cell_definition(Cell_Definition& cd)
{
  osmotically_inactive_fraction = cd.custom_data["Vb"];
    // water flux
    Lp=cd.custom_data["Lp"];
    dVw=0.0;//cell water flux
    previous_dVw=0.0;
  // solute flux
    Ps.resize(parameters.ints("number_of_solutes"),0.0);
  for(int i=0; i<parameters.ints("number_of_solutes");i++)
  {
    std::string solute= "Ps_"+std::to_string(i);
    Ps[i]= cd.custom_data[solute];
  };
    dN.resize(parameters.ints("number_of_solutes"), 0.0);//cell mole flux of solutes
    previous_dN.resize(parameters.ints("number_of_solutes"), 0.0);//previous mole flux of solutes
}
Cryo_Concentrations::Cryo_Concentrations()
{
    use_virial = false;
    exterior_osmolality = 0.0;//total exterior osmolality salt+CPA (mole/kg)
    interior_osmolality = 0.0;// total internal osmolality salt+CPA
    interior_molarity.resize(parameters.ints("number_of_solutes"), 0.0 );
    exterior_molarity.resize(parameters.ints("number_of_solutes"), 0.0 );
    interior_component_molality.resize(parameters.ints("number_of_solutes"), 0.0 );
    exterior_component_molality.resize(parameters.ints("number_of_solutes"), 0.0 );
  for(int i =0; i<parameters.ints("number_of_solutes");i++){
     
    std::string density_name=microenvironment.density_names[i];
    // if(density_name=="NaCl")
    // {
      std::string solute= "initial_molarity_"+std::to_string(i);
      interior_molarity[i] = parameters.doubles(solute); //interior molarity is set from user parameters
    // }

  }

}
Cryocell_State::Cryocell_State()
{
    is_cryocell=true;
    previous_radius = 0.0;
    solute_volume = 0.0;
    solute_moles.resize(parameters.ints("number_of_solutes"), 0.0);
    next_solute_moles.resize(parameters.ints("number_of_solutes"), 0.0);
  next_water_volume = 0.0;//for ABM 2nd order
      int num_of_solutes=parameters.ints("number_of_solutes"); 
      double total_solute_volume=0.0;
    //voxel uptakes
    water_uptake = 0.0;//um^3
    solute_uptake.resize(parameters.ints("number_of_solutes"), 0.0);
    solute_uptake_per_voxel.resize(parameters.ints("number_of_solutes"), 0.0);
    water_uptake_per_voxel=0.0;
    uptake.resize(parameters.ints("number_of_solutes"), 0.0);//molar uptake/secretion of solutes
    uptake_voxels={}; //voxels changing from uptake
   
    water_volume = 0.0;

}
void Cryocell_State::sync_to_cell_definition(Cell_Definition& cd, Cryo_Parameters& cryo_p){

  
    previous_radius = 0.0;
  surface_area = 4*3.14159*cd.custom_data["initial_cell_radius"]*cd.custom_data["initial_cell_radius"];
  temperature = cd.custom_data["initial_temp"];
    

      double osmotically_inactive_volume= cd.custom_data["initial_cell_volume"]*cd.custom_data["Vb"];
      int num_of_solutes=parameters.ints("number_of_solutes"); 
    //voxel uptakes
    water_uptake = 0.0;//um^3
    solute_uptake.resize(parameters.ints("number_of_solutes"), 0.0);
    solute_uptake_per_voxel.resize(parameters.ints("number_of_solutes"), 0.0);
    water_uptake_per_voxel=0.0;
    uptake.resize(parameters.ints("number_of_solutes"), 0.0);//molar uptake/secretion of solutes
    uptake_voxels={}; //voxels changing from uptake
   
}

void Cryocell_State::sync_moles_and_volume(Cell_Definition& cd, Cryo_Parameters& cryo_p, Cryo_Concentrations& cc){
    double total_solute_volume=0.0;
    solute_volume = 0.0;
    toxicity=0.0;
    solute_moles.resize(parameters.ints("number_of_solutes"), 0.0);
    next_solute_moles.resize(parameters.ints("number_of_solutes"), 0.0);
    next_water_volume = 0.0;//for ABM 2nd order
    water_volume = (cd.custom_data["initial_cell_volume"]*(1.0-cd.custom_data["Vb"]));
    for(int i=0; i<parameters.ints("number_of_solutes"); i++) //get initial moles from initial molarity
    {
        solute_moles[i]=cc.interior_molarity[i]*water_volume;
    }
    for (size_t i = 0; i < parameters.ints("number_of_solutes"); i++) //get initial solute volume
    {
        std::string density_name=microenvironment.density_names[i];
        total_solute_volume+=moles_to_volume(solute_moles[i], density_name);
    }
    
    solid_volume=cd.custom_data["initial_cell_volume"]*cryo_p.osmotically_inactive_fraction;
    // std::cout<< "custom_data: "<<cd.custom_data["initial_cell_volume"]<<"\n"; 
    // std::cout<< "WATER VOLUME: "<< water_volume<<"\n";   
    // std::cout<< "SOLUTE VOLUME: "<< total_solute_volume<<"\n";   
}
Cryocell::Cryocell(){

  cell_voxels.resize(1,-1);
  neighbor_voxels.resize(1,-1);
  net_force.resize(3, 0.0);
  mass=this->phenotype.volume.total;
    //ghost voxel is place holder for when there are no neighbor_voxels
    
    // ghost_voxel.mesh_index=-1;
    // ghost_voxel.center=this->position;
    // ghost_voxel.volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;

  // std::cout<<"CONSTRUCTED:!! "<< this->position<<this->cell_voxels<<"\n\n\n";
}

void Cryocell::sync_spring_connections(){
  this->spring_connections.m_pCell=this;
}
Cell* instantiate_Cryocell()
{
  Cryocell* cNew=new Cryocell;
  Cell* pNew=static_cast<Cell*>(cNew);
  
  return pNew;

}
Cell* create_Cryocell(Cell_Definition& cd){
  /*Cryocell* cCell=static_cast<Cryocell*>(create_cell());*/
  /*cCell();*/
  Cryocell* cNew=new Cryocell;
  // Cell* pNew;
  all_cryocells.push_back(cNew);
  
  Cell* pNew=static_cast<Cell*>(cNew);
	(*all_cells).push_back( pNew ); 
  cNew->index=(*all_cells).size()-1;

	if( BioFVM::get_default_microenvironment() )
	{
		cNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}

	// All the phenotype and other data structures are already set 
	// by virtue of the default Cell constructor. 
	
	cNew->type = cd.type; 
	cNew->type_name = cd.name; 
	
	cNew->custom_data = cd.custom_data; 
	cNew->parameters = cd.parameters; 
	cNew->functions = cd.functions; 	
	cNew->phenotype = cd.phenotype; 
	if (cNew->phenotype.intracellular)
		cNew->phenotype.intracellular->start();

	cNew->is_movable = cd.is_movable; //  true;
	cNew->is_out_of_domain = false;
	cNew->displacement.resize(3,0.0); // state? 
	
	cNew->assign_orientation();
	cNew->set_total_volume( cNew->phenotype.volume.total );
  cNew->functions.instantiate_cell=NULL; //don't need this atm this thing does all of it

  cNew->cryo_parameters.sync_to_cell_definition(cd);
  cNew->cryocell_state.sync_to_cell_definition(cd, cNew->cryo_parameters);
  cNew->cryocell_state.sync_moles_and_volume(cd, cNew->cryo_parameters, cNew->cryo_concentrations );
  cNew->sync_spring_connections(); 
  return cNew;
}

void Cryocell::update_cell_voxels(){
  std::vector<int> new_voxels{};
  std::vector<int> general_box{};

  std::vector<double> radius(3,this->phenotype.geometry.radius);
  double voxel_length=default_microenvironment_options.dx;
  general_voxel_bounding_box(&general_box, this->position, radius,voxel_length, this->get_microenvironment()->mesh);
  get_intersecting_voxels(this,general_box,&new_voxels);
  // #pragma omp critical
  // {

    this->cell_voxels.clear();
    this->cell_voxels.assign(new_voxels.begin(),new_voxels.end());
    this->cryocell_state.uptake_voxels.assign(new_voxels.begin(),new_voxels.end());///these are the same for now, set in two places
    
  // }
  return;

}
void Cryocell::update_neighbor_voxels(){ //not currently used possible feature for ellipsoid cells
  // 
  // if neighbors.size()==0, this->neighbor_voxels=-1
  // loop through neighbors and get overlapping voxels
  for(int i=0; i<this->state.neighbors.size(); i++){
    Cryocell* cCell=static_cast<Cryocell*>(this->state.neighbors[i]);
    std::vector<int> new_neighbor_voxels{};
    // intersecting_neighbor_voxels(this,cCell,this->cell_voxels,cCell->cell_voxels, &new_neighbor_voxels);
    #pragma omp critical
    {
      this->neighbor_voxels.clear();
      this->neighbor_voxels.assign(new_neighbor_voxels.begin(),new_neighbor_voxels.end());
    }
  }
  return;
}
void update_all_cells_voxels()
{ //TODO: this is slower with parallel??
  // #pragma omp parallel
    // #pragma omp for nowait
    for(int i=0; i<(all_cryocells.size());i++)
    {
      Cryocell* cCell=(all_cryocells)[i];
      cCell->update_cell_voxels();
    }
  /*#pragma omp parallel*/
    /*for(int i=0; i<(all_cryocells.size());i++)*/
    /*{*/
      /*Cryocell* cCell=(all_cryocells)[i];*/
      /*cCell->update_neighbor_voxels();*/
    /*}*/

  //if cell is multivoxel run update_cell_voxells in parallel, then synchronize
  // then run update_neighbor_voxels in parallel
  return;
}
 
void get_concentration_at_boundary()
{
  int num_of_solutes=microenvironment.number_of_densities();
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];   
      std::vector<double> average(num_of_solutes, 0.0);
      std::vector<double> sum(num_of_solutes, 0.0);
      
      for (size_t j = 0; j < cCell->cell_voxels.size(); j++)
      {
        int voxel_index= cCell->cell_voxels[j];
        for (int k =0; k<num_of_solutes; k++)
        {
          std::vector <double> voxel_center=microenvironment.nearest_density_vector(voxel_index);
          if(voxel_center[k]<1e-12)
          {
            sum[k]+=0.0;
          }
          else {
            sum[k]+=voxel_center[k];
          }
        }
      }
      #pragma omp critical
      {
        for (int k =0; k<num_of_solutes; k++)
        {
          average[k]=sum[k]/cCell->cell_voxels.size();
        }
        cCell->cryo_concentrations.exterior_molarity.assign(average.begin(),average.end());
      }
      
    }
  }
  return;
}//currently only 2D, 3D should probably use boundary voxels
void get_exterior_molalities(){

  int num_of_solutes=microenvironment.number_of_densities();
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];
      std::vector<double> exterior_molarity=cCell->cryo_concentrations.exterior_molarity;
      std::vector<double>* exterior_component_molality=&(cCell->cryo_concentrations.exterior_component_molality);
      std::vector<double> temp_molalities(num_of_solutes,0.0);
      double temp_osmolality=0.0;
      for (int k =0; k<num_of_solutes; k++)
      {     
        std::string density_string =cCell->get_microenvironment()->density_names[k]; 
        double molality=molarity_to_molality(exterior_molarity[k],density_string);
        temp_molalities[k]=molarity_to_molality(exterior_molarity[k], density_string);          
      }

      #pragma omp critical
      {
        exterior_component_molality->assign(temp_molalities.begin(),temp_molalities.end());
      }
    }
  }
  return;
}

void get_exterior_osmolalities(){

  int num_of_solutes=microenvironment.number_of_densities();
  #pragma omp parallel 
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];
      double* exterior_osmolality=&(cCell->cryo_concentrations.exterior_osmolality);
      double temp_osmolality=0.0;
      std::vector<double> exterior_component_molality=(cCell->cryo_concentrations.exterior_component_molality);
      if(num_of_solutes<2){
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        temp_osmolality=binary_virial(exterior_component_molality[0] ,sol_name1);
      }
      else if(num_of_solutes<3){
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        std::string sol_name2=cCell->get_microenvironment()->density_names[1];
        temp_osmolality=ternary_virial(exterior_component_molality[0], exterior_component_molality[1], sol_name1,sol_name2);
      }
      else if(num_of_solutes<4)
      {
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        std::string sol_name2=cCell->get_microenvironment()->density_names[1];
        std::string sol_name3=cCell->get_microenvironment()->density_names[2];

        temp_osmolality=binary_virial(exterior_component_molality[0],sol_name1)+ternary_virial(exterior_component_molality[1],exterior_component_molality[2],sol_name2,sol_name3);

      }
      else{
        for (int k =0; k<num_of_solutes; k++)
        {
          temp_osmolality+=exterior_component_molality[k];  
        }            
      }
      #pragma omp critical
      {
        *exterior_osmolality=temp_osmolality;
      }
    }
  }
  return;
}

void get_interior_molalities(){

  int num_of_solutes=microenvironment.number_of_densities();
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];
      std::vector<double> interior_molarity=cCell->cryo_concentrations.interior_molarity;
      std::vector<double>* interior_component_molality=&(cCell->cryo_concentrations.interior_component_molality);
      std::vector<double> temp_molalities(num_of_solutes,0.0);
      double temp_osmolality=0.0;
      for (int k =0; k<num_of_solutes; k++)
      {     
        std::string density_string =cCell->get_microenvironment()->density_names[k]; 
        double molality=molarity_to_molality(interior_molarity[k],density_string);
        temp_molalities[k]=molarity_to_molality(interior_molarity[k], density_string);          
      }

      #pragma omp critical
      {
        interior_component_molality->assign(temp_molalities.begin(),temp_molalities.end());
      }
    }
  }
  return;
}

void get_interior_osmolalities(){

  int num_of_solutes=microenvironment.number_of_densities();
  #pragma omp parallel 
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];
      double* interior_osmolality=&(cCell->cryo_concentrations.interior_osmolality);
      double temp_osmolality=0.0;
      std::vector<double> interior_component_molality=(cCell->cryo_concentrations.interior_component_molality);
      if(num_of_solutes<2){
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        temp_osmolality=binary_virial(interior_component_molality[0] ,sol_name1);
      }
      else if(num_of_solutes<3){
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        std::string sol_name2=cCell->get_microenvironment()->density_names[1];
        temp_osmolality=ternary_virial(interior_component_molality[0], interior_component_molality[1], sol_name1,sol_name2);
      }
      else if(num_of_solutes<4)
      {
        std::string sol_name1=cCell->get_microenvironment()->density_names[0];
        std::string sol_name2=cCell->get_microenvironment()->density_names[1];
        std::string sol_name3=cCell->get_microenvironment()->density_names[2];

        temp_osmolality=binary_virial(interior_component_molality[0],sol_name1)+ternary_virial(interior_component_molality[1],interior_component_molality[2],sol_name2,sol_name3);

      }
      else{
        for (int k =0; k<num_of_solutes; k++)
        {
          temp_osmolality+=interior_component_molality[k];  
        }            
      }
      #pragma omp critical
      {
        *interior_osmolality=temp_osmolality;
      }
    }
  }
  return;
}

void update_exterior_concentrations(){
  get_exterior_molalities();
  get_exterior_osmolalities();
}

void update_interior_concentrations(){
  get_interior_molalities();
  get_interior_osmolalities();
}

void calculate_derivatives(){
  #pragma omp parallel 
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell=all_cryocells[i];
      double Lp = (cCell->cryo_parameters.Lp);
      double surface_area = (cCell->cryocell_state.surface_area);
      double temperature = (cCell->cryocell_state.temperature);
      double exterior_osmolality = (cCell->cryo_concentrations.exterior_osmolality);
      double interior_osmolality = (cCell->cryo_concentrations.interior_osmolality);
      double* previous_dVw = &(cCell->cryo_parameters.previous_dVw);
      double* dVw= &(cCell->cryo_parameters.dVw);
      std::vector <double> Ps = (cCell->cryo_parameters.Ps);
      std::vector <double> exterior_molarity = (cCell->cryo_concentrations.exterior_molarity);
      std::vector <double> interior_molarity = (cCell->cryo_concentrations.interior_molarity);
      std::vector <double>* previous_dN = &(cCell->cryo_parameters.previous_dN);
      std::vector <double>* dN= &(cCell->cryo_parameters.dN);

      dVw_Osmolality(Lp, surface_area, GAS_CONSTANT, temperature, exterior_osmolality, 
                     interior_osmolality, previous_dVw, dVw);// calculate dVw/dt
      dN_molarity(Ps, surface_area, exterior_molarity, interior_molarity, previous_dN, dN);
      
    }
  }
}

void advance_osmosis(double dt){
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      std::vector <double> dN = (cCell->cryo_parameters.dN);
      std::vector <double> previous_dN = (cCell->cryo_parameters.previous_dN);
      double dVw = (cCell->cryo_parameters.dVw);
      double previous_dVw = (cCell->cryo_parameters.previous_dVw);
      double water_volume = (cCell->cryocell_state.water_volume);
      double* next_water_volume = &(cCell->cryocell_state.next_water_volume);
      std::vector <double> solute_moles = (cCell->cryocell_state.solute_moles);
      std::vector <double>* next_solute_moles = &(cCell->cryocell_state.next_solute_moles);
      if (PhysiCell_globals.current_time<dt) {
        Forward_Euler(next_water_volume,water_volume, dVw, dt);
        Forward_Euler_vec(next_solute_moles,solute_moles, dN, dt);
        //cCell->cryocell_state.next_water_volume=cCell->cryocell_state.water_volume+(dVw*dt);//forward_euler
      }
      else {
        /*Forward_Euler(next_water_volume,water_volume, dVw, dt);*/
        /*Forward_Euler_vec(next_solute_moles,solute_moles, dN, dt);*/
        Adams_Bashforth_2(next_water_volume,water_volume, dVw, previous_dVw, dt);
        Adams_Bashforth_2_vec(next_solute_moles,solute_moles, dN, previous_dN, dt);
      }
    }
  }

  return;
}

void calculate_uptakes(double dt){
  //water_uptake = next_water_volume-water_volume;
  //solute_uptake= next_solute_moles-solute_moles;
  #pragma omp parallel 
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      std::vector <int>* uptake_voxels = &(cCell->cryocell_state.uptake_voxels);
      std::vector <int> cell_voxels = (cCell->cell_voxels);
      std::vector <double> *solute_uptake = &(cCell->cryocell_state.solute_uptake);
      double water_volume = (cCell->cryocell_state.water_volume);
      double next_water_volume = (cCell->cryocell_state.next_water_volume);
      std::vector <double> solute_moles = (cCell->cryocell_state.solute_moles);
      std::vector <double> next_solute_moles = (cCell->cryocell_state.next_solute_moles);
      double *water_uptake = &(cCell->cryocell_state.water_uptake);
      #pragma omp critical
      {
        *solute_uptake=next_solute_moles-solute_moles;
        *water_uptake=next_water_volume-water_volume;
        /**uptake_voxels = cell_voxels;*/

      }
    }
  }
  return; 
}

void update_next_step(double dt){
// set current to next
// set previous to current 
// set uptake voxels to my voxels- might update in later version
  #pragma omp parallel 
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];

      std::vector <double> dN = (cCell->cryo_parameters.dN);
      std::vector <double>* previous_dN = &(cCell->cryo_parameters.previous_dN);
      double dVw = (cCell->cryo_parameters.dVw);
      double* previous_dVw = &(cCell->cryo_parameters.previous_dVw);
      double* water_volume = &(cCell->cryocell_state.water_volume);
      double next_water_volume = (cCell->cryocell_state.next_water_volume);
      std::vector<double> *interior_molarity = &(cCell->cryo_concentrations.interior_molarity);
      std::vector <double>* solute_moles = &(cCell->cryocell_state.solute_moles);
      std::vector <double> next_solute_moles = (cCell->cryocell_state.next_solute_moles);
      #pragma omp critical
      {
        /*std::cout<<"water_volume: "<< *water_volume<<"\n";*/
        for(int i=0;i<cCell->cryocell_state.solute_moles.size();i++)
        {
            (*interior_molarity)[i]=(*solute_moles)[i]/(*water_volume);
            /*std::cout<<"Solute moles: "<< (*solute_moles)[i]<<"\n";*/
        }
        *water_volume = next_water_volume;
        *solute_moles = next_solute_moles;
        *previous_dN = dN;
        *previous_dVw = dVw;

      }
    }
  }

  return;
}

void calculate_per_voxel_uptake()
{
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dx*default_microenvironment_options.dx;
      Cryocell* cCell = all_cryocells[i];
      int num_of_solutes=cCell->get_microenvironment()->number_of_densities(); 
      double uptake_voxel_num=cCell->cell_voxels.size();
      double water_volume_uptake= cCell->cryocell_state.water_uptake;
      double water_uptake_per_voxel=water_volume_uptake/(double)uptake_voxel_num;
      std::vector<double> solute_uptake_per_voxel = cCell->cryocell_state.solute_uptake_per_voxel;
      std::vector<double> solute_moles_uptakes= cCell->cryocell_state.solute_uptake; 
      for (size_t i =0; i<num_of_solutes;i++) {
        solute_uptake_per_voxel[i]=solute_moles_uptakes[i]/(double)uptake_voxel_num;
        if(abs(solute_uptake_per_voxel[i])<1e-12)//deal with very tiny values
        {
          solute_uptake_per_voxel[i]=0;
        }
      }
      #pragma omp critical
      {
        cCell->cryocell_state.solute_uptake_per_voxel.assign(solute_uptake_per_voxel.begin(),solute_uptake_per_voxel.end());
        cCell->cryocell_state.water_uptake_per_voxel=water_uptake_per_voxel;
      }
    }// ofs<< solute_uptake_per_voxel[i]<<", ";
  }

  return;
}
/*
void uptake_in_one_voxel(int &voxel, double& water_uptake_per_voxel, std::vector<double>& solute_uptake_per_voxel)
{
  //calculate moles of each solute in the voxel so we can reduce it by the known change in moles
  double voxel_volume=1000;//default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
  std::vector <double> temp_density_vec=microenvironment.density_vector(voxel);
  std::vector <double> moles_in_voxel(solute_uptake_per_voxel.size(),0.0);
  //the volume contribution of solutes compared to the voxel volume is small so we do not include it
  for(int i=0; i<(solute_uptake_per_voxel).size(); i++)
  {
    moles_in_voxel[i]=(temp_density_vec[i]*voxel_volume);//fmole/um^3* um^3
  }
  // voxels don't track water, but we know that the movement of water as bulk solution is at least twice as fast as the movement of solute 
  // since voxels have concentrations similar to their neighbors (assuming they arent on the boundary) we make the amount of water available
  // equal to that in a voxels moore neighborhood or the 3D analog this is similar to the infinite bath argument, only here it is finite but very large
  // a more accurate version would be to track water and/or subtract the volume of solutes
  std::vector <double> new_moles=moles_in_voxel-solute_uptake_per_voxel;
  double effective_water_from_diffusion=24*voxel_volume; //~water from each voxel in my neighborhood (including diagnols)+ my own
  if(default_microenvironment_options.simulate_2D){
    effective_water_from_diffusion=8*voxel_volume; //~water from each voxel in my neighborhood (including diagnols)+ my own
  }

  double new_water_volume=effective_water_from_diffusion-water_uptake_per_voxel;
  if (new_water_volume<0)
  {
    std::cout<<"WARNING!! HIGH WATER UPTAKE!!\n";
    throw std::invalid_argument("Water in voxel went negative!");
  }
  //handle numerical issue with tiny values when differences in concentration are tiny but extreme in early simulation
  for(int i=0; i<new_moles.size(); i++)
  {
    if(abs(new_moles[i])<1e-12)
    {
      new_moles[i]=moles_in_voxel[i];
    }
  }
  std::vector<double> new_density{ 0.0, 0.0, 0.0 };

  for(int i=0; i<new_moles.size();i++)
  {
    new_density[i]=new_moles[i]/new_water_volume;
    if(new_density[i]<0) //cannot take moles that aren't there
    {
      std::cout<<"NEGATIVE MOLES! \n";
      throw std::invalid_argument("Took more solute moles than voxel holds!");
    }
  }

  std::ofstream ofs;
  ofs.open("./output/uptake_in_one_voxel.csv", std::ofstream::out|std::ofstream::app);
  ofs<<PhysiCell_globals.current_time<<", "<< voxel<<", water_uptake_per_voxel: "<< water_uptake_per_voxel<<"\n";
  ofs<<"solute_uptake_per_voxel: "<< solute_uptake_per_voxel <<", ";
  ofs<<"new moles: "<<new_moles<<", "<<"new_water_volume: "<<new_water_volume <<"\n";
  for(int i=0; i<microenvironment.number_of_densities();i++)
  {
    std::cout<<"old density: "<<microenvironment.density_vector(voxel)[i]<<"\n";
  }
  microenvironment.density_vector(voxel)=new_density;//safety shmafety
  for(int i=0; i<microenvironment.number_of_densities();i++)
  {
    std::cout<<"new density: "<<microenvironment.density_vector(voxel)[i]<<"\n";
  }
  ofs<<"\n\n\n";
  ofs.close();
  return;
}
*/
void advance_uptake()
{
  #pragma omp parallel 
  {  
    #pragma omp for nowait
    for(int i = 0; i < all_cryocells.size(); i++) {
      Cryocell* cCell = all_cryocells[i];
      std::vector<int>* uptake_voxels = &(cCell->cryocell_state.uptake_voxels);
      std::vector<double> solute_uptake_per_voxel=cCell->cryocell_state.solute_uptake_per_voxel;
      double water_uptake_per_voxel=cCell->cryocell_state.water_uptake_per_voxel;
      Microenvironment* m= (cCell->get_microenvironment());
      for (int j=0; j<uptake_voxels->size();j++) {
        double voxel_volume=m->voxels((*uptake_voxels)[j]).volume;
        std::vector <double> temp_density_vec=m->density_vector((*uptake_voxels)[j]);
        std::vector <double> moles_in_voxel(solute_uptake_per_voxel.size(),0.0);
        //the volume contribution of solutes compared to the voxel volume is small so we do not include it
        for(int i=0; i<(solute_uptake_per_voxel).size(); i++)
        {
          moles_in_voxel[i]=(temp_density_vec[i]*voxel_volume);//fmole/um^3* um^3
        }
        // voxels don't track water, but we know that the movement of water as bulk solution is at least twice as fast as the movement of solute 
        // since voxels have concentrations similar to their neighbors (assuming they arent on the boundary) we make the amount of water available
        // equal to that in a voxels moore neighborhood or the 3D analog this is similar to the infinite bath argument, only here it is finite but very large
        // a more accurate version would be to track water and/or subtract the volume of solutes
        std::vector <double> new_moles=moles_in_voxel-solute_uptake_per_voxel;
        double effective_water_constant=24;
        if(default_microenvironment_options.simulate_2D){
          effective_water_constant=8; //~water from each voxel in my neighborhood (including diagnols)+ my own
        }

        double effective_water_from_diffusion=effective_water_constant*voxel_volume; //~water from each voxel in my neighborhood (including diagnols)+ my own
        double new_water_volume=(effective_water_from_diffusion-water_uptake_per_voxel)/effective_water_constant;
        // double new_water_volume=(effective_water_from_diffusion-water_uptake_per_voxel)/effective_water_constant;
        if (new_water_volume<0)
        {
          std::cout<<"WARNING!! HIGH WATER UPTAKE!!\n";
          throw std::invalid_argument("Water in voxel went negative!");
        }
        //handle numerical issue with tiny values when differences in concentration are tiny but extreme in early simulation
        for(int i=0; i<new_moles.size(); i++)
        {
          if(abs(new_moles[i])<1e-12)
          {
            new_moles[i]=moles_in_voxel[i];
          }
        }
        std::vector<double> new_density(new_moles.size(), 0.0);

        for(int i=0; i<new_moles.size();i++)
        {
          new_density[i]=(new_moles[i]/new_water_volume);
          /*new_density[i]=(new_moles[i]/new_water_volume);*/
          if(new_density[i]<0) //cannot take moles that aren't there
          {
            std::cout<<"NEGATIVE MOLES! \n";
            throw std::invalid_argument("Took more solute moles than voxel holds!");
          }
        }

        #pragma omp critical
        {
          /*std::cout<<"moles_in_voxel: "<< moles_in_voxel<<"\n";*/
          /*std::cout<<"new density: "<< new_density<<"\n\n";*/
          m->density_vector((*uptake_voxels)[j])=new_density;//safety shmafety
        }

      }
    }
  }
  return;
}

void uptake(double dt)
{

  /*calculate_per_voxel_uptake();*/
  advance_uptake();
  return;

}

  //the following function runs in main.cpp
void two_p_forward_step(double dt){

/*Cryocell* cCell=static_cast<Cryocell*>(pCell);*/
//cCell->update_cell_voxels(); has to happen at the end of time step after movement
//read in locally

//(voxels needed to be updated)
//
//get external Cryo_Concentrations
//do conversions
//calculate dN dVw 
//calculate next water and next solutes 
//calculate uptakes 
//update internal values (convert if needed)
//set current values to previous values
// uptake from voxels
//
get_concentration_at_boundary();
update_exterior_concentrations();//update and get exterior concentrations
update_interior_concentrations();
calculate_derivatives();
advance_osmosis(dt);
calculate_uptakes(dt); //handled externally?
calculate_per_voxel_uptake();// handled externally?
advance_uptake();
update_next_step(dt);
return;
}
void update_multivoxel_neighboorhood(){
 
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      Cell* pCell =static_cast<Cryocell*>(cCell);
      find_multivoxel_neighbors(pCell,&cCell->all_neighbors);
    }
  return;
  
}
void update_initial_neighbors(){

    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      Cell* pCell =static_cast<Cryocell*>(cCell);
      find_multivoxel_neighbors(pCell,&cCell->initial_neighbors);
    }
  return;
}
void update_springs(){

    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      for(int j=0; j<cCell->initial_neighbors.size(); j++)
      {
        bool is_TZP=false;
        Cell* pCell= cCell->initial_neighbors[j];
        double spring_constant=pCell->custom_data["k_granulosa"];
        if(cCell->type_name=="oocyte" || pCell->type_name=="oocyte")
        {
          spring_constant=pCell->custom_data["k_oocyte"];
          is_TZP=true;
        }
        double rest_length=norm(pCell->position-cCell->position)-pCell->phenotype.geometry.radius-cCell->phenotype.geometry.radius;
        Spring sprg(pCell, rest_length, spring_constant, is_TZP);
        cCell->spring_connections.add_spring(&sprg);
      }
    }
}
void update_all_spring_forces(){
    for(int i=0; i<all_cryocells.size(); i++)
    {
      Cryocell* cCell = all_cryocells[i];
      for(int j=0; j<cCell->spring_connections.neighbor_springs.size(); j++)
      {
        Spring* sprg = cCell->spring_connections.neighbor_springs[j];
        cCell->spring_connections.calculate_spring_force(*sprg);

      }
    }

  return;
}
void sum_spring_forces(Cryocell* cCell){
  std::vector<double> sum_forces(3,0.0);
  for(int i=0; i<cCell->spring_connections.neighbor_springs.size();i++)
  {
      sum_forces+=cCell->spring_connections.neighbor_springs[i]->m_force;
  }
  cCell->net_force+=sum_forces;
}
void sum_youngs_modulus(Cryocell* cCell){

  std::vector<double> sum_forces(3,0.0);
  for(int i=0; cCell->all_neighbors.size();i++)
  {
    std::vector<double> temp_force(3,0.0);
    Cell* pNeighbor=cCell->all_neighbors[i];
    cell_to_cell_youngs_modulus(cCell, pNeighbor, &temp_force);
    sum_forces+=temp_force;
  }
  cCell->net_force+=sum_forces;
}
void cell_to_cell_youngs_modulus( Cryocell* pMe, Cell* pOther, std::vector<double> *return_force){
  //assumes only Hertzian, no friction and smooth spherical surfaces 
  std::vector<double> displacement= pMe->position-pOther->position;//force exerted on me
  double penetration_depth= norm(displacement)-pOther->phenotype.geometry.radius-pMe->phenotype.geometry.radius;
  if(penetration_depth>0)
  {return;}
  double constant = std::pow((3*3.14159),(2.0/3.0))/2;
  double E1=pMe->custom_data["youngs_modulus"];
  double E2=pOther->custom_data["youngs_modulus"];
  double sigma_1=pMe->custom_data["poisson_ratio"];
  double sigma_2=pOther->custom_data["poisson_ratio"];
  double V1= (1-sigma_1*sigma_1)/(3.14159*E1);
  double V2= (1-sigma_2*sigma_2)/(3.14159*E2);
  double D1= 1/pMe->phenotype.geometry.radius*2;
  double D2= 1/pOther->phenotype.geometry.radius*2;
  penetration_depth*=-1;
  double Vterm=std::pow((V1+V2),(2.0/3.0));
  double Dterm=std::pow((D1+D2),(2.0/3.0));
  double force_magnitude=penetration_depth/(constant*Vterm*Dterm);
  force_magnitude=std::pow(force_magnitude,(3.0/2.0));
  force_magnitude=force_magnitude/norm(displacement);
  std::vector<double> force= force_magnitude*displacement;
  if(std::fabs(force_magnitude)<1e-16)
  {
    force={0.0, 0.0, 0.0};
  }
  *return_force= force;
  
  return;
}

void update_net_force(){
  for(int i=0; i<all_cryocells.size(); i++)
  {
    Cryocell* cCell = all_cryocells[i];
    sum_spring_forces(cCell);
    sum_youngs_modulus(cCell);
  } 
}

void Cell::update_position( double dt )
{
	// BioFVM Basic_Agent::update_position(dt) returns without doing anything. 
	// So we remove this to avoid any future surprises. 
	// 
	// Basic_Agent::update_position(dt);
		
	// use Adams-Bashforth 
	static double d1; 
	static double d2; 
	static bool constants_defined = false; 
	if( constants_defined == false )
	{
		d1 = dt; 
		d1 *= 1.5; 
		d2 = dt; 
		d2 *= -0.5; 
		constants_defined = true; 
	}
	
	// new AUgust 2017
	if( default_microenvironment_options.simulate_2D == true )
	{ velocity[2] = 0.0; }
	
	std::vector<double> old_position(position); 
	axpy( &position , d1 , velocity );  
	axpy( &position , d2 , previous_velocity );  
	// overwrite previous_velocity for future use 
	// if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
		// std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;
	
	return; 
}
void update_position_from_net_force(Cryocell* cCell, double dt){
  
	previous_velocity = velocity; 
	
	velocity[0]=0; velocity[1]=0; velocity[2]=0;
	if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
	{
		updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	}
	else
	{
		updated_current_mechanics_voxel_index=-1;
		
		is_out_of_domain = true; 
		is_active = false; 
		is_movable = false; 
	}
}
void two_p_update_volume() //TODO: check if parallel is faster
{

  // #pragma omp parallel
  // {
    // #pragma omp for nowait
    for(int i=0; i<all_cryocells.size(); i++)
    {
      double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
      Cryocell* cCell = all_cryocells[i];
      double osmotically_inactive_volume= cCell->custom_data["initial_cell_volume"]*cCell->custom_data["Vb"];
      // std::cout<< "osmotically_inactive_volume: "<< osmotically_inactive_volume<< "\n";
      int num_of_solutes=cCell->get_microenvironment()->number_of_densities(); 
      double total_solute_volume=0.0;
      for (size_t i = 0; i < num_of_solutes; i++)
      {
        std::string density_name=microenvironment.density_names[i];
           total_solute_volume+=moles_to_volume(cCell->cryocell_state.solute_moles[i], density_name);
      }
      // #pragma omp critical
      // {
        cCell->cryocell_state.solute_volume=total_solute_volume;
      
        // std::cout<< "water volume: "<< cCell->cryocell_state.water_volume<<"\n";
        cCell->set_total_volume(cCell->cryocell_state.water_volume+osmotically_inactive_volume+total_solute_volume); 
        // std::cout<< "total volume: "<< cCell->phenotype.volume.total<<"\n";
      // }
    }// ofs<< solute_uptake_per_voxel[i]<<", ";
  // }


  return;
}
