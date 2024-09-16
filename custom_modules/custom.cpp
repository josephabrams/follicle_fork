/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
// #include "./multivoxel/multivoxel_functions.h"
// #include <ios>
#include <vector>
#include <cmath>
#include "./cryomodule/cryocell.h"
// #include "multivoxel/multivoxel_neighborhood.h"
#include "spring_class/spring_class.h"
#include "tissue_construction/tissue_construction.h"
void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	//    Put any modifications to default cell definition here if you 
	//    want to have "inherited" by other cell types. 
	//
	//    This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = NULL;
	// cell_defaults.functions.update_velocity = NULL; 
    /*physimess_update_cell_velocity;*/

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL;  
	cell_defaults.functions.custom_cell_rule = custom_function; 
	// cell_defaults.functions.contact_function = custom_contact_function; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
  // cell_defaults.custom_data.add_variable("initial_volume","um^3",0.0);
	initialize_cell_definitions_from_pugixml(); 
	build_cell_definitions_maps(); 
	setup_signal_behavior_dictionaries(); 	
	setup_cell_rules(); 
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	// cell_defaults.functions.custom_cell_rule = custom_function; 
	// cell_defaults.functions.contact_function = contact_function; 
  // should be instantiated correctly on creation	
  // cell_defaults.functions.instantiate_cell=instantiate_Cryocell;
	// cell_defaults.functions.custom_cell_rule = custom_function; 
  display_cell_definitions( std::cout ); 
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}
void setup_tissue( void )
{
  double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax-Xmin; 
	double Yrange = Ymax-Ymin; 
	double Zrange = Zmax-Zmin; 
	
	// create some of each type of cell 
	
  // load_cells_from_pugixml();
  //
  //

  std::vector <double> initial_velocity={0.0,0.0,0.0};
  double initial_granulosa_radius=parameters.doubles("average_granulosa_radius"); 
  double initial_oocyte_radius=parameters.doubles("average_oocyte_radius");
	double initial_overlap=.1;//representation of the initial packing density -- taken from physicell simulations of cancer spheroids
	double follicle_radius =initial_oocyte_radius+3*(initial_granulosa_radius*2)+initial_granulosa_radius; // 3 to 4 layers
	
  double cell_spacing = initial_granulosa_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
  std::vector<double> center_position{0.0, 0.0, 0.0};
  // std::vector<std::vector<double>> cell_positions_1= create_spheroid_2D(initial_cell_radius, sphere_radius);//
  // std::vector<std::vector<double>> cell_positions= create_spherical_shell(cell_spacing, follicle_radius, initial_oocyte_radius);//
  std::vector<std::vector<double>> cell_positions= twoD_symmetric_test_cells(65);//
  
  std::cout<<"THERE ARE "<< cell_definitions_by_index.size()<< " TYPES OF AGENTS!\n";
  for( int k=0; k < cell_definitions_by_index.size() ; k++ ) 
  {
      Cell_Definition *pCD = cell_definitions_by_index[k];
      std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
      if(pCD->name=="oocyte")
      {
          Cell* pC_oocyte; 
          pC_oocyte=create_Cryocell(*pCD); 
          pC_oocyte->assign_position( center_position );
          pC_oocyte->set_radius(pC_oocyte->custom_data["initial_cell_radius"]);
          pC_oocyte->velocity=initial_velocity;
          
          // std::cout<< pC_oocyte->custom_data["initial_cell_radius"]<<"\n";
      }
      if(pCD->name=="granulosa") {
          for(int i=0; i<cell_positions.size(); i++)
          {
            Cell* pC_granulosa; 
            pC_granulosa=create_Cryocell(*pCD); 
            pC_granulosa->assign_position( cell_positions[i] );
            pC_granulosa->set_radius(initial_granulosa_radius);
            pC_granulosa->velocity=initial_velocity;
            Cryocell* cCell=static_cast<Cryocell*>(pC_granulosa);
              void (*func)(Cell*, Phenotype&, double);
            // std::cout<<cCell<< " and  "<< pC_granulosa;
            // for(int i=0; i<pC_granulosa->custom_data.variables.size();i++)
            // { 
            //   std::cout<< pC_granulosa->custom_data[i]<<"\n";
            //
            // }
          }
      // 
      }
  } 
            
	
  return; 
}
void custom_contact_function(Cell* pCell, Phenotype& phenotype,Cell* pCell_neighbor, Phenotype& neighbor_phenotype, double dt)
{

  Cryocell* cCell=static_cast<Cryocell*>(pCell);
  cCell->spring_connections.spring_contact_function(pCell, phenotype,dt);
  return;
}
std::vector<std::string> paint_by_volume( Cell* pCell ){

  Cryocell* cCell=static_cast<Cryocell*>(pCell);
	std::vector< std::string > output( 0);
	int color = (int) round( ((cCell->phenotype.volume.total) / 10) * 255 );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,0)", 255 - color);
	output.push_back( std::string("red") );
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );
	return output;
}
std::vector<std::string> paint_by_osmolality( Cell* pCell ){

  Cryocell* cCell=static_cast<Cryocell*>(pCell);
	std::vector< std::string > output( 0);
	int color = (int) round( ((cCell->cryo_concentrations.interior_osmolality) / 10) * 255 );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,0)", 255 - color);
	output.push_back( std::string("red") );
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );
	return output;
}
std::vector<std::string> paint_by_cell_pressure( Cell* pCell ){

  Cryocell* cCell=static_cast<Cryocell*>(pCell);
	std::vector< std::string > output( 0);
	int color = (int) round( ((pCell->state.simple_pressure) / 10) * 255 );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,0)", 255 - color);
	output.push_back( std::string("red") );
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );
	return output;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ 
	if (parameters.bools("color_cells_by_pressure")){
		return paint_by_osmolality(pCell); 
		// return paint_by_cell_pressure(pCell); 
	} else {
		return paint_by_number_cell_coloring(pCell);
	}
}
std::vector<std::string> my_coloring_function_for_substrate( double concentration, double max_conc, double min_conc )
{
  return paint_by_density_percentage( concentration,  max_conc,  min_conc); 
}
/*std::vector<std::string> my_coloring_function( Cell* pCell )*/
/*{ return paint_by_number_cell_coloring(pCell); }*/

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // std::vector<int> test_box;
  return; }
void test_function()
{
  return;
}
void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
  if(pCell->custom_data["is_cryocell"]==1 && pCell->type_name=="oocyte")
  {
    #pragma omp critical
    {
      
      // std::cout<<"Volume:" <<pCell->phenotype.volume.total<<"\n";
      Cryocell* cCell=static_cast<Cryocell*>(pCell);
      // std::cout<<"cell type: "<<cCell->type_name<<"\n";
      // std::cout<<"volume: "<<cCell->phenotype.volume.total<<"\n";
      // std::cout<<"number of cryocells: "<<all_cryocells.size()<<"\n";
      // std::cout<<"water volume: "<<cCell->cryocell_state.water_volume<<"\n";
      // std::cout<<"interior molarity: "<< cCell->cryo_concentrations.interior_molarity<<"\n";
      // std::cout<<"interior molality: "<< cCell->cryo_concentrations.interior_component_molality<<"\n";
      // std::cout<<"interior osmolality: "<< cCell->cryo_concentrations.interior_osmolality<<"\n";
      // std::cout<<"exterior molarity: "<< cCell->cryo_concentrations.exterior_molarity<<"\n";
      // std::cout<<"surface_area: "<< cCell->cryocell_state.surface_area<<"\n";
      // std::cout<<"exterior osmolality: "<< cCell->cryo_concentrations.exterior_osmolality<<"\n";
      // std::cout<<"Lp: "<< cCell->cryo_parameters.Lp<<"\n";
      // std::cout<<"Ps: "<< cCell->cryo_parameters.Ps<<"\n";
      // std::cout<<"Next water: "<< cCell->cryocell_state.next_water_volume<<"\n";
      // std::cout<<"Next solute: "<< cCell->cryocell_state.next_solute_moles<<"\n";
      // std::cout<<"number of uptake voxels: "<<cCell->cryocell_state.uptake_voxels.size()<<"\n";
      // std::cout<<"dN: "<< cCell->cryo_parameters.dN<<"\n\n";
      // std::cout<<"dVw: "<<cCell->cryo_parameters.dVw<<"\n\n\n";
      // std::cout<<"solute_uptake: "<<cCell->cryocell_state.solute_uptake<<"\n";
      // std::cout<<"water_uptake: "<<cCell->cryocell_state.water_uptake<<"\n";
      // std::cout<<"solute_uptake_per_voxel: "<<cCell->cryocell_state.solute_uptake_per_voxel<<"\n";
      // std::cout<<"water_uptake_per_voxel: "<<cCell->cryocell_state.water_uptake_per_voxel<<"\n\n\n";
      std::string plot1="neighbor_plot-";
      // python_plot_cell_and_voxels(cCell, dt,cCell->cell_voxels,plot1);
      // std::string plot2="uptake-voxels-";
      // python_plot_cell_and_voxels(pCell, dt,cCell->cryocell_state.uptake_voxels,plot2);

	    // double Zmin = microenvironment.mesh.bounding_box[2]; 
	    // double Zmax = microenvironment.mesh.bounding_box[5]; 
      double z_height=microenvironment.mesh.dz/2;
      double Zmin= 4*z_height*-1; 
      double Zmax= 4*z_height; 
      double zz=Zmin;
      std::vector<Cell*>nn=cCell->all_neighbors;
      std::cout<<"NEIGHBORS SIZE: "<< nn.size()<<"\n";
      std::cout<<"SPRING CONNECTIONS SIZE: "<< cCell->spring_connections.neighbor_springs.size()<<"\n";
      while(zz<Zmax)
      {
        // python_plot_cell_and_voxels_single_layer(cCell,dt, cCell->cell_voxels, plot1, zz);
          python_plot_cell_with_Neighbors(cCell, dt, cCell->cell_voxels, plot1, zz, nn);       
        zz+=z_height;
      }
    }
  }
/*{*/
/*  std::vector<int> test_box{};*/
/*  std::vector<int> test_box2{};*/
/**/
/*  std::vector<int> test_box3{};*/
/*    // diffusion_bounding_box(pCell, &test_box);*/
/*  std::vector<double> radius(3,pCell->phenotype.geometry.radius);*/
/*  double voxel_length=default_microenvironment_options.dx;*/
/*  general_voxel_bounding_box(&test_box, pCell->position, radius,voxel_length, pCell->get_microenvironment()->mesh);*/
/*  std::vector<int> return_box{}; */
/*  get_intersecting_voxels(pCell,test_box,&return_box);*/
/*  std::string plot2="intersecting-neighbour-voxels-";*/
/*  std::vector<int> neighbor_voxels{};*/
/*  Cell* me=(*all_cells)[0];*/
/*  Cell* neighbor=(*all_cells)[1];*/
/**/
/*  general_voxel_bounding_box(&test_box2, me->position, radius,voxel_length, me->get_microenvironment()->mesh);*/
/*  general_voxel_bounding_box(&test_box3, neighbor->position, radius,voxel_length, neighbor->get_microenvironment()->mesh);*/
/*  intersecting_neighbor_voxels(me, neighbor, test_box2,test_box3, &neighbor_voxels);*/
/**/
  /*#pragma omp critical*/
/*  python_plot_two_cells_and_voxels(me,neighbor, dt,neighbor_voxels,plot2);*/
return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


bool arrest(Cell* pCell, Phenotype &phenotype, double dt){
  return true;
}
void custom_arrest_function(double arrest_time, double dt){
  #pragma omp parrallel
  {
    if(PhysiCell_globals.current_time>arrest_time)
    {
      #pragma omp for nowait
      for(int i=0; i<all_cryocells.size(); i++)
      {
        Cell* pCell=static_cast<Cell*>(all_cryocells[i]);
        int j= pCell->phenotype.cycle.data.current_phase_index;
        
        for(int k=0; k< pCell->phenotype.cycle.model().phase_links.size(); k++)
        {
          pCell->phenotype.cycle.model().phase_links[j][k].arrest_function=arrest;
        }
      }
    }
  }
  return;
}
void force_of_youngs_modulus( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt, std::vector<double>* return_force )
{
  //this version assumes pOther is rigid like the probe of an AFM
  std::vector<double> displacement= pMe->position-pOther->position;//force exerted on me
  double penetration_depth= norm(displacement)-pOther->phenotype.geometry.radius-pMe->phenotype.geometry.radius;
  if(penetration_depth>0)
  {return;}
  
    penetration_depth*=-1;
    double R1= pOther->phenotype.geometry.radius+pMe->phenotype.geometry.radius;
    double R2= pMe->phenotype.geometry.radius*pMe->phenotype.geometry.radius;
    double E=pMe->custom_data["youngs_modulus"];
    double force_magnitude= std::pow(penetration_depth,1.5)*std::sqrt(R1/R2)*E*(16.0/9.0);
  force_magnitude=force_magnitude/norm(displacement);
  std::vector<double> force= force_magnitude*displacement;
  if(penetration_depth>0.1*pMe->phenotype.geometry.radius)
  {
    std::cout<< "WARNING! HERTZ YOUNGS MODULUS IS NOT LINEAR AT THIS DEPTH\n";
  }
  *return_force= force;
  return;
}

void caculate_position_from_acceleration(std::vector<double> &old_position, std::vector<double>&current_position, std::vector<double> &net_acceleration, double dt, std::vector<double> *new_position){
  //new_position=2*current_position-old_position+acceleration*dt^2
  if(std::fabs(norm(net_acceleration))<1e-16)
  {
    return;
  }
  else {
    std::vector<double> temp_position=2*current_position;
    temp_position=temp_position+(*new_position)-old_position;
    net_acceleration=dt*dt*net_acceleration;
    temp_position=temp_position+net_acceleration;
    (*new_position)=temp_position;
    
  }

}
