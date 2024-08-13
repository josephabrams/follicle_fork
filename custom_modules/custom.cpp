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
/*#include "./multivoxel/multivoxel_functions.h"*/
#include <ios>
#include <vector>
#include <cmath>
#include "./cryomodule/cryocell.h"
void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	if (PhysiCell::parameters.bools("fibre_custom_degradation"))
		cell_defaults.functions.instantiate_cell = instantiate_physimess_cell_custom_degrade;	
	else
		cell_defaults.functions.instantiate_cell = instantiate_physimess_cell;	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity; 
    /*physimess_update_cell_velocity;*/

	cell_defaults.functions.update_migration_bias = NULL; 
	// cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	// cell_defaults.functions.custom_cell_rule = NULL; 
	// cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
  // cell_defaults.custom_data.add_variable("initial_volume","um^3",0.0);
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

  
	for (auto* pCD: *getFibreCellDefinitions()){
		pCD->functions.instantiate_cell = instantiate_physimess_fibre;
		pCD->functions.plot_agent_SVG = fibre_agent_SVG;
		pCD->functions.plot_agent_legend = fibre_agent_legend;
    pCD->functions.custom_cell_rule=custom_function;
	  pCD->functions.update_velocity = physimess_update_cell_velocity;
  }
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
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
	
  load_cells_from_pugixml();
  //
  //
  bool isFibreFromFile = false;
    
//    for( int i=0; i < (*all_cells).size(); i++ ){

  //      if (isFibre((*all_cells)[i]))
    //    {
            /* fibre positions are given by csv
               assign fibre orientation and test whether out of bounds */
      //      isFibreFromFile = true;
			//static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->assign_fibre_orientation();
			
        //} 
    //}

    /* agents have not been added from the file but do want them
       create some of each agent type */

    if(!isFibreFromFile){
        Cell* pC;
        std::vector<double> position = {0, 0, 0};
        
  double initial_cell_radius=parameters.doubles("average_hepg2_radius"); 
  double initial_inner_radius=parameters.doubles("average_pocket_radius");
	double initial_overlap=.1;//representation of the initial packing density -- taken from physicell simulations of cancer spheroids
          std::vector<double> shift{30.0, 30.0, 0.0};
	double sphere_radius=initial_cell_radius*1; 
	
  double cell_spacing = initial_cell_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
	
    std::vector<std::vector<double>> cell_positions_1= create_spheroid_2D(initial_cell_radius, sphere_radius);//
    std::cout<<"THERE ARE "<< cell_definitions_by_index.size()<< " TYPES OF AGENTS!\n";
    Cell_Definition *pCD0=cell_definitions_by_index[0];
    Cell_Definition *pCD1=cell_definitions_by_index[1];
    Cell_Definition *pCD2=cell_definitions_by_index[2];
        for( int k=0; k < cell_definitions_by_index.size() ; k++ ) {

            Cell_Definition *pCD = cell_definitions_by_index[k];
            // std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
            
            if (!isFibre(pCD))
            {
                
                /*for (int n = 0; n < parameters.ints("number_of_cells"); n++) {*/
//                for(int n=0; n<cell_positions_1.size(); n++)
  //              {
                  /*position[0] = 0 + n*120; */
                  /*position[1] = 0 + n*120; */
                  /*position[2] = 0;*/
                  /*pC=create_Cryocell(*pCD); */
                  /*pC = create_cell( *pCD ); */
                  /*pC->assign_position( (cell_positions_1[n]+shift) );*/
                  /*pC->set_radius(pC->custom_data["initial_cell_radius"]);*/
                    /*position[0] = Xmin + UniformRandom() * Xrange;*/
                    /*position[1] = Ymin + UniformRandom() * Yrange;*/
                    /*position[2] = Zmin + UniformRandom() * Zrange;*/

                    /*pC = create_cell(*pCD);*/
                                        
                    /*pC->assign_position(position);*/
    //            }
               // for(int n=0; n<cell_positions_1.size(); n++)
               // {
                  /*position[0] = 0 + n*120; */
                  /*position[1] = 0 + n*120; */
                  /*position[2] = 0;*/
                  pC=create_Cryocell(*pCD); 
                  /*pC = create_cell( *pCD ); */
                  pC->assign_position( position );
                  pC->set_radius(pC->custom_data["initial_cell_radius"]);
                    /*position[0] = Xmin + UniformRandom() * Xrange;*/
                    /*position[1] = Ymin + UniformRandom() * Yrange;*/
                    /*position[2] = Zmin + UniformRandom() * Zrange;*/

                    /*pC = create_cell(*pCD);*/
                                        
                    /*pC->assign_position(position);*/
                //}
            } 
            
            else if(pCD->name==(*pCD1).name) 
            {
                double disk_edge=200;
                std::cout<<(*pCD1).name<<" being placed!\n";
                for ( int nf = 0 ; nf < parameters.ints("number_of_fibres") ; nf++ ) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = 0;
                    if(norm(position)>sphere_radius+100 && norm(position)<= 300)
                    {

                    pC = create_cell(*pCD1);

                      PhysiMeSS_Fibre* pF = static_cast<PhysiMeSS_Fibre*>(pC);
            pF->set_length(pC->custom_data["fibre_length"]);
                      pF->assign_fibre_orientation();
                      pF->check_out_of_bounds(position);
                      
                      std::cout<<"length: "<<pF->mLength <<"\n";
                      pF->assign_position(position);
                       
                    }
                }
            }

            /*else if(pCD->name==(*pCD2).name) */
            /*{*/
                
                /*std::cout<< (*pCD2).name<<" being placed!\n";*/
                /*for ( int nf = 0 ; nf < parameters.ints("number_of_fibres")*15 ; nf++ ) {*/

                    /*position[0] = Xmin + UniformRandom()*Xrange;*/
                    /*position[1] = Ymin + UniformRandom()*Yrange;*/
                    /*position[2] = 0;*/
                    /*if(norm(position)>sphere_radius+20 && norm(position)<= sphere_radius+50)*/
                    /*{*/

                      /*pC = create_cell(*pCD2);*/
                      /*PhysiMeSS_Fibre* pF = static_cast<PhysiMeSS_Fibre*>(pC);*/
            
                      /*pF->set_length(pC->custom_data["fibre_length"]);*/
                      /*pF->assign_fibre_orientation();*/
                      /*pF->check_out_of_bounds(position);*/
                      
                      /*std::cout<<"length: "<<pF->mLength <<"\n";*/
                      /*pF->assign_position(position);*/
                    /*}*/
                /*}*/
            /*}*/
        }
    }

    remove_physimess_out_of_bounds_fibres();
	/*Cell* pC;*/
	
	/*for( int k=0; k < cell_definitions_by_index.size() ; k++ )*/
	/*{*/
		/*Cell_Definition* pCD = cell_definitions_by_index[k]; */
		/*std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; */
		/*for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )*/
		/*{*/
			/*std::vector<double> position = {0,0,0}; */
			// position[0] = Xmin + UniformRandom()*Xrange; 
			// position[1] = Ymin + UniformRandom()*Yrange; 
			// position[2] = Zmin + UniformRandom()*Zrange; 
			
			/*position[0] = 0 + n*20; */
			/*position[1] = 0 + n*20; */
      /*position[2] = 0;*/
      /*pC=create_Cryocell(*pCD); */
			/*pC = create_cell( *pCD ); */
			/*pC->assign_position( position );*/
      /*pC->set_radius(pC->custom_data["initial_cell_radius"]);*/
      
	
    /*}*/
	/*}*/
	/*std::cout << std::endl; */
	
    for( int i=0; i < (*all_cells).size(); i++ ){

        if (isFibre((*all_cells)[i]))
        {
            /* fibre positions are given by csv
               assign fibre orientation and test whether out of bounds */
      //      isFibreFromFile = true;
      //
			static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->register_fibre_voxels();
			static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->find_agent_voxels();
			static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->find_agent_neighbors();
      if(std::abs(uniform_random())>0.8)
      {
			  static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->add_crosslinks();
      }
        } 
    }
	// // load cells from your CSV file (if enabled)
	// // load_cells_from_pugixml(); 
  //
  return; 
}

std::vector<std::string> paint_by_cell_pressure( Cell* pCell ){

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
		return paint_by_cell_pressure(pCell); 
	} else {
		return paint_by_number_cell_coloring(pCell);
	}
}
std::vector<std::string> my_coloring_function_for_substrate( double concentration, double max_conc, double min_conc )
{ return paint_by_density_percentage( concentration,  max_conc,  min_conc); }
/*std::vector<std::string> my_coloring_function( Cell* pCell )*/
/*{ return paint_by_number_cell_coloring(pCell); }*/

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
  // std::vector<int> test_box;
  return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
  /*if(isFibre(pCell))*/
  /*{*/
    /*#pragma omp critical*/
    /*{*/
      /*std::cout<<"isFibre: "<<isFibre(pCell)<<"\n";*/
      /*PhysiMeSS_Fibre* fCell=static_cast<PhysiMeSS_Fibre*>(pCell);*/
      /*fCell->find_agent_neighbors();*/
      /*fCell->add_crosslinks();*/
      /*std::cout<<"Number of crosslinks: "<< fCell->fibres_crosslinkers.size()<<"\n";*/

    /*}*/
  /*}*/
  /*if(pCell->custom_data["is_cryocell"]==1)*/
  /*{*/
  /*#pragma omp critical*/
  /*{*/
    
    /*std::cout<<"Volume:" <<pCell->phenotype.volume.total<<"\n";*/
    /*Cryocell* cCell=static_cast<Cryocell*>(pCell);*/
    /*std::cout<<"volume: "<<cCell->phenotype.volume.total<<"\n";*/
    /*std::cout<<"number of cryocells: "<<all_cryocells.size()<<"\n";*/
    /*std::cout<<"water volume: "<<cCell->cryocell_state.water_volume<<"\n";*/
    /*std::cout<<"interior molarity: "<< cCell->cryo_concentrations.interior_molarity<<"\n";*/
    /*std::cout<<"interior molality: "<< cCell->cryo_concentrations.interior_component_molality<<"\n";*/
    /*std::cout<<"interior osmolality: "<< cCell->cryo_concentrations.interior_osmolality<<"\n";*/
    /*std::cout<<"exterior molarity: "<< cCell->cryo_concentrations.exterior_molarity<<"\n";*/
    /*std::cout<<"surface_area: "<< cCell->cryocell_state.surface_area<<"\n";*/
    /*std::cout<<"exterior osmolality: "<< cCell->cryo_concentrations.exterior_osmolality<<"\n";*/
    /*std::cout<<"Lp: "<< cCell->cryo_parameters.Lp<<"\n";*/
    /*std::cout<<"Ps: "<< cCell->cryo_parameters.Ps<<"\n";*/
    /*std::cout<<"Next water: "<< cCell->cryocell_state.next_water_volume<<"\n";*/
    /*std::cout<<"Next solute: "<< cCell->cryocell_state.next_solute_moles<<"\n";*/
    /*std::cout<<"number of uptake voxels: "<<cCell->cryocell_state.uptake_voxels.size()<<"\n";*/
  /*std::cout<<"dN: "<< cCell->cryo_parameters.dN<<"\n\n";*/
      /*std::cout<<"dVw: "<<cCell->cryo_parameters.dVw<<"\n\n\n";*/
    /*std::cout<<"solute_uptake: "<<cCell->cryocell_state.solute_uptake<<"\n";*/
    /*std::cout<<"water_uptake: "<<cCell->cryocell_state.water_uptake<<"\n";*/
    /*std::cout<<"solute_uptake_per_voxel: "<<cCell->cryocell_state.solute_uptake_per_voxel<<"\n";*/
    /*std::cout<<"water_uptake_per_voxel: "<<cCell->cryocell_state.water_uptake_per_voxel<<"\n\n\n";*/
      /*std::string plot1="intersecting-voxels-";*/
      /*python_plot_cell_and_voxels(cCell, dt,cCell->cell_voxels,plot1);*/
      /*std::string plot2="uptake-voxels-";*/
      /*python_plot_cell_and_voxels(pCell, dt,cCell->cryocell_state.uptake_voxels,plot2);*/

  /*}*/
  /*}*/
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


Cell* instantiate_physimess_cell() { return new PhysiMeSS_Cell; }
Cell* instantiate_physimess_fibre() { return new PhysiMeSS_Fibre; }
Cell* instantiate_physimess_cell_custom_degrade() { return new PhysiMeSS_Cell_Custom_Degrade; }


void PhysiMeSS_Cell_Custom_Degrade::degrade_fibre(PhysiMeSS_Fibre* pFibre)
{
	// Here this version of the degrade function takes cell pressure into account in the degradation rate
    double distance = 0.0;
    pFibre->nearest_point_on_fibre(position, displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    
    
        // Fibre degradation by cell - switched on by flag fibre_degradation
        double stuck_threshold = PhysiCell::parameters.doubles("fibre_stuck_time");
        double pressure_threshold = PhysiCell::parameters.doubles("fibre_pressure_threshold");
        if (PhysiCell::parameters.bools("fibre_degradation") && (stuck_counter >= stuck_threshold
                                                        || state.simple_pressure > pressure_threshold)) {
            // if (stuck_counter >= stuck_threshold){
            //     std::cout << "Cell " << ID << " is stuck at time " << PhysiCell::PhysiCell_globals.current_time
            //                 << " near fibre " << pFibre->ID  << std::endl;;
            // }
            // if (state.simple_pressure > pressure_threshold){
            //     std::cout << "Cell " << ID << " is under pressure of " << state.simple_pressure << " at "
            //                 << PhysiCell::PhysiCell_globals.current_time << " near fibre " << pFibre->ID  << std::endl;;
            // }
            displacement *= -1.0/distance;
            double dotproduct = dot_product(displacement, phenotype.motility.motility_vector);
            if (dotproduct >= 0) {
                double rand_degradation = PhysiCell::UniformRandom();
                double prob_degradation = PhysiCell::parameters.doubles("fibre_degradation_rate");
                if (state.simple_pressure > pressure_threshold){
                    prob_degradation *= state.simple_pressure;
                }
                if (rand_degradation <= prob_degradation) {
                    //std::cout << " --------> fibre " << (*other_agent).ID << " is flagged for degradation " << std::endl;
                    // (*other_agent).parameters.degradation_flag = true;
                    pFibre->flag_for_removal();
                    // std::cout << "Degrading fibre agent " << pFibre->ID << " using flag for removal !!" << std::endl;
                    stuck_counter = 0;
                }
            }
        }
    // }
}
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

void custom_excretion_function(double excretion_time, double excretion_distrance,double dt){
    std::cout<< "EXCRETION! 1 \n";
	  /*if( fabs( PhysiCell_globals.current_time - excretion_time ) < 0.01 * dt )*/
    /*{*/
      std::cout<< "EXCRETION! 2 \n";
        for(int i=0; i<(*all_cells).size(); i++)
        {
          Cell* pCell=(*all_cells)[i];
          if(!isFibre(pCell))
          {
            std::cout<< "NOT FIBER!\n";
            excrete_ECM(pCell, dt, excretion_distrance);
          }
        }
      
    /*}*/
  return;
}
void excrete_ECM(Cell* pCell, double dt, double excretion_distrance){
  
    Cell_Definition *pCD2=cell_definitions_by_index[2];
    std::cout<< (*pCD2).name<<" being placed!\n";
    for ( int nf = 0 ; nf < 2 ; nf++ ) {
      std::vector<double> position= pCell->position;
      if(UniformRandom()<0.25)
      {
        position[0] = position[0]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[0];
        position[1] = position[1]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[1];
        position[2] = 0;
      }
      else if(UniformRandom()>0.25 && UniformRandom()<0.5){
        position[0] = position[0]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[0]*-1;
        position[1] = position[1]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[1];
        position[2] = 0;
      }
      else if(UniformRandom()>0.50 && UniformRandom()<0.75){
        position[0] = position[0]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[0];
        position[1] = position[1]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[1]*-1;
        position[2] = 0;
      }
      else {
        position[0] = position[0]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[0]*-1;
        position[1] = position[1]+ 25*excretion_distrance*pCell->phenotype.geometry.radius*UniformOnUnitCircle()[1]*-1;
        position[2] = 0;
      

      
      }
      if(norm(position)<(norm(pCell->position)+30) && norm(position)>(norm(pCell->position)+pCell->phenotype.geometry.radius))
      {
          
        #pragma omp critical 
        {
          Cell* pC = create_cell(*pCD2);
          PhysiMeSS_Fibre* pF = static_cast<PhysiMeSS_Fibre*>(pC);
          pF->set_length(pC->custom_data["fibre_length"]);
          pF->assign_fibre_orientation();
          pF->check_out_of_bounds(position);
          /*std::cout<<"length: "<<pF->mLength <<"\n";*/
          pF->assign_position(position);
        }
      }
    }
  
  return;
}
