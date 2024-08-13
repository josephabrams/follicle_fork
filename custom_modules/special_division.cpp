#include "../core/PhysiCell.h"
using namespace PhysiCell;
using namespace BioFVM;
Cell* special_divide(Cell* pCell )
{
	
	// make sure ot remove adhesions 
	remove_all_attached_cells(); 
	remove_all_spring_attachments(); 

	for( int nn = 0 ; nn < custom_data.variables.size() ; nn++ )
	{
		if( custom_data.variables[nn].conserved_quantity == true )
		{ custom_data.variables[nn].value *= 0.5; }
	}
	for( int nn = 0 ; nn < custom_data.vector_variables.size() ; nn++ )
	{
		if( custom_data.vector_variables[nn].conserved_quantity == true )
		{ custom_data.vector_variables[nn].value *= 0.5; }
	}

	
	Cell* child = create_cell(functions.instantiate_cell);
	child->copy_data( this );	
	child->copy_function_pointers(this);
	child->parameters = parameters;
	
	// evenly divide internalized substrates 
	// if these are not actively tracked, they are zero anyway 
	*internalized_substrates *= 0.5; 
	*(child->internalized_substrates) = *internalized_substrates ; 
	
	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );
	
	// randomly place the new agent close to me, accounting for orientation and 
	// polarity (if assigned)
		
	// May 30, 2020: 
	// Set cell_division_orientation = LegacyRandomOnUnitSphere to 
	// reproduce this code 
	/*
	double temp_angle = 6.28318530717959*UniformRandom();
	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
	
	double radius= phenotype.geometry.radius;
	std::vector<double> rand_vec (3, 0.0);
	
	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
	rand_vec[2]= cos( temp_phi );
	
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
	
	if( norm(rand_vec) < 1e-16 )
	{
		std::cout<<"************ERROR********************"<<std::endl;
	}
	normalize( &rand_vec ); 
	rand_vec *= radius; // multiply direction times the displacement 
	*/
	
	std::vector<double> rand_vec = cell_division_orientation(); 
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;	
	rand_vec *= phenotype.geometry.radius;

	child->assign_position(position[0] + rand_vec[0],
						   position[1] + rand_vec[1],
						   position[2] + rand_vec[2]);
						 
	//change my position to keep the center of mass intact 
	// and then see if I need to update my voxel index
	static double negative_one_half = -0.5; 
	axpy( &position, negative_one_half , rand_vec ); // position = position - 0.5*rand_vec; 

	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
	if( !get_container()->underlying_mesh.is_position_valid(position[0], position[1], position[2]))
	{
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
	}	
	 
	update_voxel_in_container();
	phenotype.volume.divide(); 
	child->phenotype.volume.divide();
	child->set_total_volume(child->phenotype.volume.total);
	set_total_volume(phenotype.volume.total);
	
	// child->set_phenotype( phenotype ); 
	child->phenotype = phenotype; 

    if (child->phenotype.intracellular){
        child->phenotype.intracellular->start();
		child->phenotype.intracellular->inherit(this);
	}
// #ifdef ADDON_PHYSIDFBA
// 	child->fba_model = this->fba_model;
// #endif


	// changes for new phenotyp March 2022
	state.damage = 0.0; 
	state.total_attack_time = 0; 
	child->state.damage = 0.0; 
	child->state.total_attack_time = 0.0; 

	return child;
}
