#include "./custom_coloring.h"

using namespace PhysiCell;
using namespace BioFVM;
std::string formatted_seconds_to_HHMMSS( double seconds )
{
	static std::string output; 
	output.resize( 1024 ); 
  int nSeconds = rint(seconds);	
	int nHours = nSeconds / 3600; 
	nSeconds -= nHours*3600; 
	int nMinutes = nSeconds / 60; 
	
	// int nHours = (int) floor( (nMinutes+1e-6) / 60.0 ); // nMinutes / 60;
	double dSeconds = seconds - 60*( nMinutes + 60*nHours ); 
	if( dSeconds < 0 )
	{ dSeconds = 0.0; }
	sprintf( (char*) output.c_str(),"%d hours, %d minutes , and %2.2f seconds",nHours,nMinutes,dSeconds);
	
	return output ;
}

void Custom_SVG_plot( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::vector<std::string> (*substrate_coloring_function)(double, double, double) , void (cell_counts_function)(char*))
{
	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	if(PhysiCell_settings.enable_substrate_plot == true && (*substrate_coloring_function) != NULL){

		double legend_padding = 200.0; // I have to add a margin on the left to visualize the bar plot and the values

		Write_SVG_start( os, plot_width + legend_padding, plot_height + top_margin );

		// draw the background 
		Write_SVG_rect( os , 0 , 0 , plot_width + legend_padding, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

	}
	else{

		Write_SVG_start( os, plot_width , plot_height + top_margin );

		// draw the background 
		Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

	}
	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

  std::string time_label = formatted_seconds_to_HHMMSS( temp_time );
	/*std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); */
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	if (cell_counts_function != NULL){
		cell_counts_function(szString);
	} else {
		sprintf( szString , "%u agents" , total_cell_count ); 
	}
	
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 


	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 

	// used for the legend
	double max_conc;
	double min_conc;
 // color in the background ECM
	if(PhysiCell_settings.enable_substrate_plot == true && (*substrate_coloring_function) != NULL)
	{
		double dz_stroma = M.mesh.dz;

		std::string sub = PhysiCell_settings.substrate_to_monitor;
		int sub_index = M.find_density_index(sub); // check the substrate does actually exist
		if(sub_index == -1){
			std::cout << "ERROR SAMPLING THE SUBSTRATE: COULD NOT FIND THE SUBSTRATE " << sub << std::endl; //if not print error message
		}
		else
		{
			if(PhysiCell_settings.limits_substrate_plot){
			 max_conc = PhysiCell_settings.max_concentration;
			 min_conc = PhysiCell_settings.min_concentration;
			}
			else{
			 max_conc = M.density_vector(5)[sub_index];
			 min_conc = M.density_vector(5)[sub_index];	 // so here I am sampling the concentration to set a min and a mx
			//look for the max and min concentration among all the substrates
			for (int n = 0; n < M.number_of_voxels(); n++)
			{
				double concentration = M.density_vector(n)[sub_index];
				if (concentration > max_conc)
					max_conc = concentration;
				if (concentration < min_conc)
					min_conc = concentration;
			}
			};

			//check that max conc is not zero otherwise it is a big problem!
			if(max_conc == 0){

				max_conc = 1.0;

			};
		
			for (int n = 0; n < M.number_of_voxels(); n++)
			{
				auto current_voxel = M.voxels(n);
				int z_center = current_voxel.center[2];
				double z_displ = z_center -  dz_stroma/2; 
				
				double z_compare = z_displ;

				if (default_microenvironment_options.simulate_2D == true){
				z_compare = z_center;
				};

				if (z_slice == z_compare){			//this is to make sure the substrate is sampled in the voxel visualized (so basically the slice)
					int x_center = current_voxel.center[0];
					int y_center = current_voxel.center[1];
					
					double x_displ = x_center -  dx_stroma/2;
					double y_displ = (y_center - dy_stroma) +  dy_stroma/2;

					double concentration = M.density_vector(n)[sub_index];

					std::vector< std::string > output = substrate_coloring_function(concentration, max_conc, min_conc );

					Write_SVG_rect( os , x_displ - X_lower , y_displ - Y_lower, dx_stroma, dy_stroma , 0 , "none", output[0] );
				}

			}

		}
	}
/* 
 if( ECM.TellRows() > 0 )
 {
  // find the k corresponding to z_slice
  
  
  
  Vector position; 
  *position(2) = z_slice; 
  

  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3
  
  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
  {
   // *position(1) = *Y_environment(j); 
   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size; 
   
   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
   {
    // *position(0) = *X_environment(i); 
    *position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size; 
	
    double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );	
	double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment , position );
	if( isnan( BV ) )
	{ BV = 0.0; }

	vector<string> Colors;
	Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
	Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower, 
	voxel_size , voxel_size , 1 , Colors[0], Colors[0] );
   
   }
  }
 
 }
*/
	os << "  </g>" << std::endl; 
 
	// Now draw vessels

	/*
	 std::vector<std::string> VesselColors = hematoxylin_and_eosin_stroma_coloring( 0,1 );

	 os << " <g id=\"BloodVessels\">" << endl; 
	 extern vector<BloodVesselSegment*> BloodVesselSegments; 
	 Vector Offset; 
	 *Offset(0) = X_lower; 
	 *Offset(1) = Y_lower-top_margin;
	*/
 

 
	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl; 
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			os << "   <g id=\"cell" << pC->ID << "\" " 
			<< "type=\"" << pC->type_name << "\" "; // new April 2022  
			if( pC->phenotype.death.dead == true )
			{ os << "dead=\"true\" " ; } 
			else
			{ os << "dead=\"false\" " ; } 
			os << ">" << std::endl; 
			
			pC->functions.plot_agent_SVG(os, pC, z_slice, cell_coloring_function, X_lower, Y_lower);

			os << "   </g>" << std::endl; 

		}
		
	}
	os << "  </g>" << std::endl; 
	
	// plot intersecting BM points
	/* 
	 for( int i=0 ; i < BasementMembraneNodes.size() ; i++ )
	 {
		// vector<string> Colors = false_cell_coloring( pC ); 
		BasementMembraneNode* pBMN = BasementMembraneNodes[i]; 
		double thickness =0.1; 
		
		if( fabs( *(pBMN->Position)(2) - z_slice ) < thickness/2.0 ) 
		{
		 string bm_color ( "rgb(0,0,0)" );
		 double r = thickness/2.0; 
		 double z = fabs( *(pBMN->Position)(2) - z_slice) ; 

		 os << " <g id=\"BMN" << pBMN->ID << "\">" << std::endl; 
		 Write_SVG_circle( os,*(pBMN->Position)(0)-X_lower, *(pBMN->Position)(1)+top_margin-Y_lower, 10*thickness/2.0 , 0.5 , bm_color , bm_color ); 
		 os << " </g>" << std::endl;
		}
		// pC = pC->pNextCell;
	 }
	*/ 
	
	// end of the <g ID="tissue">
	os << " </g>" << std::endl; 
 
	// draw a scale bar
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	// convert from micron to mm
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// plot runtime 
	szString = new char [1024]; 
	RUNTIME_TOC(); 
	std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );

	if(substrate_coloring_function != NULL){

		// add legend for the substrate

		double conc_interval = (max_conc - min_conc) / 10; // setting the interval for the values in the legend.
	
		szString = new char [1024]; 
		double upper_left_x = plot_width + 25.0;
		double sub_rect_height = (plot_height - 25.0) / 10.0;
		for(int i = 0; i <= 9; i++){ //creating 10 rectangoles for the bar, each one with a different shade of color.

			double concentration_sample = min_conc + (conc_interval * (9-i)); // the color depends on the concentration, starting from the min concentration to the max (which was sampled before)

			std::vector< std::string > output = substrate_coloring_function(concentration_sample, max_conc, min_conc );

			double upper_left_y = sub_rect_height * i; // here I set the position of each rectangole

			Write_SVG_rect(os, upper_left_x, top_margin + upper_left_y, 25.0, sub_rect_height, 0.002 * plot_height , "none", output[0]); //drawing each piece of the barplot

			if(i%2 != 0){ // of course I am not printing each value of the barplot, otherwise is too crowded, so just one each 2

				sprintf( szString , " %.2g", concentration_sample);
				Write_SVG_rect(os, upper_left_x + 25, top_margin + upper_left_y + sub_rect_height - (0.001 * plot_height), 3, 0.002 * plot_height, 0 , "rgb(0,0,0)", "rgb(0,0,0)");
				Write_SVG_text( os , szString, upper_left_x + 28, top_margin + upper_left_y + sub_rect_height, font_size , 
					PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); // misterious values set with a trial and error approach due to OCD. But now the legend is coherent at pixel level
			}
		}

		sprintf( szString , "%.2g", max_conc);
			
		Write_SVG_rect(os, upper_left_x + 25, top_margin - (0.001 * plot_height), 3, 0.002 * plot_height, 0 , "rgb(0,0,0)", "rgb(0,0,0)");
		Write_SVG_text( os , szString, upper_left_x + 28, top_margin, font_size , 
			PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); // misterious values set with a trial and error approach due to OCD. But now the legend is coherent at pixel level

		delete [] szString;
	}
	
	Write_SVG_rect(os, 25.0 + plot_width, top_margin, 25.0, plot_height - 25, 0.002 * plot_height , "black", "none"); // nice black contour around the legend
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}

void Cryo_agent_SVG(std::ofstream& os, PhysiCell::Cell* pC, double z_slice, std::vector<std::string> (*cell_coloring_function)(Cell*), double X_lower, double Y_lower) {

	double r = pC->phenotype.geometry.radius ; 
	double rn = pC->phenotype.geometry.nuclear_radius ; 
	double z = fabs( (pC->position)[2] - z_slice) ; 

	std::vector<std::string> Colors = cell_coloring_function( pC ); 
	
	// figure out how much of the cell intersects with z = 0
	double plot_radius = sqrt( r*r - z*z );

	// then normal cell, plot sphere if it intersects z = 0;
	Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
						plot_radius , 0.5, Colors[1], Colors[0] );
	// plot the nucleus if it, too intersects z = 0;
	if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
	{
		plot_radius = sqrt( rn*rn - z*z );
		Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
							plot_radius, 0.5, Colors[3],Colors[2]);
	}
}

