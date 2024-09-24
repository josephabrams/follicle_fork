
#include "./multivoxel_neighborhood.h"
#include "multivoxel_functions.h"
#include <algorithm>
using namespace PhysiCell;
using namespace BioFVM;
void find_multivoxel_neighbors(Cell* pCell, std::vector<Cell*> *return_neighbors){
  
  std::vector<int> general_box{};

  std::vector<double> radius(3,pCell->phenotype.geometry.radius);
  double voxel_length=default_microenvironment_options.dx;
  
  std::vector<double> voxel_size{microenvironment.mesh.dx,microenvironment.mesh.dy,microenvironment.mesh.dz};
  std::vector<double> neighboorhood_radius= radius+voxel_size;
  general_voxel_bounding_box(&general_box, pCell->position, neighboorhood_radius,voxel_length, pCell->get_microenvironment()->mesh);
  
  std::vector<Cell *> agents_in_voxel={};
  // std::cout<<"microenvironment mesh voxels size: "<< pCell->get_microenvironment()->mesh.voxels.size()<<"\n"; 
  // std::cout<<"microenvironment connected_voxel_indices size: "<< pCell->get_microenvironment()->mesh.connected_voxel_indices.size()<<"\n"; 
  // std::cout<<"microenvironment connected_voxel_indices size: "<< pCell->get_microenvironment()->cartesian_indices(i)<<"\n"; 
  for(int i=0; i<general_box.size(); i++)
  {
    //convert to agent grid voxel instead of microenvironment voxel
    int grid_voxel= pCell->get_container()->underlying_mesh.nearest_voxel_index(microenvironment.voxels(general_box[i]).center);
    // std::cout<<"GENARAL BOX VOXEL: "<< general_box[i]<<"\n";
    // std::cout<<"AGENT GRID VOXEL: "<< grid_voxel <<"\n"; 
    //   std::cout<< "microenvironment position: "<< microenvironment.voxels(general_box[i]).center<<"\n";
    //   std::cout<< "agent_grid position: "<< pCell->get_container()->underlying_mesh.voxels[grid_voxel].center<<"\n";
    //   std::cout<< "agent_grid size: "<< pCell->get_container()->agent_grid[grid_voxel].size()<<"\n";
  // // #pragma omp private(agents_in_voxel)
      for (int j = 0; j < pCell->get_container()->agent_grid[grid_voxel].size(); j++)//agent grid holds the cells in each cartesian voxel 
      {
        Cell* temp_agent_ptr=pCell->get_container()->agent_grid[grid_voxel][j];
        // if not pCell and close enough to be a neighbor 
        // if ( temp_agent_ptr!= pCell )
        if ( temp_agent_ptr!= pCell && norm(temp_agent_ptr->position-pCell->position)<=(pCell->phenotype.geometry.radius+pCell->custom_data["maximum_interaction_distance"])) 
        {
          agents_in_voxel.push_back(temp_agent_ptr);
          // std::cout<< i << " Cell: "<<
          // agents_in_voxel[i]->position<<std::endl;
        }
          // std::cout<<"AGENT GRID VOXEL: "<< grid_voxel <<"\n"; 
          // std::cout<< "agent_grid size: "<< pCell->get_container()->agent_grid[grid_voxel].size()<<"\n";
      }
  //
  }

  std::sort(agents_in_voxel.begin(), agents_in_voxel.end());
  auto sorted_duplicates=std::unique(agents_in_voxel.begin(),agents_in_voxel.end());
  agents_in_voxel.erase(sorted_duplicates,agents_in_voxel.end());
  // std::cout<<"pCell: "<<pCell<<"\n";
  // std::cout<<"type: "<<pCell->type_name<<"\n";
  // std::cout<< "possible neighbor SIZE: "<< agents_in_voxel.size()<<"\n";
  // #pragma omp critical
  // {
    return_neighbors->assign(agents_in_voxel.begin(), agents_in_voxel.end());
  // std::cout<< "RETURN NEIGHBOR SIZE: "<< (*return_neighbors).size()<<"\n";

  // }
  return;
}

void find_multivoxel_neighbors_direct_contact(Cell* pCell, std::vector<Cell*> *return_neighbors){
  
  std::vector<int> general_box{};

  std::vector<double> radius(3,pCell->phenotype.geometry.radius);
  double voxel_length=default_microenvironment_options.dx;
  
  std::vector<double> voxel_size{microenvironment.mesh.dx,microenvironment.mesh.dy,microenvironment.mesh.dz};
  std::vector<double> neighboorhood_radius= radius+voxel_size;
  general_voxel_bounding_box(&general_box, pCell->position, neighboorhood_radius,voxel_length, pCell->get_microenvironment()->mesh);
  
  std::vector<Cell *> agents_in_voxel={};
  for(int i=0; i<general_box.size(); i++)
  {
    
  // #pragma omp private(agents_in_voxel)
      for (int j = 0; j < pCell->get_container()->agent_grid[general_box[i]].size(); j++)//each voxel has an agent grid that holds the cells in that voxel 
      {
        Cell* temp_agent_ptr=pCell->get_container()->agent_grid[general_box[i]][j];
        //if not pCell and close enough to be a neighbor 
        if ( temp_agent_ptr!= pCell && norm(temp_agent_ptr->position-pCell->position)<=(pCell->phenotype.geometry.radius)) 
        {
          agents_in_voxel.push_back(temp_agent_ptr);
          // std::cout<< i << " Cell: "<<
          // agents_in_voxel[i]->position<<std::endl;
        }

      }

  }
  // #pragma omp critical
  // {
    return_neighbors->insert(return_neighbors->end(), agents_in_voxel.begin(), agents_in_voxel.end());
  // }
  return;
}
// void update_multivoxel_neighbors(){
// // handled in cryocell for the moment but should eventually be here for flushed out addon
//   //for all multivoxel cells -- multivoxel neighbors
//   return;
// }


void python_plot_cell_with_Neighbors(Cell* pCell, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name, double z_height, std::vector<Cell*> &neighbors)
{
  // if(PhysiCell_globals.current_time>1 && PhysiCell_globals.current_time<1.2)
  if(PhysiCell_globals.current_time>0)
  {
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/cell-"+plot_name+std::to_string(pCell->index)+std::to_string(z_height)+"-layer.py";
    // std::vector <double> cell_plot_position={50,50,50};
    // std::vector <double> translation_vec=pCell->position-cell_plot_position;
    std::ofstream ofs;
    ofs.open (file, std::ofstream::out | std::ofstream::trunc);
    ofs<<"import numpy as np\n"<<"import matplotlib.pyplot as plt\n"<<"import matplotlib.patches as mpatches\n";
    ofs<<"xy_artists = [\n";
    // std::vector<int> bounding_box{0};
    std::vector<double> radial_dimensions{pCell->phenotype.geometry.radius,pCell->phenotype.geometry.radius,pCell->phenotype.geometry.radius};
    std::vector<double> cell_position=pCell->position;
    double voxel_size= default_microenvironment_options.dx;
    BioFVM::Cartesian_Mesh the_mesh=pCell->get_microenvironment()->mesh;
    // general_voxel_bounding_box_3D(&bounding_box_by_index, cell_position, radial_dimensions,voxel_size,the_mesh);
    // diffusion_bounding_box(pCell, &bounding_box);
    for(int i=0; i<bounding_box_by_index.size(); i++)
    {

      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      std::vector <double> voxel_position=microenvironment.voxels(bounding_box_by_index[i]).center;
      // ofs<<"\t#("<<voxel_position[0]<<", "<<voxel_position[1]<<", "<< voxel_position[2] <<"), \n";
      if(std::abs(voxel_position[2]-z_height)<=0.1)
      {
        std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
        ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'), \n";
      }
    }
    double adjusted_radius=intersection_of_cell_and_plane(z_height, pCell->phenotype.geometry.radius);
    if(std::isnan(adjusted_radius) || adjusted_radius<=0.1)
    {
      adjusted_radius=0.0;
    }
    ofs<<"\tmpatches.Circle(("<<pCell->position[0]<<", "<<pCell->position[1]<<"), radius="<<adjusted_radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    for(int j=0; j<neighbors.size();j++)
    {
      Cell* nCell=neighbors[j];
      double neighbor_cell_radius=nCell->phenotype.geometry.radius;//intersection_of_cell_and_plane(z_height, nCell->phenotype.geometry.radius);
      if(std::isnan(neighbor_cell_radius) || neighbor_cell_radius<=0.1)
      {
        neighbor_cell_radius=0.0;
      }
      ofs<<"\tmpatches.Circle(("<<nCell->position[0]<<", "<<nCell->position[1]<<"), radius="<<neighbor_cell_radius<<",alpha=1, ec=\"blue\", fc=\'blue\'),\n";
    }
    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  }
  return;
}
