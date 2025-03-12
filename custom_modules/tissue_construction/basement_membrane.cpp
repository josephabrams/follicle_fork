#include "./basement_membrane.h"
#include <algorithm>
#include <cmath>
#include <string>
// #include "../multivoxel/multivoxel_functions.h"
using namespace BioFVM;
using namespace PhysiCell;
//should I make the basement membrane a class?
std::vector<int> basement_membrane_voxels={};
std::vector<Cell*> basement_neighbors={};
std::vector<Cell*> basement_initial_neighbors={};
std::vector<Cell*> outter_neighbors={};

double BASEMENT_K=0.0;
void spherical_bounding_region(std::vector <double> &center_point, double &radius,std::vector <int> *return_bounding_region)
{
  std::vector<int> bounding_box_by_index={};
  std::vector<int> spherical_interior_voxels={};
  std::vector<double> half_dimension={radius,radius,radius};
  double voxel_length=default_microenvironment_options.dx;
  general_voxel_bounding_box(&bounding_box_by_index, center_point, half_dimension,voxel_length, microenvironment.mesh);
  //gets voxels of a bounding box around an arbitrary sphere
  //NOTE!!! currently PhysiCells microenvironment voxels are both mechanics and diffusion
  //bounding box of voxels in the microenvironment
  //for speed get center voxel figure out x,y,z offset and only check local voxels
  std::vector<double> center_voxels_center=microenvironment.nearest_voxel(center_point).center;
  std::vector<double> offset=center_voxels_center-center_point;
  
  for (size_t i = 0; i < (bounding_box_by_index).size(); i++)
  {
      
      std::vector<double> test_voxel_center=microenvironment.mesh.voxels[bounding_box_by_index[i]].center;
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
    
      if(default_microenvironment_options.simulate_2D)
      {
        // #pragma omp critical
        test_corners.resize(4,std::vector<double>(3,0.0));
      }
      get_voxel_corners(test_voxel_center,test_corners);
            
      int sum=0;
      if(norm(test_voxel_center-center_point)<=radius)
      {
        spherical_interior_voxels.push_back(bounding_box_by_index[i]);
      }
  }
  #pragma omp critical
  {
    std::sort(spherical_interior_voxels.begin(),spherical_interior_voxels.end());
    std::unique(spherical_interior_voxels.begin(),spherical_interior_voxels.end());
    return_bounding_region->assign(spherical_interior_voxels.begin(),spherical_interior_voxels.end());
  }
  return;
}
void get_basement_membrane_voxels(std::vector <double> &center_point, double &inner_radius, double &outter_radius, std::vector<int> *basement_voxels )
{
  //make two circular voxel boundarys and subtract the inner from the outter, the remaining make up a spherical basement membrane
  std::vector<int> remaining;
  std::vector<int> inner={};
  std::vector<int> outter={};
  spherical_bounding_region(center_point, inner_radius, &inner);
  spherical_bounding_region(center_point, outter_radius, &outter);
  // std::sort(inner.begin(),inner.end());
  // std::sort(outter.begin(),outter.end());
  std::set_difference(outter.begin(),outter.end(), inner.begin(),inner.end(), std::inserter(remaining, remaining.end()));
    // #pragma omp critical
    // {
  basement_voxels->assign(remaining.begin(), remaining.end());
  // basement_membrane_voxels=(*basement_voxels);
  std::cout<<"BASEMENT MEMBRANE VOXELS: "<< basement_membrane_voxels<<"\n";
  // }
  return;
}

void python_plot_BM(double inner_radius, double outter_radius){
  std::vector<double> center{0.0, 0.0, 0.0};
  python_plot_BM(center, inner_radius, outter_radius, basement_membrane_voxels);
}
void python_plot_BM( std::vector<double> center, double inner_radius, double outter_radius, std::vector<int>& bm_voxels)
{

  // if(PhysiCell_globals.current_time<0)
  // {
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/basement_membrane_voxels.py";
    std::ofstream ofs;
    ofs.open (file, std::ofstream::out | std::ofstream::trunc);
    ofs<<"import numpy as np\n"<<"import matplotlib.pyplot as plt\n"<<"import matplotlib.patches as mpatches\n";
    ofs<<"xy_artists = [\n";
    for(int i=0; i<bm_voxels.size(); i++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      
      std::vector <double> voxel_position=microenvironment.voxels(bm_voxels[i]).center;
      if(std::fabs(voxel_position[2])<default_microenvironment_options.dx)
      {
        std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
        ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'),\n";
      
      }
    }
    double voxel_size= default_microenvironment_options.dx;
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<inner_radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<outter_radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  // }
  return;
}
void python_plot_BM_cells(double inner_radius, double outter_radius){
  std::vector<double> center{0.0, 0.0, 0.0};
  
    double region_middle=(inner_radius+outter_radius)/2;
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/basement_membrane_cells.py";
    std::ofstream ofs;
    ofs.open (file, std::ofstream::out | std::ofstream::trunc);
    ofs<<"import numpy as np\n"<<"import matplotlib.pyplot as plt\n"<<"import matplotlib.patches as mpatches\n";
    ofs<<"xy_artists = [\n";
    for(int i=0; i<basement_membrane_voxels.size(); i++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      
      std::vector <double> voxel_position=microenvironment.voxels(basement_membrane_voxels[i]).center;
      if(std::fabs(voxel_position[2])<default_microenvironment_options.dx)
      {
        std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
        ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'),\n";
      
      }
    }
    double voxel_size= default_microenvironment_options.dx;
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<inner_radius<<",alpha=0.2, ec=\"black\", fc=\'none\', linewidth=2.0 ),\n";
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<outter_radius<<",alpha=0.2, ec=\"black\", fc=\'none\', linewidth=2.0 ),\n";
    
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<region_middle<<",alpha=0.2, ec=\"red\", fc=\'None\', linewidth=5.0),\n";
    for(int j=0; j<all_point_springs.size(); j++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      Point_Spring* PS_ptr=all_point_springs[j];
      std::vector <double> cell_position=PS_ptr->m_me->position;
      if(std::fabs(cell_position[2])<default_microenvironment_options.dx)
      {
        ofs<<"\tmpatches.Circle(("<<cell_position[0]<<", "<<cell_position[1]<<"), radius="<<PS_ptr->m_me->phenotype.geometry.radius<<",alpha=0.9, ec=\"black\", fc=\'black\'),\n";
      
      }
    }



    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  // }
  return;




}
void python_plot_BM_all_cells(double inner_radius, double outter_radius){
  std::vector<double> center{0.0, 0.0, 0.0};
  
    double region_middle=(inner_radius+outter_radius)/2;
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/basement_membrane_all_cells"+std::to_string(PhysiCell_globals.current_time)+".py";
    std::ofstream ofs;
    ofs.open (file, std::ofstream::out | std::ofstream::trunc);
    ofs<<"import numpy as np\n"<<"import matplotlib.pyplot as plt\n"<<"import matplotlib.patches as mpatches\n";
    ofs<<"xy_artists = [\n";
    for(int i=0; i<basement_membrane_voxels.size(); i++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      
      std::vector <double> voxel_position=microenvironment.voxels(basement_membrane_voxels[i]).center;
      if(std::fabs(voxel_position[2])<default_microenvironment_options.dx)
      {
        std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
        ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'),\n";
      
      }
    }
    double voxel_size= default_microenvironment_options.dx;
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<inner_radius<<",alpha=0.2, ec=\"black\", fc=\'none\', linewidth=2.0 ),\n";
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<outter_radius<<",alpha=0.2, ec=\"black\", fc=\'none\', linewidth=2.0 ),\n";
    
    ofs<<"\tmpatches.Circle(("<<center[0]<<", "<<center[1]<<"), radius="<<region_middle<<",alpha=0.2, ec=\"red\", fc=\'None\', linewidth=5.0),\n";
    for(int j=0; j<all_point_springs.size(); j++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      Point_Spring* PS_ptr=all_point_springs[j];
      std::vector <double> cell_position=PS_ptr->m_me->position;
      if(std::fabs(cell_position[2])<default_microenvironment_options.dx)
      {
        ofs<<"\tmpatches.Circle(("<<cell_position[0]<<", "<<cell_position[1]<<"), radius="<<PS_ptr->m_me->phenotype.geometry.radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
      
      }
    }

    for(int k=0; k<basement_initial_neighbors.size(); k++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      Cell* pCell=basement_initial_neighbors[k];
      std::vector <double> cell_position=pCell->position;
      if(std::fabs(cell_position[2])<default_microenvironment_options.dx)
      {
        ofs<<"# index: "<<pCell->index<<"\n";
        ofs<<"\tmpatches.Circle(("<<cell_position[0]<<", "<<cell_position[1]<<"), radius="<<pCell->phenotype.geometry.radius<<",alpha=0.2, ec=\"blue\", fc=\'blue\'),\n";
      
      }
    }

    for(int m=0; m<basement_neighbors.size(); m++)
    {
      // std::cout<<"bounding_box size: "<< bounding_box_by_index.size()<<"\n\n";
      Cell* pCell=basement_neighbors[m];
      std::vector <double> cell_position=pCell->position;
      if(std::fabs(cell_position[2])<default_microenvironment_options.dx)
      {
        ofs<<"# index: "<<pCell->index<<"\n";
        ofs<<"\tmpatches.Circle(("<<cell_position[0]<<", "<<cell_position[1]<<"), radius="<<pCell->phenotype.geometry.radius<<",alpha=0.9, ec=\"purple\", fc=\'purple\'),\n";
         
      }
    }


    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  // }
  return;




}
//
void find_multivoxel_neighbors_region(std::vector<int> &voxel_region, std::vector<Cell*> *return_neighbors, double outter_radius, double inner_radius){
  
  std::vector<Cell *> agents_in_voxel={};
  double region_middle=(inner_radius+outter_radius)/2;
  // std::cout<<"microenvironment mesh voxels size: "<< pCell->get_microenvironment()->mesh.voxels.size()<<"\n"; 
  // std::cout<<"microenvironment connected_voxel_indices size: "<< pCell->get_microenvironment()->mesh.connected_voxel_indices.size()<<"\n"; 
  // std::cout<<"microenvironment connected_voxel_indices size: "<< pCell->get_microenvironment()->cartesian_indices(i)<<"\n"; 
  for(int i=0; i<voxel_region.size(); i++)
  {
    //convert to agent grid voxel instead of microenvironment voxel
  
    int grid_voxel= (*all_cells)[0]->get_container()->underlying_mesh.nearest_voxel_index(microenvironment.voxels(voxel_region[i]).center);
    // std::cout<<"GENARAL BOX VOXEL: "<< general_box[i]<<"\n";
    // std::cout<<"AGENT GRID VOXEL: "<< grid_voxel <<"\n"; 
    //   std::cout<< "microenvironment position: "<< microenvironment.voxels(general_box[i]).center<<"\n";
    //   std::cout<< "agent_grid position: "<< pCell->get_container()->underlying_mesh.voxels[grid_voxel].center<<"\n";
    //   std::cout<< "agent_grid size: "<< pCell->get_container()->agent_grid[grid_voxel].size()<<"\n";
  // // #pragma omp private(agents_in_voxel)
      for (int j = 0; j < (*all_cells)[0]->get_container()->agent_grid[grid_voxel].size(); j++)//agent grid holds the cells in each cartesian voxel 
      {
        Cell* temp_agent_ptr=(*all_cells)[0]->get_container()->agent_grid[grid_voxel][j];
        // if not pCell and close enough to be a neighbor 
        // if ( temp_agent_ptr!= pCell )
        double position_length= norm(temp_agent_ptr->position);
        if (position_length+temp_agent_ptr->phenotype.geometry.radius>inner_radius && position_length-temp_agent_ptr->phenotype.geometry.radius<= outter_radius) 
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
  return_neighbors->assign(agents_in_voxel.begin(), agents_in_voxel.end());

  // }
  return;
}
void get_initial_BM_neighbors(double outter_radius, double inner_radius)
{
  find_multivoxel_neighbors_region( basement_membrane_voxels, &basement_initial_neighbors, outter_radius, inner_radius);
  return;
}
void create_BM_springs(double outter_radius, double inner_radius)
{
  for(int i=0; i<basement_initial_neighbors.size();i++)
  {
    Cell* BM_neighbor=basement_initial_neighbors[i];  
    double rest_length=((inner_radius+outter_radius)/2) - (norm(BM_neighbor->position)+BM_neighbor->phenotype.geometry.radius);
    double spring_constant=BASEMENT_K;//BM_neighbor->custom_data["spring_k"];
    std::cout<<"SPRING CONSTANT: "<< spring_constant<<"\n";
    create_point_spring(BM_neighbor, rest_length, spring_constant); 
  }
  return;
}
void update_BM_neighbors(double outter_radius, double inner_radius)
{

  std::vector<Cell*> remaining_neighbors={};
  std::vector<Cell*> current_neighbors={};
  find_multivoxel_neighbors_region(basement_membrane_voxels, &current_neighbors, outter_radius, inner_radius);
  std::set_difference(current_neighbors.begin(),current_neighbors.end(), basement_initial_neighbors.begin(),basement_initial_neighbors.end(), std::inserter(remaining_neighbors, remaining_neighbors.end()));
  basement_neighbors=remaining_neighbors;
  //
  // std::string file= "./output/update_BM_neighbors"+std::to_string(PhysiCell_globals.current_time)+".csv";
  // std::ofstream ofs;
  // ofs.open (file, std::ofstream::out | std::ofstream::trunc);
  // for(int i=0; i<basement_initial_neighbors.size(); i++)
  // {
  //   ofs<< "initial: "<<basement_initial_neighbors[i]->index<<"\n";
  // }
  // for(int j=0; j<current_neighbors.size(); j++)
  // {
  //   ofs<<"current: "<< current_neighbors[j]->index<<"\n";
  // }
  // for(int k=0; k<basement_neighbors.size(); k++)
  // {
  //   ofs<<"bm: "<< basement_neighbors[k]->index<<"\n";
  // }
  // ofs.close();
  // std::cout<<"initial_neighbors_size: "<< basement_initial_neighbors.size()<<"\n";

  // std::cout<<"basement_neighbors size: "<< basement_neighbors.size()<<"\n";
}
//not thread safe to run in parallel - could be made so by updating all spring_lengths seperately
void advance_BM_springs(double outter_radius, double inner_radius)
{
  for(int i=0; i<all_point_springs.size(); i++)
  {
    Point_Spring* PS_ptr= all_point_springs[i];

    double spring_length=((inner_radius+outter_radius)/2) - (norm(PS_ptr->m_me->position)+(PS_ptr->m_me->phenotype.geometry.radius));
    PS_ptr->m_force_normal=(1/norm(PS_ptr->m_me->position))*PS_ptr->m_me->position;
    PS_ptr->m_spring_length=spring_length;
    PS_ptr->calculate_spring_force();
    double bm_midpoint=((inner_radius+outter_radius)/2);
    
    if(norm(PS_ptr->m_me->position)+PS_ptr->m_me->phenotype.geometry.radius>= bm_midpoint )
    {
      PS_ptr->m_force=-1.0*PS_ptr->m_force;
    }
    Cryocell* cMe=static_cast<Cryocell*>(PS_ptr->m_me);
    PS_ptr->update_force_vector(&cMe->net_force);
  }
  //if a cell overlaps the midpoint exert force inward, if a cell passes beyond the BM it will no longer be a neighbor
  for(int j=0; j<basement_neighbors.size(); j++)
  {
    // std::cout<< "WARNING!!NON-INITIAL NEIGHBOR MADE CONTACT WITH THE BM!\n";
    double bm_midpoint=((inner_radius+outter_radius)/2);
    Cell* pCell=basement_neighbors[j];
    Cryocell* cCell=static_cast<Cryocell*>(pCell);
    if(norm(cCell->position)+cCell->phenotype.geometry.radius>= bm_midpoint )
    {
      std::cout<< "BM CROSSED!!!\n";
      double spring_length= ((inner_radius+outter_radius)/2) - (norm(cCell->position)+cCell->phenotype.geometry.radius);
      double delta_x=std::fabs(spring_length); //force points from me to neighbor
      std::vector<double> unit_vec=(1/norm(cCell->position))*cCell->position;
      double spring_constant=BASEMENT_K;//parameters.doubles("bm_spring_constant");//pCell->custom_data["spring_k"];
      std::vector<double> force=(1.0*spring_constant*(delta_x))*unit_vec;
      cCell->net_force= cCell->net_force+force;
    }

  }
}

void linker_test(){
  std::cout<<"linker worked!"<<"\n";
}
