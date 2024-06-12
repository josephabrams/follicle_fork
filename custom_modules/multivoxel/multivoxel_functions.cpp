
#include "./multivoxel_functions.h"
using namespace BioFVM;
using namespace PhysiCell;
Voxel ghost_voxel;
std::vector<Cell*> multivoxel_cells;
void general_voxel_bounding_box(std::vector<int> *return_bounding_box, std::vector<double> center, std::vector<double> half_dimensions, double voxel_length, BioFVM::Cartesian_Mesh &a_mesh)  
{
  
 //  double Xmin = microenvironment.mesh.bounding_box[0]; 
	// double Ymin = microenvironment.mesh.bounding_box[1]; 
	// double Zmin = microenvironment.mesh.bounding_box[2]; 
	//
	// double Xmax = microenvironment.mesh.bounding_box[3]; 
	// double Ymax = microenvironment.mesh.bounding_box[4]; 
	// double Zmax = microenvironment.mesh.bounding_box[5]; 
  std::vector<int> bounding_box_by_index;
  std::vector<double> voxel_dimensions(3,voxel_length);
  std::vector<double> lower_point{0.0,0.0,0.0};
  std::vector<double> upper_point{0.0,0.0,0.0};
  lower_point=center-(half_dimensions)-voxel_dimensions;
  upper_point=center+(half_dimensions)+voxel_dimensions;
  
    for(int i=0; i<3; i++){
      if(a_mesh.bounding_box[i]>lower_point[i]){
        lower_point[i]=a_mesh.bounding_box[i]+(0.5*voxel_length);
      }

      if(a_mesh.bounding_box[i+3]<upper_point[i]){
        upper_point[i]=a_mesh.bounding_box[i+3]-(0.5*voxel_length);
      }
    }
  // std::vector<double> voxel_start= a_mesh.nearest_voxel(lower_point).center;
  // std::vector<double> voxel_end= a_mesh.nearest_voxel(upper_point).center;
  std::vector<double> voxel_start= microenvironment.nearest_voxel(lower_point).center;
  std::vector<double> voxel_end= microenvironment.nearest_voxel(upper_point).center;
  if(voxel_start==voxel_end)
  {
    #pragma omp critical
    {

      return_bounding_box->assign(voxel_start.begin(),voxel_end.end());
      // return_bounding_box->resize(1);
      // return_bounding_box->at(0)=a_mesh.nearest_voxel_index(voxel_start);
      // std::cout<<"single box: "<< return_bounding_box;
    }
    //contained entirely in 1 voxel my exterior is my voxel
    return;
  }
  else
  {
    // #pragma omp parallel for collapse(3)
    // {
      // std::vector<int> bounding_box_by_index_private;
      // #pragma omp for nowait
      for (double x=voxel_start[0];x<voxel_end[0]; x+=voxel_length){
        for (double y=voxel_start[1];y<voxel_end[1]; y+=voxel_length){
          if(default_microenvironment_options.simulate_2D){
              // std::cout<<"("<<x<<", "<<y<<")"<<"\n";
              std::vector<double> voxel_position {x,y,0};
              bounding_box_by_index.push_back(a_mesh.nearest_voxel_index(voxel_position)); 
          }
          else{
            for (double z=voxel_start[2];z<voxel_end[2]; z+=voxel_length){
              // std::cout<<"("<<x<<", "<<y<<", "<<z<< ")"<<"\n";
              std::vector<double> voxel_position {x,y,z};
              bounding_box_by_index.push_back(a_mesh.nearest_voxel_index(voxel_position)); 
            }
          }
        }
      }
      #pragma omp critical
      {

        // std::vector<int>::iterator it;
        // it = std::unique (bounding_box_by_index.begin(), bounding_box_by_index.end());   // 10 20 30 20 10 ?  ?  ?  ?
        // bounding_box_by_index.resize( std::distance(bounding_box_by_index.begin(),it) );
      //std::cout<<"bounding_box_by_index : "<<bounding_box_by_index<<"\n\n";
        // bounding_box_by_index.insert(bounding_box_by_index.end(), bounding_box_by_index_private.begin(), bounding_box_by_index_private.end());
        return_bounding_box->assign(bounding_box_by_index.begin(),bounding_box_by_index.end());
        // std::cout<<"return_bounding_box ("<< (*return_bounding_box).size()<<"): "<<*return_bounding_box<<"\n";
        if((*return_bounding_box).size()==0){
          // std::cout<<"return_bounding_box ("<< (*return_bounding_box).size()<<"): "<<*return_bounding_box<<"\n\n";
          // if(voxel_end[0]<voxel_start[0] || voxel_end[1]<voxel_start[1] || voxel_end[2]<voxel_end[2]){
            // std::cout<<"center "<< center <<"\n";
            // std::cout<<"start voxel: "<< voxel_start <<"\n";
            // std::cout<<"end voxel: "<< voxel_end <<"\n\n";
            // std::cout<<"lower_point: "<< lower_point <<"\n";
            // std::cout<<"upper_point: "<< upper_point <<"\n\n\n\n";
          // }
        }
      }
    return;
  }
}



void get_voxel_corners(std::vector<double> &voxel_center, std::vector<std::vector<double>> &return_corners )
{
  //find all 8 corners of the voxel cube
  double zz=default_microenvironment_options.dz/2.0;
  double yy=default_microenvironment_options.dy/2.0;
  double xx=default_microenvironment_options.dx/2.0;
  std::vector <std::vector <double>> corners(4,std::vector<double>(3,0.0));
  int count=0;
  for (int i = -1; i < 2; i+=2)
  {
    for (int j = -1; j < 2; j+=2)
    {
      if(!default_microenvironment_options.simulate_2D){
        #pragma omp critical
        corners.resize(8,std::vector<double>(3,0.0));
        for (int k = -1; k < 2; k+=2)
        {
          std::vector<double> temp_point={xx*i,yy*j,zz*k};
          corners[count]=voxel_center+temp_point;
          count++;
        }
      }
      else {

          std::vector<double> temp_point={xx*i,yy*j,0};
          corners[count]=voxel_center+temp_point;
          count++;
      }
    }   
  }
  #pragma omp critical
  {
    for(int i =0; i<count; i++)
    {
      return_corners[i].assign(corners[i].begin(),corners[i].end());
    }
    // return_corners[0].assign(corners[0].begin(),corners[0].end());
    // return_corners[1].assign(corners[1].begin(),corners[1].end());
    // return_corners[2].assign(corners[2].begin(),corners[2].end());
    // return_corners[3].assign(corners[3].begin(),corners[3].end());
    // return_corners[4].assign(corners[4].begin(),corners[4].end());
    // return_corners[5].assign(corners[5].begin(),corners[5].end());
    // return_corners[6].assign(corners[6].begin(),corners[6].end());
    // return_corners[7].assign(corners[7].begin(),corners[7].end());
  }
  return;
}

//return edges as sets of x, y and z
void get_voxel_edges(std::vector<double> &voxel_center, std::vector<std::vector<double>> &x_edges,std::vector<std::vector<double>> &y_edges,std::vector<std::vector<double>> &z_edges)
{
  //find all 8 corners of the voxel cube
  double zz=default_microenvironment_options.dz/2.0;
  double yy=default_microenvironment_options.dy/2.0;
  double xx=default_microenvironment_options.dx/2.0;
  std::vector <std::vector <double>> edges(6,std::vector<double>(3,0.0));
  int count=0;
  for (int i = -1; i < 2; i+=2)
  {
      std::vector<double> temp_point={xx*i, 0.0, 0.0};
      edges[count]=voxel_center+temp_point;
      temp_point={0.0, yy*i, 0.0};
      edges[count+2]=voxel_center+temp_point;
      if(default_microenvironment_options.simulate_2D)
      {
        temp_point={0.0, 0.0, zz*i};
        edges[count+4]=voxel_center+temp_point;
      }
      count++;
  }
  #pragma omp critical
  {
    for(int i =0; i<2; i++)
    {
      x_edges[i].assign(edges[i].begin(),edges[i].end());
      y_edges[i].assign(edges[i+2].begin(),edges[i+2].end());
      z_edges[i].assign(edges[i+4].begin(),edges[i+4].end());

    }

  }
  return;
}
bool edge_intersect(Cell* pCell, std::vector <std::vector<double>> &edges, int& face_count){
  bool is_interesect=false;
  #pragma omp critical
  { 
    std::vector <double> inner_edge(3,0.0);
    std::vector <double> outter_edge(3,0.0);

    if(std::abs(norm(pCell->position-edges[0]))<std::abs(norm(pCell->position-edges[1])))
    {
      inner_edge=edges[0];
      outter_edge=edges[1];
    }
    else {
      inner_edge=edges[1];
      outter_edge=edges[0];
    }

    if(norm(inner_edge-pCell->position)<=pCell->phenotype.geometry.radius && norm(outter_edge-pCell->position)>pCell->phenotype.geometry.radius){
      is_interesect=true;
      face_count+=1;
    }
  }
  return is_interesect;
}
bool corner_intersect(Cell* pCell, std::vector<std::vector<double>> &test_corners, int& sum){
  bool is_intersect=false;
  #pragma omp critical
  {
    for(size_t j =0; j<test_corners.size();j++)
    {
      // std::cout<<j<<" "<<test_corners[j]<<" distance: "<< norm(test_corners[j]-pCell->position)<<" pos: "<< pCell->position<<std::endl;
      if(norm(test_corners[j]-pCell->position)<pCell->phenotype.geometry.radius)
      {
        is_intersect=true;
        sum+=1;
      }
        
    }
  }
  return is_intersect;
}

void get_intersecting_voxels(Cell* pCell,std::vector<int>& bounding_voxels,std::vector<int>* return_intersecting_voxel_indicies)
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    // double voxel_length=default_microenvironment_options.dx;
    // std::vector<double> radial_dimensions(3,pCell->phenotype.geometry.radius);
    // general_voxel_bounding_box(return_intersecting_voxel_indicies, pCell->position, radial_dimensions,voxel_length, pCell->get_microenvironment()->mesh);
    // diffusion_bounding_box(pCell,&bounding_voxels);
    // std::vector<int> intersecting_voxels={};
    // std::cout<<"BOUNDING BOX SIZE: "<< bounding_voxels.size()<<"\n";
    // if(bounding_voxels.size()==1)
    // { 
      // (*return_intersecting_voxel_indicies)[0]=bounding_voxels[0];
      // return;
    // }
    std::vector<int> intersecting_voxels{};
    // #pragma omp private(intersecting_voxels);
    for (size_t i = 0; i < (bounding_voxels).size(); i++)
    {
      
      // std::vector<int> intersecting_voxels_private{};
      std::vector<double> test_voxel_center=pCell->get_microenvironment()->mesh.voxels[(bounding_voxels)[i]].center;
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      std::vector<std::vector <double>> test_x_edges(2,std::vector<double>(3,0.0));
      std::vector<std::vector <double>> test_y_edges(2,std::vector<double>(3,0.0));
      std::vector<std::vector <double>> test_z_edges(2,std::vector<double>(3,0.0));

      if(default_microenvironment_options.simulate_2D)
      {
        #pragma omp critical
        test_corners.resize(4,std::vector<double>(3,0.0));
      }
      // std::vector<std::vector <double>> test_faces(6,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
      get_voxel_edges(test_voxel_center,test_x_edges,test_y_edges,test_z_edges);
      int sum=0;
      int face_count=0;
      if(sum==0 && face_count==0&&test_voxel_center==pCell->get_microenvironment()->nearest_voxel(pCell->position).center)
      {
                    intersecting_voxels.push_back(bounding_voxels[i]);
                    sum=1;
      }

      else if(sum==0&&face_count==0&&corner_intersect(pCell,test_corners,sum))
      {
       intersecting_voxels.push_back(bounding_voxels[i]);
      }
      else if(sum==0 && face_count==0 && edge_intersect(pCell, test_x_edges,face_count))
      {        
       intersecting_voxels.push_back(bounding_voxels[i]); 
      }
      else if(sum==0 && face_count==0 && edge_intersect(pCell, test_y_edges,face_count))
      {        
       intersecting_voxels.push_back(bounding_voxels[i]);
      }
      else if(!default_microenvironment_options.simulate_2D && sum==0 &&face_count==0 && edge_intersect(pCell, test_z_edges,face_count)){
        intersecting_voxels.push_back(bounding_voxels[i]);
      }
    }
    #pragma omp critical
    {   
      // std::vector<int>::iterator it;
      // it = std::unique (intersecting_voxels.begin(), intersecting_voxels.end());   // 10 20 30 20 10 ?  ?  ?  ?
      // intersecting_voxels.resize( std::distance(intersecting_voxels.begin(),it) ); 
        
      return_intersecting_voxel_indicies->assign(intersecting_voxels.begin(),intersecting_voxels.end());
    // std::cout<< "returing voxels size: "<< return_intersecting_voxel_indicies->size()<<"\n\n";
    }
    
    return;
}



void get_exterior_voxels(Cell* pCell, std::vector<double>* return_exterior_voxel_indicies)//some but not all the coners of the voxel are within the cell and the voxel center is outside the cell
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels={};

    double voxel_length=default_microenvironment_options.dx;
    std::vector<double> radial_dimensions(3,pCell->phenotype.geometry.radius);
    general_voxel_bounding_box(&bounding_voxels, pCell->position, radial_dimensions,voxel_length, pCell->get_microenvironment()->mesh);
    // diffusion_bounding_box(pCell,&bounding_voxels);
    std::vector<int> exterior_voxels={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
    
     
      int sum=0;
      // #pragma omp private(sum,exterior_voxels)
      for(size_t j =0; j<test_corners.size();j++)
      {
        //std::cout<<j<<" "<<test_corners[j]<<" distance: "<< norm(test_corners[j]-pCell->position)<<std::endl;
        if(norm(test_corners[j]-pCell->position)<pCell->phenotype.geometry.radius)
        {
          sum+=1;
        }
      }
      
      if(sum!=8 && sum !=0 && norm(test_voxel_center-pCell->position)>pCell->phenotype.geometry.radius)
      {
        exterior_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_exterior_voxel_indicies->assign(exterior_voxels.begin(), exterior_voxels.end());
    }

    return;
}



void get_interior_voxels(Cell* pCell, std::vector<int>* return_interior_voxel_indicies) //voxel center is inside the cell
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector<int> bounding_voxels={};

    double voxel_length=default_microenvironment_options.dx;
    std::vector<double> radial_dimensions(3,pCell->phenotype.geometry.radius);
    general_voxel_bounding_box(&bounding_voxels, pCell->position, radial_dimensions,voxel_length, pCell->get_microenvironment()->mesh);
    // diffusion_bounding_box(pCell,&bounding_voxels);
    
    std::vector<int> interior_voxels={};
    // #pragma omp private(interior_voxels)

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
     
      int sum=0;
      // #pragma omp private(interior_voxels)

      if(norm(test_voxel_center-pCell->position)<=pCell->phenotype.geometry.radius)
      {
        interior_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_interior_voxel_indicies->assign(interior_voxels.begin(), interior_voxels.end());
    }

    return;
}


void intersecting_neighbor_voxels(Cell* pCell, Cell* pNeighbor, std::vector<int> my_bounding_voxels, std::vector<int> neighbor_bounding_voxels, std::vector<int> *return_voxels){
  std::vector<int> my_voxels{};
  std::vector<int> neighbor_voxels{};

  get_intersecting_voxels(pCell, my_bounding_voxels, &my_voxels);
  get_intersecting_voxels(pNeighbor, neighbor_bounding_voxels, &neighbor_voxels);
  std::vector<int> intersecting_voxels{};
  // #pragma omp private(intersecting_voxels)
  for(int i=0; i<my_voxels.size();i++){
    for(int j=0; j<neighbor_voxels.size();j++){
      if(my_voxels[i]==neighbor_voxels[j])
      {
        intersecting_voxels.push_back(my_voxels[i]);
      }
    }
  }
  #pragma omp critical
  { 
    // if(intersecting_voxels.size()==1){
      // intersecting_voxels.push_back(0);    }
      return_voxels->assign(intersecting_voxels.begin(), intersecting_voxels.end());
  }
}


void check_out_of_bounds(Cell* pCell, int fail_count)
{
  std::vector<double> position=pCell->position;
  double length=pCell->phenotype.geometry.radius;
  double Xmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[0]; 
	double Ymin = BioFVM::get_default_microenvironment()->mesh.bounding_box[1]; 
	double Zmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[2]; 

	double Xmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[3]; 
	double Ymax = BioFVM::get_default_microenvironment()->mesh.bounding_box[4]; 
	double Zmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
    
    double xs = position[0] - length;
    double xe = position[0] + length;
    double ys = position[1] - length;
    double ye = position[1] + length;
    double zs = 0.0;
    double ze = 0.0;
    if (default_microenvironment_options.simulate_2D) {
    }
    else if (!default_microenvironment_options.simulate_2D) {
        zs = position[2] - length;
        ze = position[2] + length;
    }

    /* check whether a cell leaves the domain and if so initialise fibre again
                assume user placed the centre of fibre within the domain so reinitialise orientation,
                break after 10 failures
                It needs re-writing at some stage to handle the 3D case properly */

    if (PhysiCell::parameters.bools("multivoxel_off_boundary")) {
        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
            fail_count = 10;
        }
    }
    else{
        if (default_microenvironment_options.simulate_2D) {
            while (fail_count < 10) {
                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                    fail_count++;
                    xs = position[0] - length;
                    xe = position[0] + length;
                    ys = position[1] - length;
                    ye = position[1] + length;
                }
                else {
                    break;
                }
            }
        }

        if (!default_microenvironment_options.simulate_2D) {
            while (fail_count < 10) {
                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                    zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                    fail_count++;
                    xs = position[0] - length;
                    xe = position[0] + length;
                    ys = position[1] - length;
                    ye = position[1] + length;
                    zs = position[2] - length;
                    ze = position[2] + length;
                }
                else {
                    break;
                }
            }
        }
    }
}

void python_plot_cell_and_voxels(Cell* pCell, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name)
{
  if(PhysiCell_globals.current_time>1 && PhysiCell_globals.current_time<1.2)
  {
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/cell-"+plot_name+std::to_string(pCell->index)+"-singe_cell_plot.py";
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
      std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
      ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'),\n";
    }
    ofs<<"\tmpatches.Circle(("<<pCell->position[0]<<", "<<pCell->position[1]<<"), radius="<<pCell->phenotype.geometry.radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  }
  return;
}


void python_plot_two_cells_and_voxels(Cell* pCell, Cell* pNeighbor, double dt, std::vector<int> &bounding_box_by_index, std::string plot_name)
{

  if(PhysiCell_globals.current_time>1 && PhysiCell_globals.current_time<1.2)
  /*if(PhysiCell_globals.current_time>dt)*/
  {
    double voxel_length=default_microenvironment_options.dx;
    std::string file= "./output/cell-"+plot_name+std::to_string(pCell->index)+"-singe_cell_plot.py";
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
      std::vector <double> bottom_corner={voxel_position[0]-(voxel_length/2),voxel_position[1]-(voxel_length/2)};      
      ofs<<"\tmpatches.Rectangle(("<<bottom_corner[0]<<", "<<bottom_corner[1]<<"), "<<voxel_length<<","<< voxel_length << ", alpha=0.5, ec=\"red\", fc=\'green\'),\n";
    }
    ofs<<"\tmpatches.Circle(("<<pCell->position[0]<<", "<<pCell->position[1]<<"), radius="<<pCell->phenotype.geometry.radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    // ofs<<"]\n";
    ofs<<"\tmpatches.Circle(("<<pNeighbor->position[0]<<", "<<pNeighbor->position[1]<<"), radius="<<pNeighbor->phenotype.geometry.radius<<",alpha=0.2, ec=\"black\", fc=\'black\'),\n";
    ofs<<"]\n";
    ofs<<"fig,ax=plt.subplots()\n"<<"for i in xy_artists:\n"<<"\tax.add_patch(i)\n"<<"ax.autoscale_view()\n"<<"ax.set_aspect('equal', 'box')\n"<<"plt.show()";
    ofs.close();
  }
  return;
}
