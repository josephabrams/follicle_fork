#include "./multivoxel_uptake.h"

void calculate_per_voxel_uptake(Cell* pCell, std::vector<int> &cell_voxels, double &water_volume_uptake, std::vector<double> &solute_moles_uptakes)
{
  double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
  double uptake_voxel_num=cell_voxels.size();
  double reduce_water_volume=water_volume_uptake/uptake_voxel_num;
  std::vector<double> specific_volumes(solute_moles_uptakes.size(),0.0);
  std::vector<double> solute_uptake_per_voxel(solute_moles_uptakes.size(),0.0);
  for (size_t i =0; i<solute_moles_uptakes.size();i++) {
    solute_uptake_per_voxel[i]=solute_moles_uptakes[i]/(double)uptake_voxel_num;
    if(abs(solute_uptake_per_voxel[i])<1e-12)//deal with very tiny values
    {
      solute_uptake_per_voxel[i]=0;
    }
    // ofs<< solute_uptake_per_voxel[i]<<", ";
  }

  //reduce voxel concentration by cell uptake
  for (size_t i =0; i<cell_voxels.size();i++) {
    uptake_in_one_voxel(cell_voxels[i],reduce_water_volume,solute_uptake_per_voxel);
  }
  // ofs<<"\n";
  // ofs.close();
  return;
}



//water uptake in um^3, solute uptake in fmoles
void uptake_in_one_voxel(int &voxel, double& water_uptake_per_voxel, std::vector<double>& solute_uptake_per_voxel)
{
  //calculate moles of each solute in the voxel so we can reduce it by the known change in moles
  double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
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
      new_moles[i]=0;
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
// //
void custom_update_voxels(){
  return;
}
