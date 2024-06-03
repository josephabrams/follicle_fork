#include "twop.h"
#include "../../core/PhysiCell_cell.h"
using namespace PhysiCell;
using namespace BioFVM;

void get_external_concentration(Cell *pCell, std::vector <double> *average_concentrations, std::vector<int> &solute_indecies, std::vector<int> &voxels ){
  std::vector <double> concentration(solute_indecies.size(),0.0);
  for(int i=0; i<solute_indecies.size(); i++){
    double sum=0.0;
    for(int j=0; j<voxels.size(); j++){
      if(microenvironment.nearest_density_vector(voxels[j])[solute_indecies[i]]<1e-12)
      {
        sum+=0.0;
      }
      else {
        sum+=microenvironment.nearest_density_vector( voxels[j])[solute_indecies[i]];
      }
    }
    #pragma omp critical
    {
      concentration[i]=sum/voxels.size();
    }
  }
  #pragma omp critical
  {
   average_concentrations->assign(concentration.begin(),concentration.end()); 
  }
  return;
}


// int find_density_index( std::string name ); 

// void update_exterior_concentrations(Spring_Cell* SPcell)
// {
//   // this function passes the exterior concentrations to the Spring_Cell and specifies what solutes are being simulated
//   // std::cout<<"Simulation selected"<< SPcell->simulation_selected<<"\n";
//   if(SPcell->simulation_selected==1)
//   {
//     // std::cout<<"Salt molarity: "<< concentration_at_boundary(SPcell->m_my_pCell,0)<<"\n";
//     SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);// read in the averaged exterior molarity in the voxels along the cell boundary
//     SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); //convert the molarity into molality using the best fit polynomial
//     
//     // std::cout<<"Eg molarity: "<< concentration_at_boundary(SPcell->m_my_pCell,1)<<"\n";
//     SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
//     
//     // std::cout<<"Test molal: "<< molarity_to_molality((SPcell->exterior_molarity[1]),"EG")<<"\n";
//     SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"EG");
//     
//    // std::cout<<"Test virial: "<< ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG")<<"\n";
//     SPcell->exterior_osmolality=ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG");
//     // SPcell->exterior_component_molality[0]*1.68+SPcell->exterior_component_molality[1];//ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG");
//   }
//   else if(SPcell->simulation_selected==2)
//   { 
//     SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
//     SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
//     SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
//     SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"GLY"); 
//     SPcell->exterior_osmolality=ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","GLY");
//   }
//   else if(SPcell->simulation_selected==3)
//   { 
//     SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
//     SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
//     SPcell->exterior_osmolality=binary_virial(SPcell->exterior_component_molality[0],"NaCl");
//   }
//   else if(SPcell->simulation_selected==4)
//   { 
//     SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
//     SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
//     SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
//     SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"EG"); 
//     SPcell->exterior_molarity[2]=concentration_at_boundary(SPcell->m_my_pCell,2);
//     SPcell->exterior_component_molality[2]=molarity_to_molality(SPcell->exterior_molarity[2],"GLY");
//     SPcell->exterior_osmolality=binary_virial(SPcell->exterior_component_molality[0],"NaCl")+ternary_virial(SPcell->exterior_component_molality[1],SPcell->exterior_component_molality[2],"EG","GLY");
//   }
//   return;
// }
//
//
// void uptake(Spring_Cell* SPcell)
// {
//   // std::ofstream ofs;
//   // ofs.open ("./output/uptake_function.csv", std::ofstream::out | std::ofstream::app);
//   
//   double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
//   double uptake_voxel_num=SPcell->uptake_voxels.size();
//   double reduce_water_volume=SPcell->water_uptake/uptake_voxel_num;
//   // ofs<<PhysiCell_globals.current_time<<", "<<voxel_volume<<", "<<uptake_voxel_num<<", "<<reduce_water_volume<<", ";
//   std::vector<double> specific_volumes;
//   std::vector<double> solute_uptake_per_voxel;
//   // std::cout<<"Number of voxels: "<< SPcell->uptake_voxels.size()<<"\n"; 
//   // std::cout<<"Uptakes: "<<SPcell->solute_uptake<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
//   // std::cout<<"reduce_water_volume: "<<reduce_water_volume<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
//   // ofs<<SPcell->index<<", "<<SPcell->solute_uptake;
//   solute_uptake_per_voxel.resize(SPcell->solute_uptake.size(),0.0);
//   specific_volumes.resize(SPcell->solute_uptake.size(),0.0);
//   for (size_t i =0; i<SPcell->solute_moles.size();i++) {
//     solute_uptake_per_voxel[i]=SPcell->solute_uptake[i]/(double)uptake_voxel_num;
//     if(abs(solute_uptake_per_voxel[i])<1e-12)//deal with very tiny values
//     {
//       solute_uptake_per_voxel[i]=0;
//     }
//     // ofs<< solute_uptake_per_voxel[i]<<", ";
//   }
//   for(size_t i=0; i<SPcell->solute_moles.size();i++){
//     specific_volumes[i]=SPcell->solute_specific_volume[i];
//
//     // ofs<< specific_volumes[i]<<", ";
//   } 
//   // std::cout<<"uptake per voxel: "<<solute_uptake_per_voxel<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
//   //reduce voxel concentration by cell uptake
//   for (size_t i =0; i<SPcell->uptake_voxels.size();i++) {
//     uptake_in_one_voxel(SPcell->uptake_voxels[i],reduce_water_volume,&solute_uptake_per_voxel, SPcell->solute_specific_volume);
//   }
//   // ofs<<"\n";
//   // ofs.close();
//   return;
// }
//
//
//
//
// void uptake_in_one_voxel(int voxel, double water_uptake_per_voxel, std::vector<double>* solute_uptake_per_voxel, std::vector<double> specific_volumes )
// {
//   //current version of code assumes "infinite water bath" based on the phenomological 2P model
//   double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
//   std::vector<double> temp_density_vec=microenvironment.density_vector(voxel);
//   std::vector <double> moles_in_voxel;
//   std::vector <double> debug_moles_output;
//   moles_in_voxel.resize((*solute_uptake_per_voxel).size(),0.0);
//   debug_moles_output.resize((*solute_uptake_per_voxel).size(),0.0);
//   for(int i=0; i<(*solute_uptake_per_voxel).size(); i++)
//   {
//     moles_in_voxel[i]=(temp_density_vec[i]*voxel_volume);//fmole/um^3* um^3 
//   }
//   // std::ofstream ofs;
//   // ofs.open("./output/uptake_in_one_voxel.csv", std::ofstream::out|std::ofstream::app);
//   // ofs <<PhysiCell_globals.current_time<<", "<<microenvironment.voxels(voxel).mesh_index<<", "<<microenvironment.voxels(voxel).center[0]<<", "<<microenvironment.voxels(voxel).center[1]<<", "<<microenvironment.voxels(voxel).center[2]<<", "<<microenvironment.density_vector(voxel)[1]<<", ";
//    
//   std::vector <double> new_moles=moles_in_voxel-(*solute_uptake_per_voxel);
//   double effective_water_from_diffusion=11*voxel_volume-water_uptake_per_voxel; //~water from each voxel in my neighborhood (including diagnols)+ my own
//   double new_water_volume=voxel_volume-(water_uptake_per_voxel/11);//represent approximate change in bulk water which has much higher diffusivity than solutes (at least double)
//   if (effective_water_from_diffusion<0)
//   {
//     std::cout<<"WARNING!! HIGH WATER UPTAKE!!\n";
//     //new_water_volume=new_water_volume+effective_water_from_diffusion;
//   }
//   debug_moles_output=new_moles;
//   for(int i=0; i<new_moles.size(); i++)
//   {
//     if(abs(new_moles[i])<1e-6)
//     {
//       new_moles[i]=0;
//     }
//   }
//   for(int i=0; i<new_moles.size();i++)
//   {
//     if((new_moles[i]/new_water_volume)<0) //cannot take moles that aren't there
//     {
//       std::cout<<"NEGATIVE MOLES! "<<moles_in_voxel[i]<<", "<<(*solute_uptake_per_voxel)[i] <<", "<<i<<"\n";
//       // new_moles[i]=0;
//     }
//   }
//   if(new_water_volume<=0)
//   {
//     std::cout<<"WATER UPTAKE WENT NEGATIVE!!\n";
//     // new_water_volume=0.0000001;
//   }
//   // ofs<< new_moles[1]<<", "<<moles_in_voxel[1]<<", "<<solute_uptake_per_voxel[1]<<", "<<new_water_volume<<", "<< voxel_volume<<", "<<water_uptake_per_voxel<<", "; 
//   // std::ofstream ofs;
//   // ofs.open("./output/uptake_in_one_voxel.csv", std::ofstream::out|std::ofstream::app);
//   // ofs<<PhysiCell_globals.current_time<<", "<<voxel<<", "<<moles_in_voxel<< microenvironment.density_vector(voxel);
//   // std::cout<<"new water volume: "<< new_water_volume<<"\n";
//   for(int i=0; i<(*solute_uptake_per_voxel).size();i++)
//   {
//     // ofs<<"old density: "<< microenvironment.density_vector(voxel)[i]<<",";
//     // std::cout<<"HOLEY MOLEY: "<< new_moles[i]<<"\n";
//     // std::cout<<"old old: "<<microenvironment.density_vector(voxel)[i]<<"\n";
//     // std::cout<<"old: "<<temp_density_vec[i]<<"\n";
//     microenvironment.density_vector(voxel)[i]=new_moles[i]/new_water_volume;
//     // std::cout<<"new: "<<microenvironment.density_vector(voxel)[i]<<"\n";
//   }
//   // ofs<<microenvironment.density_vector(voxel)[1]<<"\n ";
//   // ofs<<microenvironment.density_vector(voxel)<<new_moles<<new_water_volume <<(*solute_uptake_per_voxel)<< "\n";
//   // ofs.close();
//
//   return;
// }
// //
