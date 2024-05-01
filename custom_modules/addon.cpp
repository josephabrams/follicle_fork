#include "./addon.h"
#include "addon_factory.h"
#include "base_addon.h"
// #include "debug_log.h"
std::vector <Addon*> Addon_list{};
// std::vector <Addon*>& Addon::Addon_list_ref=Addon_list;

//----- Addon class for managing a single addon class type
Addon::Addon(Addon_Factory* custom_class_type):m_addon_factory{custom_class_type}{
  return;
}

Addon::~Addon(){
  std::cout<< "Addon Destructor Called!\n";
  //delete all the instances located in map 
  // this could be very slow let the stack handle Addon instances at the end of simulations
  for (int i=0; i<class_instances_by_pCell.size();i++)
  {  
    Base_Addon_Class* destroy_ptr=class_instances_by_pCell[i];
    if(!(destroy_ptr->is_blank)){
      delete class_instances_by_pCell[i];
    }
  }
  for (int i=0; i<detached_instances_by_pCell.size();i++)
  {  
    Base_Addon_Class* destroy_ptr=detached_instances_by_pCell[i];
    if(!(destroy_ptr->is_blank)){
      delete detached_instances_by_pCell[i];
    }
  }
  delete m_blank_ptr; 
  delete m_addon_factory;// cleanup addon factory 1 factory per Addon
 return; 
}

void Addon::spawn_instance(Cell* pCell){
  Base_Addon_Class* ptr= m_addon_factory->Create_Addon_Instance(pCell);
  ptr->is_blank=false;
  class_instances_by_pCell[pCell->index]=ptr;
  check_pCell_safety(pCell);
  return;
}

Base_Addon_Class* Addon::get_instance(Cell* pCell){//could overload [] if code gets really cludgey
  return class_instances_by_pCell[pCell->index];
  
}
void Addon::copy_instance_to_daughter(Base_Addon_Class* instance){
  this->spawn_instance(instance->m_daughter);
  this->class_instances_by_pCell[instance->m_pCell->index]->on_division();
  this->class_instances_by_pCell[instance->m_daughter->index]->on_division();
  instance->m_daughter=nullptr; //reset 
  return;
}

void Addon::detach_instance(Cell* pCell){
  #pragma omp critical
  {
    Base_Addon_Class* detach_ptr=class_instances_by_pCell[pCell->index];
    Base_Addon_Class* blank_ptr= detached_instances_by_pCell[pCell->index];
    if(!(detach_ptr->is_blank))
    {
      detached_instances_by_pCell[pCell->index]=detach_ptr;
      class_instances_by_pCell[pCell->index]->is_blank=true;
    }
  }
  return;
}
void Addon::reattach_instance(Cell* pCell){
  #pragma omp critical
  {
    Base_Addon_Class* attach_ptr=class_instances_by_pCell[pCell->index];
    Base_Addon_Class* detach_ptr= detached_instances_by_pCell[pCell->index];
    if((attach_ptr->is_blank) && !detach_ptr->is_blank)
    {
      detach_ptr->is_blank=true;
      attach_ptr->is_blank=false;
    }
  }
  return;
}

void Addon::check_pCell_safety(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  if(!(instance->is_blank)){
  if(!(pCell->is_active) || pCell->is_out_of_domain){
    
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
  }
  else if(pCell->phenotype.death.dead){
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
    //do dead stuff if you've overwritten what PhysiCell will do by default
  }
  else if(pCell->phenotype.flagged_for_removal){
    // #pragma critical
    // {
      instance->pCell_is_safe=false;
      detach_instance(pCell);
    // }
  }
  else if(pCell->phenotype.flagged_for_division){
     // #pragma critical
    // {
      instance->pCell_is_safe=false;
      // detach_instance(pCell);// maybe you want to detach, maybe not,
      //instance->m_daughter should be set here need some way to get child if you want to attach new addon to it
      // if (instance->m_daughter!=nullptr)
      // {
        // copy_instance_to_daughter(instance);// be sure addon has on_division specified if you want this to work correctly 
      // }
    // }
  } 
  else{ instance->pCell_is_safe=true;}
    instance->is_updated=false;
  }
  return;
}

void Addon::update_custom_class(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  // std::cout<<"Instance: "<< instance<<"\n";
  // std::cout<<"Safety: "<< instance->pCell_is_safe<<"\n";
  if(!(instance->is_updated) && instance->pCell_is_safe && !(instance->is_blank)){
    instance->update_state();
  }
  check_pCell_safety(instance->m_pCell);
  return;
}
Addon* create_Addon(Addon_Factory* custom_class_type ){//do not run in parrallel section
    set_debug();
    Addon* ptr=new Addon(custom_class_type);
  #pragma omp critical
  {
    Base_Addon_Class* null_instance=ptr->m_addon_factory->Create_blank();
    std::cout<< "all cells size: "<< (*all_cells).size()<<"\n";
    ptr->class_instances_by_pCell.resize((*all_cells).size()*2,null_instance);
    ptr->detached_instances_by_pCell.resize((*all_cells).size()*2,null_instance);
    ptr->m_blank_ptr=null_instance;
    std::cout<< "sizes: "<< ptr->class_instances_by_pCell.size()<< "\n";
    Addon_list.push_back(ptr);
  }
    return ptr; 
}
Addon* create_Addon(Addon_Factory* custom_class_type, int final_cell_count ){//do not run in parrallel section
    set_debug();
    Addon* ptr=new Addon(custom_class_type);
  #pragma omp critical
  {
    Base_Addon_Class* null_instance=ptr->m_addon_factory->Create_blank();
    std::cout<< "all cells size: "<< final_cell_count<<"\n";
    ptr->class_instances_by_pCell.resize(final_cell_count,null_instance);
    ptr->detached_instances_by_pCell.resize(final_cell_count,null_instance);
    ptr->m_blank_ptr=null_instance;
    std::cout<< "sizes: "<< ptr->class_instances_by_pCell.size()<< "\n";
    Addon_list.push_back(ptr);
  }
    return ptr; 
}

void clean_up_Addons(){
  // std::cout<<"Cleanup Called! on " << Addon::Addon_list.size()<<" objects \n";
  for(int i=0; i<Addon_list.size();i++){
    Addon* ptr=Addon_list[i];
    delete ptr;
  }
  Addon_list.clear();
}
  

