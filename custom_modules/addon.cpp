#include "./addon.h"
#include "debug_log.h"
#include <string>
bool DEBUG_MODE=true;
auto DEBUG_LOG = Log::get_instance();
void set_debug(){
  auto log_file = "./debug_log.txt";
  auto log_mode = Log_Destination::FILE;
  DEBUG_LOG->set_destination(log_file, log_mode);
}
  //---- generic base class for addon
void Base_Addon_Class::initialize(Cell* pCell)  {
  if(pCell->is_active && !pCell->is_out_of_domain)
  {
    this->is_updated=false;
    this->m_daughter=nullptr;
    this->m_pCell=pCell;
    this->pCell_is_safe=true;

    if(DEBUG_MODE){
      DEBUG_LOG->Log_this("Base_Addon_Class initialized!");  
  }

  }
  else
  {
   throw std::invalid_argument("assigning Addon to invalid pCell"); 
  }
  return;
}

Base_Addon_Class::~Base_Addon_Class(){
  return;
  if(DEBUG_MODE){
    DEBUG_LOG->Log_this("Base_Addon_Class Created!");  
  }
}

//------ factory method for creating instances
Base_Addon_Class* Addon_Factory::Create_Addon_Instance(Cell* pCell){
  Base_Addon_Class* ptr = this->Factory_Method(pCell);
  return ptr;
}

//----- Addon class for managing a single addon class type
Addon::Addon(Addon_Factory* custom_class_type):m_addon_factory{custom_class_type}{
  if(DEBUG_MODE){
    DEBUG_LOG->Log_this("Addon Created!");  
  }

  return;
}

Addon::~Addon(){
  //delete all the instances located in map 
  // this could be very slow let the stack handle Addon instances at the end of simulations
  if(DEBUG_MODE){
    DEBUG_LOG->Log_this("Addon Destructed!");  
  }
  
  for ( auto it = class_instances_by_pCell.begin(); it != class_instances_by_pCell.end(); ++it )
  {  delete it->second;}
  for ( auto it = detached_instances_by_pCell.begin(); it != detached_instances_by_pCell.end(); ++it )
  {  delete it->second;}
 return; 
}

void Addon::spawn_instance(Cell* pCell){
  Base_Addon_Class* ptr= m_addon_factory->Create_Addon_Instance(pCell);
  std::pair<int,Base_Addon_Class*> new_instance (pCell->index,ptr);
  class_instances_by_pCell.insert(new_instance);
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
  Base_Addon_Class* detach_ptr=class_instances_by_pCell[pCell->index];
  std::pair<int,Base_Addon_Class*> detached_instance (pCell->index,detach_ptr);
  detached_instances_by_pCell.insert(detached_instance);
  class_instances_by_pCell.erase(pCell->index);
  return;
}

void Addon::check_pCell_safety(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  if(!pCell->is_active || pCell->is_out_of_domain){
    
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
  }
  else if(pCell->phenotype.death.dead){
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
    //do dead stuff if you've overwritten what PhysiCell will do by default
  }
  else if(pCell->phenotype.flagged_for_removal){
    #pragma critical
    {
      instance->pCell_is_safe=false;
      detach_instance(pCell);
    }
  }
  else if(pCell->phenotype.flagged_for_division){
     #pragma critical
    {
      instance->pCell_is_safe=false;
      detach_instance(pCell);// maybe you want to detach, maybe not,
      //instance->m_daughter should be set here need some way to get child if you want to attach new addon to it
      if (instance->m_daughter!=nullptr)
      {
        copy_instance_to_daughter(instance);// be sure addon has on_division specified if you want this to work correctly 
      }
    }
  } 
  else{ instance->pCell_is_safe=true;}
  if(DEBUG_MODE){
    std::string message= "pCell safety checked and found to be: "+ std::to_string(instance->pCell_is_safe);
    DEBUG_LOG->Log_this(message);  
  }

  return;
}

void Addon::update_custom_class(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  if(!(instance->is_updated) || instance->pCell_is_safe){
    instance->update_state();
    if(DEBUG_MODE){
      DEBUG_LOG->Log_this("Custom Class State Updated!");  
    }
  }
  check_pCell_safety(instance->m_pCell);
  return;
}


