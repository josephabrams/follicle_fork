#include "base_addon.h"


  //---- generic base class for addon
void Base_Addon_Class::initialize(Cell* pCell, Addon* wrapper)  {
  if(pCell->is_active && !pCell->is_out_of_domain)
  {
    this->is_blank=false;
    this->is_updated=false;
    this->m_daughter=nullptr;
    this->m_pCell=pCell;
    this->pCell_is_safe=true;
    this->m_wrapper=wrapper;
  }
  else
  {
   throw std::invalid_argument("assigning Addon to invalid pCell"); 
  }
  return;
}
void Base_Addon_Class::initialize_blank(){
  this->is_blank=true;
  return;
}

Base_Addon_Class::~Base_Addon_Class(){
  return;

}

