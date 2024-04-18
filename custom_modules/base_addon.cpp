#include "base_addon.h"


  //---- generic base class for addon
void Base_Addon_Class::initialize(Cell* pCell)  {
  if(pCell->is_active && !pCell->is_out_of_domain)
  {
    this->is_updated=false;
    this->m_daughter=nullptr;
    this->m_pCell=pCell;
    this->pCell_is_safe=true;
  }
  else
  {
   throw std::invalid_argument("assigning Addon to invalid pCell"); 
  }
  return;
}

Base_Addon_Class::~Base_Addon_Class(){
  return;

}

