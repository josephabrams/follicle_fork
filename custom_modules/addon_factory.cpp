#include "./addon_factory.h"
//------ factory method for creating instances
Base_Addon_Class* Addon_Factory::Create_Addon_Instance(Cell* pCell){
  Base_Addon_Class* ptr = this->Factory_Method(pCell);
  return ptr;
}

Addon_Factory::~Addon_Factory(){
  std::cout<< "Factory Destroyed!"<<"\n";
}

