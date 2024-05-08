#include "./addon_factory.h"
//------ factory method for creating instances
Base_Addon_Class* Addon_Factory::Create_Addon_Instance(Cell* pCell, Addon* wrapper){
  // m_addon= wrapper;
  Base_Addon_Class* ptr = this->Factory_Method(pCell, wrapper);
  return ptr;
}
Base_Addon_Class* Addon_Factory::Create_blank(){
  Base_Addon_Class* ptr = this->blank_factory_method();
  return ptr;
}

Addon_Factory::~Addon_Factory(){
  std::cout<< "Factory Destroyed!"<<"\n";
}

