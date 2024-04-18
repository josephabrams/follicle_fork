
#include "./addon_example.h"
Example_Addon_Class::Example_Addon_Class(Cell* pCell){
  this->initialize(pCell);//requires initialization of Base_Addon_Class
  this->example_int = pCell->index*2;
  
  return;
}


void Example_Addon_Class::update_state() {
    this->is_updated=true;
  }

void Example_Addon_Class::on_division() {
    this->example_int = this->example_int/2;
    return;
}

Base_Addon_Class* Example_Addon_Creator::Factory_Method(Cell* pCell) {
  return new Example_Addon_Class(pCell);
}
