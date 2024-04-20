
#include "./addon_example.h"
#include <fstream>
Example_Addon_Class::Example_Addon_Class(Cell* pCell){
  this->initialize(pCell);//requires initialization of Base_Addon_Class
  this->example_int = 0;
  this->is_blank=false;
  std::ofstream of2;
  of2.open("class_state.txt",std::ofstream::app);
  of2 << "constructed\n";
  of2.close();
  
  return;
}
Example_Addon_Class::Example_Addon_Class(){
  this->initialize_blank();//requires initialization of Base_Addon_Class
  std::ofstream of2;
  of2.open("class_state.txt",std::ofstream::app);
  of2 << "made a blank\n";
  of2.close();
  
  return;
}


void Example_Addon_Class::update_state() {
    std::ofstream of2;
    of2.open("class_state_2.txt",std::ofstream::app);
    of2<<"example int: "<<this->example_int<<"\n";
    of2.close();
    this->is_updated=true;
    this->example_int=10;
  }

void Example_Addon_Class::on_division() {
    this->example_int = this->example_int/2;
    return;
}

Base_Addon_Class* Example_Addon_Creator::Factory_Method(Cell* pCell) {
  return new Example_Addon_Class(pCell);
}

Base_Addon_Class* Example_Addon_Creator::blank_factory_method() {
  return new Example_Addon_Class();
}



