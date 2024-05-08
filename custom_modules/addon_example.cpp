
#include "./addon_example.h"
#include <fstream>
// #include <cstddef>
Example_Addon_Class::Example_Addon_Class(Cell* pCell, Addon* wrapper){
  this->initialize(pCell, wrapper);//requires initialization of Base_Addon_Class
  this->example_int = 0;
  this->is_blank=false;

  // std::ofstream ofs ("test.txt", std::ofstream::app);

  // ofs <<"Comparisons \n";

  // ofs.close();
  
  return;
}
Example_Addon_Class::Example_Addon_Class(){
  this->initialize_blank();//requires initialization of Base_Addon_Class
  
  return;
}


void Example_Addon_Class::update_state() {
    this->is_updated=true;
    this->example_int=10;
  }

void Example_Addon_Class::on_division() {
    this->example_int = this->example_int/2;
    return;
}

Base_Addon_Class* Example_Addon_Creator::Factory_Method(Cell* pCell, Addon* wrapper) {
  return new Example_Addon_Class(pCell, wrapper);
}

Base_Addon_Class* Example_Addon_Creator::blank_factory_method() {
  return new Example_Addon_Class();
}



