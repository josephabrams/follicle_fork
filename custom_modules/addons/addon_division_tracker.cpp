
#include "./addon_division_tracker.h"
#include <fstream>
// #include <cstddef>
Division_Tracker::Division_Tracker(Cell* pCell, Addon* wrapper){
  this->initialize(pCell,wrapper);//requires initialization of Base_Addon_Class
  int synchronize_division=  (this->m_pCell->index+(*all_cells).size()+2)*2;
  this->division_value = synchronize_division;
  this->m_pCell->custom_data.add_variable("division_value","unitless",synchronize_division);
  int index=this->m_pCell->custom_data.find_variable_index("division_value");
  this->m_pCell->custom_data.variables[index].conserved_quantity=true;
  this->division_found=false;
  this->is_blank=false;
   
  std::ofstream ofs ("test.txt",std::ofstream::app);

  ofs <<"Division Tracker \n";

  ofs.close();
  
  return;
}
Division_Tracker::Division_Tracker(){
  this->initialize_blank();//requires initialization of Base_Addon_Class
  
  return;
}

void Division_Tracker::find_offspring(int search_value){
  for(int i=0; i<(*all_cells).size(); i++)
  {
    Cell* pCell=(*all_cells)[i];
    if(pCell!= m_pCell && pCell->custom_data["division_value"]==search_value)
    {
      this->m_daughter=pCell;
      this->m_wrapper->spawn_instance(pCell);
      // #pragma omp critical
      // {
        // new_cells.push_back(pCell);
      // }
      return;
    }
    
  }
  return;
}
void Division_Tracker::update_state() {
  if (this->division_value != this->m_pCell->custom_data["division_value"])
  {
    #pragma omp critical
    {
      division_found=true;
      this->find_offspring(this->m_pCell->custom_data["division_value"]);
      std::ofstream ofs;
      ofs.open("test.txt", std::ofstream::app);
      ofs <<this->m_pCell<< " split into  "<< this->m_daughter<< "with index value: "<< this->m_daughter->index<<"\n";
      ofs <<this->m_pCell->index<< ", "<< this->division_value<< ", "<<this->m_pCell->custom_data["division_value"] <<"\n";
      ofs.close();
      //reset
      // this->m_daughter=nullptr; //dont reset this for lineage tracking
      //
      int resynchronize_division=  (this->m_pCell->index+(*all_cells).size()+2)*2;
      this->division_value=resynchronize_division;
      this->m_pCell->custom_data["division_value"]=resynchronize_division;
    }
  }
  else{
    std::ofstream ofs;
    ofs.open("test.txt", std::ofstream::app);
    int index=this->m_pCell->custom_data.find_variable_index("division_value");
    bool test_value=this->m_pCell->custom_data.variables[index].conserved_quantity=true;
    ofs <<"No divisions found: "<< this->m_pCell<< ", "<< this->division_value<< ", "<<this->m_pCell->custom_data["division_value"] <<", "<<test_value<<"\n";
    ofs.close();
  }
  this->is_updated=true;
  return;
  }
void Division_Tracker::on_division() {
  #pragma omp critical
  {
    for(int i=0; i<new_cells.size(); i++)
    {
      Cell* pCell= new_cells[i];
      this->m_wrapper->spawn_instance(new_cells[i]);
    }
    new_cells.clear();
  }
  return;
}

Base_Addon_Class* Division_Tracker_Creator::Factory_Method(Cell* pCell, Addon* wrapper) {
  return new Division_Tracker(pCell,wrapper);
}

Base_Addon_Class* Division_Tracker_Creator::blank_factory_method() {
  return new Division_Tracker();
}



