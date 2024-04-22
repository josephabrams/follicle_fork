#ifndef __ADDON_H__
#define __ADDON_H__
#include <array>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include "./quaternion.h"
#include <stdexcept>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./debug_log.h"
#include "addon_factory.h"
#include "base_addon.h"
#include <string>
#include <vector>
using namespace BioFVM;
using namespace PhysiCell;
constexpr int NUMBER_OF_ADDONS=10;

void set_debug();



class Addon{
private:
  
public:
  Base_Addon_Class* m_blank_ptr;
  // static std::vector <Addon*>& Addon_list_ref;
  // #pragma omp threadprivate(Addon_list_ref)
  Addon(Addon const&) = delete; //cannot copy
  Addon& operator=(Addon const&) = delete; // cannot assign
  Addon_Factory* m_addon_factory;
  std::vector <Base_Addon_Class*> class_instances_by_pCell;//store instance
  std::vector <Base_Addon_Class*> detached_instances_by_pCell;//store detached instance
  Addon(Addon_Factory* custom_class_type);
  //Addon_factory* class_one = new Class_One_Creator();
  ~Addon();
  void spawn_instance(Cell* pCell);//make an instance using the factory and then run the creator
    //m_addon_factory->Create_Addon_Instance
  Base_Addon_Class* get_instance(Cell* pCell);
  void copy_instance_to_daughter(Base_Addon_Class* instance);
  void detach_instance(Cell* pCell);
  void reattach_instance(Cell* pCell);
  void check_pCell_safety(Cell* pCell);// if pCell goes out of domain or is destructed deal with addons, maybe an observer is needed
  void update_custom_class(Cell* pCell); //run the class
};
extern std::vector<Addon*> Addon_list;
Addon* create_Addon(Addon_Factory* custom_class_type);
Addon* create_Addon(Addon_Factory* custom_class_type, int final_cell_count );
void clean_up_Addons();

#endif
