#ifndef __BASE_ADDON_H__
#define __BASE_ADDON_H__
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./debug_log.h"
using namespace BioFVM;
using namespace PhysiCell;

class Base_Addon_Class{
private:
public:
  bool is_blank;
  Cell* m_pCell;
  bool pCell_is_safe;
  bool is_updated;
  bool divided;
  Cell* m_daughter;
  void initialize(Cell* pCell);
  void initialize_blank(); //empty place holder to fill class containers
  //Base_Addon_Class(Cell* pCell);
  virtual ~Base_Addon_Class();
  Base_Addon_Class* get_instance();
  virtual void update_state()=0;
  virtual void on_division()=0;
};

#endif // !__BASE_ADDON_H__
