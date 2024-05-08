#ifndef __ADDON_DIVISION_TRACKER_H__
#define __ADDON_DIVISION_TRACKER_H__
#include "./addon.h"
#include "base_addon.h"
#include <omp.h>
class Division_Tracker: public Base_Addon_Class{
  public:
    std::vector <Cell*> new_cells;
    Division_Tracker(Cell* pCell, Addon* wrapper);
    Division_Tracker();
    void find_offspring(int search_value);
    int division_value;
    bool division_found;
    void update_state() override;

    void on_division() override;
};


class Division_Tracker_Creator : public Addon_Factory{
  Base_Addon_Class* Factory_Method(Cell* pCell, Addon* wrapper) override;
  Base_Addon_Class* blank_factory_method() override;
};
#endif //__ADDON_EXAMPLE_H__
