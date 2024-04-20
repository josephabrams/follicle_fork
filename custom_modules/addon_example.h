#ifndef __ADDON_EXAMPLE_H__
#define __ADDON_EXAMPLE_H__
#include "./addon.h"
#include "base_addon.h"
class Example_Addon_Class: public Base_Addon_Class{
  public:
    Example_Addon_Class(Cell* pCell);
    Example_Addon_Class();
    int example_int;
    void update_state() override;
    void on_division() override;
};


class Example_Addon_Creator : public Addon_Factory{
  Base_Addon_Class* Factory_Method(Cell* pCell) override;
  Base_Addon_Class* blank_factory_method() override;
};
#endif //__ADDON_EXAMPLE_H__
