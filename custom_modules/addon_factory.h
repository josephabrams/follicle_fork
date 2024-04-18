#ifndef __ADDON_FACTORY_H__
#define __ADDON_FACTORY_H__
#include <vector>
#include <memory>
#include "../core/PhysiCell_cell.h"
using namespace PhysiCell;
using namespace BioFVM;
class Addon;
class Base_Addon_Class;

class Addon_Factory{
  public:
    Addon* m_addon;
    virtual Base_Addon_Class* Factory_Method(Cell* pCell)=0;
    Base_Addon_Class* Create_Addon_Instance(Cell* pCell);
    virtual ~Addon_Factory();
};

#endif // !__ADDON_FACTORY_H__

