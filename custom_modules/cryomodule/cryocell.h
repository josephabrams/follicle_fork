
#ifndef __PhysiMeSS_agent_h__
#define __PhysiMeSS_agent_h__

#include "../../core/PhysiCell.h"
#include <vector>
class Cryocell : public PhysiCell::Cell {

    private:
    public:
    
    std::vector <int> neighbour_occupied_voxels; 
    std::vector <int> cell_boundary_voxels;
    
    Cryocell();
    ~Cryocell(){};
    
    void custom_update_volume();
    void two_parameter_single_step();
    
};


#endif
