#ifndef PRC1SYSTEM_H
#define PRC1SYSTEM_H

#include <vector>

class PRC1System {
public:
    // positions are relative to the top and bottom microtubules respectively;
    // the left side of each microtubule is the 0 position
    std::vector<double> top_positions;
    std::vector<double> bottom_positions;
    
    std::vector<bool> top_is_attached;
    std::vector<bool> bottom_is_attached;

    int num_crosslinkers;
    
};

#endif