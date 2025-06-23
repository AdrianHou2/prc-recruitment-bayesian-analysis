#include <vector>

class PRC1System {
public:
    // positions where the top and bottom heads of each microtubule are attached.
    // positions are relative to the top and bottom microtubules respectively;
    // the left side of each microtubule is the 0 position
    std::vector<double> top_positions;
    std::vector<double> bottom_positions;
    
    // whether the top or bottom heads of the microtubule is attached
    std::vector<bool> top_is_attached;
    std::vector<bool> bottom_is_attached;

    // keeps track whether each of the top or bottom sites are taken
    std::vector<bool> top_sites_are_taken;
    std::vector<bool> bottom_sites_are_taken;

    // horizontal displacement of the top microtubule's left end from the bottom microtubule's left end
    double microtubule_offset;

    int num_sites;
    int num_crosslinkers;

    PRC1System(double microtubule_length, double site_spacing, double offset) {
        microtubule_offset = offset;

        top_sites_are_taken = {0};
        bottom_sites_are_taken = {0};
        double cur = 0;
        while (cur < microtubule_length) {
            top_sites_are_taken.push_back(site_spacing);
            bottom_sites_are_taken.push_back(site_spacing);
            cur += site_spacing;
        }
    };
    
};