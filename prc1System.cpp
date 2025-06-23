#include <vector>

class PRC1System {
public:
    // positions where the top and bottom heads of each microtubule are attached.
    // positions are relative to the top and bottom microtubules respectively;
    // the left side of each microtubule is the 0 position
    std::vector<int> top_positions;
    std::vector<int> bottom_positions;
    
    // whether the top or bottom heads of the microtubule is attached
    std::vector<bool> top_is_attached;
    std::vector<bool> bottom_is_attached;

    // keeps track whether each of the top or bottom sites are taken
    std::vector<double> sites;
    std::vector<bool> top_sites_are_taken;
    std::vector<bool> bottom_sites_are_taken;

    // horizontal displacement of the top microtubule's left end from the bottom microtubule's left end
    double microtubule_offset;

    int num_sites;
    int num_agents;
    int num_heads_attached_top;
    int num_heads_attached_bottom;


    PRC1System(double microtubule_length, double site_spacing, double offset) {
        microtubule_offset = offset;
        num_agents = 0;
        num_heads_attached_top = 0;
        num_heads_attached_bottom = 0;

        top_sites_are_taken = {false};
        bottom_sites_are_taken = {false};
        sites = {0};
        double cur_site = 0;
        while (cur_site < microtubule_length) {
            sites.push_back(cur_site);
            top_sites_are_taken.push_back(false);
            bottom_sites_are_taken.push_back(false);
            cur_site += site_spacing;
        }
        num_sites = sites.size();
    };
    
};