#include "prc1System.h"
#include "gillespie.h"

#include <vector>

std::vector<double> run_prc1_sim(
    double initial_binding_rate,
    double second_head_binding_rate,
    double singly_bound_unbinding_rate,
    double doubly_bound_unbinding_rate,
    const std::vector<double>& eval_times
){
    // define rate function
    auto rate_function = [
        initial_binding_rate,
        second_head_binding_rate,
        singly_bound_unbinding_rate,
        doubly_bound_unbinding_rate]

        (const PRC1System& state) {
        std::vector<double> rates = {};
        int num_crosslinkers = state.num_crosslinkers;

        // 1st head attachment rate (since there's no existing prc1 this reaction happens to,
        // we just give each existing prc1 a equal probability of this happening to them and just ignoring
        // which one it happens to)
        for (int i = 0; i < num_crosslinkers; i++) {
            rates.push_back(initial_binding_rate/num_crosslinkers);
        }

        // 2nd head attachment rates
        for (int i = 0; i < num_crosslinkers; i++) {
            rates.push_back(second_head_binding_rate);
        }

        // detachment rates
        for (int i = 0; i < num_crosslinkers; i++) {
            if (state.bottom_is_attached[i] && state.top_is_attached[i]) {
                rates.push_back(doubly_bound_unbinding_rate);
            } else {
                rates.push_back(singly_bound_unbinding_rate);
            }
        }
    };

    // define attachment function
    auto initial_attachment_function = [](PRC1System& state, int reaction_agent) {
        // attach head somehow, I'm going to sleep now
    };

    // TODO: define the rest of the reaction functions (2nd head attach and detachment)
    // and then call the function and see if it even works right (it probably doesn't yet)
}