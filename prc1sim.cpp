#include "prc1System.h"
#include "gillespie.h"

#include <vector>
#include <random>
#include <iostream>
#include <functional>

struct RateFunction {
    double initial_binding_rate;
    double second_head_binding_rate;
    double singly_bound_unbinding_rate;
    double doubly_bound_unbinding_rate;

    RateFunction(
        double initial_binding_rate_,
        double second_head_binding_rate_,
        double singly_bound_unbinding_rate_,
        double doubly_bound_unbinding_rate_
    ){
        initial_binding_rate = initial_binding_rate_;
        second_head_binding_rate = second_head_binding_rate_;
        singly_bound_unbinding_rate = singly_bound_unbinding_rate_;
        doubly_bound_unbinding_rate = doubly_bound_unbinding_rate_;
    }

    std::vector<double> operator()(const PRC1System& state) {
        std::vector<double> rates = {};
        int num_crosslinkers = state.num_agents;

        if (state.num_agents == 0) {
            return {initial_binding_rate, 0, 0};
        }

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

        return rates;
    }
};

struct InitialAttachmentFunction {
	inline static std::mt19937 gen{std::random_device{}()};
    void operator()(PRC1System& state, int reaction_agent) {
        state.num_agents++;

        // determine which site to attach to (equal probability for each site)
        int num_open_sites = state.sites.size()*2 - state.num_heads_attached_top - state.num_heads_attached_bottom;
        std::uniform_int_distribution integer_distribution(0, num_open_sites-1);
        int attachment_number = integer_distribution(gen);

        // determine which index that site lines up with
        int attachment_index = 0;
        while (attachment_index < state.num_sites && attachment_number > 0) {
            if (!state.top_sites_are_taken[attachment_index]) {
                attachment_number -= 1;
            }
            attachment_index++;
        }
        if (attachment_number < 0) {
            state.top_is_attached.push_back(true);
            state.bottom_is_attached.push_back(false);
            state.top_positions.push_back(attachment_index);
            state.bottom_positions.push_back(-1);
            state.top_sites_are_taken[attachment_index] = true;
            state.num_heads_attached_top++;
            return;
        }
        attachment_index = 0;
        while (attachment_index < state.num_sites && attachment_number > 0) {
            if (!state.bottom_sites_are_taken[attachment_index]) {
                attachment_number -= 1;
            }
            attachment_index++;
        }
        state.bottom_is_attached.push_back(true);
        state.top_is_attached.push_back(false);
        state.bottom_positions.push_back(attachment_index);
        state.top_positions.push_back(-1);
        state.bottom_sites_are_taken[attachment_index] = true;
        state.num_heads_attached_bottom++;
    }
};

// fixed but not really optimize, might be a little slow
struct SecondHeadAttachmentFunction {
    inline static std::mt19937 gen{std::random_device{}()};
    
    void operator()(PRC1System& state, int reaction_agent) {
        // error check
        if (state.top_is_attached[reaction_agent] == state.bottom_is_attached[reaction_agent]) {
            throw std::runtime_error("Invalid state: PRC1 has either both or no heads attached");
        }
        // get position
        int first_head_pos_index;
        double first_head_pos_nm;
        bool is_top_attached = state.top_is_attached[reaction_agent];
        if (is_top_attached) {
            first_head_pos_index = state.top_positions[reaction_agent];
            first_head_pos_nm = state.sites[first_head_pos_index];
        } else {
            first_head_pos_index = state.bottom_positions[reaction_agent];
            first_head_pos_nm = state.sites[first_head_pos_index];
        }
        double target_pos_nm = is_top_attached ? first_head_pos_nm - state.microtubule_offset : first_head_pos_nm + state.microtubule_offset;
        // closest site index on the opposite microtubule
        int target_index = 0;
        double min_diff = std::abs(state.sites[0] - target_pos_nm);
        for (int i = 1; i < state.num_sites; i++) {
            double diff = std::abs(state.sites[i] - target_pos_nm);
            if (diff < min_diff) {
                min_diff = diff;
                target_index = i;
            }
        }
        // nearest occupied sites (left and right) on the opposite microtubule
        int min_index = -1; // left bound 
        int max_index = state.num_sites; // right bound 
        const std::vector<bool>& opposite_sites = is_top_attached ? state.bottom_sites_are_taken : state.top_sites_are_taken;
        // left index
        for (int i = target_index; i >= 0; i--) {
            if (opposite_sites[i]) {
                min_index = i;
                break;
            }
        }
        // right index
        for (int i = target_index; i < state.num_sites; i++) {
            if (opposite_sites[i]) {
                max_index = i;
                break;
            }
        }
        // available sites 
        std::vector<int> available_sites;
        for (int i = min_index + 1; i < max_index; i++) {
            if (!opposite_sites[i]) {
                available_sites.push_back(i);
            }
        }
        // select an available site
        std::uniform_int_distribution<> dist(0, available_sites.size() - 1);
        int selected_site = available_sites[dist(gen)];
        // update state
        if (is_top_attached) {
            // attach bottom head
            state.bottom_is_attached[reaction_agent] = true;
            state.bottom_positions[reaction_agent] = selected_site;
            state.bottom_sites_are_taken[selected_site] = true;
            state.num_heads_attached_bottom++;
        } else {
            // attach top head
            state.top_is_attached[reaction_agent] = true;
            state.top_positions[reaction_agent] = selected_site;
            state.top_sites_are_taken[selected_site] = true;
            state.num_heads_attached_top++;
        }
    }
};

struct DetachmentFunction {
    static inline std::mt19937 gen{std::random_device{}()};
    static inline std::uniform_int_distribution random_bool{0,1};

    void operator()(PRC1System& state, int reaction_agent) {
;
        if (state.bottom_is_attached[reaction_agent] && state.top_is_attached[reaction_agent]) {
            if (random_bool(gen)) {
                // detach top head
                state.top_sites_are_taken[state.top_positions[reaction_agent]] = false;
                state.top_positions[reaction_agent] = -1;
                state.top_is_attached[reaction_agent] = false;
                state.num_heads_attached_top--;
            } else {   
                // detach bottom head
                state.bottom_sites_are_taken[state.bottom_positions[reaction_agent]] = false;
                state.bottom_positions[reaction_agent] = -1;
                state.bottom_is_attached[reaction_agent] = false;
                state.num_heads_attached_bottom--;
            }
        } else {
            if (state.bottom_is_attached[reaction_agent]) {
                state.bottom_sites_are_taken[state.bottom_positions[reaction_agent]] = false;
                state.num_heads_attached_bottom--;
            } else {
                state.top_sites_are_taken[state.top_positions[reaction_agent]] = false;
                state.num_heads_attached_top--;
            }
            state.bottom_is_attached.erase(state.bottom_is_attached.begin() + reaction_agent);
            state.top_is_attached.erase(state.top_is_attached.begin() + reaction_agent);
            state.bottom_positions.erase(state.bottom_positions.begin() + reaction_agent);
            state.top_positions.erase(state.top_positions.begin() + reaction_agent);
            state.num_agents--;
        }
    }
};

struct TimestepFunction {
    void operator()(PRC1System& state, double& time, double dt) {
        time += dt;
    }
};

std::vector<double> run_prc1_sim(
    double initial_binding_rate,
    double second_head_binding_rate,
    double singly_bound_unbinding_rate,
    double doubly_bound_unbinding_rate,
    const std::vector<double>& eval_times
){
	// random generator initialization
	static std::random_device rd;
	static std::mt19937 gen(rd());
    static std::uniform_int_distribution random_bool(0, 1);
    static std::uniform_real_distribution uniform_distribution(0.0, 1.0);
    
    // define rate function
    RateFunction rate_function(
        initial_binding_rate,
        second_head_binding_rate,
        singly_bound_unbinding_rate,
        doubly_bound_unbinding_rate);

    // define attachment function
    InitialAttachmentFunction initial_attachment_function;

    // FIX THIS TO BE RIGHT (it's currently just always attaching the second head to the 0 position on the other head)
    SecondHeadAttachmentFunction second_head_attachment_function;

    DetachmentFunction detachment_function;

    TimestepFunction timestep_function;

    double microtubule_legnth = 5000;
    double site_spacing = 0.2;
    double offset = 2000;
    PRC1System initial_state(microtubule_legnth, site_spacing, offset);

    std::vector<std::function<void(PRC1System&, int)>> reaction_functions = {
        initial_attachment_function,
        second_head_attachment_function,
        detachment_function
    };

    double end_time = eval_times.back();
    
    std::tuple<std::vector<PRC1System>, std::vector<double>> history =
        run_gillespie(
            initial_state,
            rate_function,
            reaction_functions,
            timestep_function,
            end_time
        );
    
    std::vector<PRC1System> state_history = std::get<0>(history);
    std::vector<double> timesteps = std::get<1>(history);

    int timestep_index = 0;
    std::vector<double> answer = {}; // awful name but idk what else to call it
    for (const double& eval_time : eval_times) {
        while (timesteps[timestep_index] <= eval_time) {
            timestep_index++;
        }
        PRC1System cur_state = state_history[timestep_index-1];
        answer.push_back(cur_state.num_agents);
    }
    return answer;
}

int main() {
    std::vector<double> eval_times = {1, 2, 3, 4, 5};
    std::vector<double> answer = run_prc1_sim(2.83723976, 620.3984088, 3.19827888, 2.45246624, eval_times);
    for (const double& n : answer) {
        std::cout << n << "\n";
    }
}