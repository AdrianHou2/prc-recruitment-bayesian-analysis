// see if I can template statistic function return type
// maybe change to lambda functions instead of functor objects? may be faster
// to compile:
// g++ -c ./prc1System.cpp ./prc1sim.cpp
// g++ -shared -o gillespielibrary.dll prc1sim.o prc1System.o

#include "prc1System.h"
#include "gillespie.h"

#include <vector>
#include <random>
#include <iostream>
#include <functional>
#include <chrono>

struct RateFunction {
    double initial_binding_rate;
    double singly_bound_unbinding_rate;
    double k0;

    RateFunction(
        double initial_binding_rate_,
        double singly_bound_unbinding_rate_,
        double k0_
    ){
        initial_binding_rate = initial_binding_rate_;
        singly_bound_unbinding_rate = singly_bound_unbinding_rate_;
        k0 = k0_;
    }

    std::vector<double> operator()(const PRC1System& state) {
        std::vector<double> rates;
        int num_crosslinkers = state.num_agents;
        rates.reserve(num_crosslinkers*3);

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

        // determine tables of reaction rates to the left and right (already attached at bottom, attaching to top)
        std::vector<double> left_rates = {0};
        std::vector<double> right_rates = {0};
        double site_difference = std::fmod(state.microtubule_offset, state.site_spacing);
        for (int i = 0; i < 500; i++) {
            double left_spacing = site_difference - (i+1) * state.site_spacing;
            double right_spacing = site_difference + (i) * state.site_spacing;
            double left_E = state.spring_constant/2 * (sqrt(pow(left_spacing, 2) + pow(state.microtubule_seperation, 2)) - state.rest_length);
            double right_E = state.spring_constant/2 * (sqrt(pow(right_spacing, 2) + pow(state.microtubule_seperation, 2)) - state.rest_length);
            double left_rate = k0 * exp(-.5 * left_E / state.k_B_T);
            double right_rate = k0 * exp(-.5 * right_E / state.k_B_T);
            left_rates.push_back(left_rate);
            right_rates.push_back(right_rate);
        }

        // cumulative rates
        std::vector<double> left_cumulative_rates(left_rates.size());
        std::vector<double> right_cumulative_rates(right_rates.size());
        std::partial_sum(left_rates.cbegin(), left_rates.cend(), left_cumulative_rates.begin());
        std::partial_sum(right_rates.cbegin(), right_rates.cend(), right_cumulative_rates.begin());

        for (int reaction_agent = 0; reaction_agent < num_crosslinkers; reaction_agent++) {
            if (!(state.bottom_is_attached[reaction_agent] && state.top_is_attached[reaction_agent])) {
                // find closest left and right indices that are attached on both microtubules
                int left_index = -1;
                int right_index = state.num_sites;
                int prc1_index = -1;
                if (state.top_is_attached[reaction_agent]) {
                    prc1_index = state.top_positions[reaction_agent];
                    for (int i = prc1_index-1; i >= 0; i--) {
                        bool found = false;
                        if (state.top_sites_are_taken[i]) {
                            for (int j = 0; j < state.num_agents; j++) {
                                if (state.top_positions[j] == i && state.bottom_is_attached[j]) {
                                    found = true;
                                    left_index = state.bottom_positions[j];
                                    break;
                                }
                            }
                            if (found) {
                                break;
                            }
                        }
                    }
                    for (int i = prc1_index+1; i < state.num_sites; i++) {
                        bool found = false;
                        if (state.top_sites_are_taken[i]) {
                            for (int j = 0; j < state.num_agents; j++) {
                                if (state.top_positions[j] == i && state.bottom_is_attached[j]) {
                                    found = true;
                                    right_index = state.bottom_positions[j];
                                    break;
                                }
                            }
                            if (found) {
                                break;
                            }
                        }
                    }
                } else {
                    prc1_index = state.bottom_positions[reaction_agent];
                    for (int i = prc1_index-1; i >= 0; i--) {
                        bool found = false;
                        if (state.bottom_sites_are_taken[i]) {
                            for (int j = 0; j < state.num_agents; j++) {
                                if (state.bottom_positions[j] == i && state.top_is_attached[j]) {
                                    found = true;
                                    left_index = state.top_positions[j];
                                    break;
                                }
                            }
                            if (found) {
                                break;
                            }
                        }
                    }
                    for (int i = prc1_index+1; i < state.num_sites; i++) {
                        bool found = false;
                        if (state.bottom_sites_are_taken[i]) {
                            for (int j = 0; j < state.num_agents; j++) {
                                if (state.bottom_positions[j] == i && state.top_is_attached[j]) {
                                    found = true;
                                    right_index = state.top_positions[j];
                                    break;
                                }
                            }
                            if (found) {
                                break;
                            }
                        }
                    }
                }
                
                if (state.bottom_is_attached[reaction_agent]) {
                    double target_pos_nm = prc1_index*state.site_spacing - state.microtubule_offset;
                    int left_available_sites = int((target_pos_nm+site_difference) / state.site_spacing  - left_index - 1);
                    int right_available_sites = right_index - int((target_pos_nm+site_difference) / state.site_spacing);
                    left_available_sites = std::min(left_available_sites, 500);
                    right_available_sites = std::min(right_available_sites, 500);
                    double left_rate = left_cumulative_rates[left_available_sites];
                    double right_rate = right_cumulative_rates[right_available_sites];

                    int non_existent_sites_left = std::max(0.0, -(target_pos_nm - state.sites.back())/state.site_spacing);
                    int non_existent_sites_right = std::max(0.0, -target_pos_nm/state.site_spacing);
                    non_existent_sites_left = std::min(non_existent_sites_left, 500);
                    non_existent_sites_right = std::min(non_existent_sites_right, 500);
                    left_rate -= left_cumulative_rates[non_existent_sites_left];
                    right_rate -= right_cumulative_rates[non_existent_sites_right];

                    left_rate = std::max(0.0, left_rate);
                    right_rate = std::max(0.0, right_rate);

                    rates.push_back(left_rate + right_rate);
                } else {
                    double target_pos_nm = prc1_index*state.site_spacing + state.microtubule_offset;
                    int left_available_sites = int((target_pos_nm-site_difference) / state.site_spacing - left_index);
                    int right_available_sites = right_index - int((target_pos_nm-site_difference) / state.site_spacing) - 1;
                    left_available_sites = std::min(left_available_sites, 500);
                    right_available_sites = std::min(right_available_sites, 500);
                    double left_rate = right_cumulative_rates[left_available_sites]; // flipped since attaching to bottom instead
                    double right_rate = left_cumulative_rates[right_available_sites];

                    int non_existent_sites_left = std::max(0.0, -(target_pos_nm - state.sites.back())/state.site_spacing);
                    int non_existent_sites_right = std::max(0.0, -target_pos_nm/state.site_spacing);
                    non_existent_sites_left = std::min(non_existent_sites_left, 500);
                    non_existent_sites_right = std::min(non_existent_sites_right, 500);
                    left_rate -= right_cumulative_rates[non_existent_sites_left];
                    right_rate -= left_cumulative_rates[non_existent_sites_right];
                    left_rate = std::max(0.0, left_rate);
                    right_rate = std::max(0.0, right_rate);

                    rates.push_back(left_rate + right_rate);
                }

                // determine the position directly above reaction_agent's attached site
                // double target_pos_nm = state.top_is_attached[reaction_agent] ? prc1_index*state.site_spacing + state.microtubule_offset
                                                                            //  : prc1_index*state.site_spacing - state.microtubule_offset;
                
                
            } else {
                rates.push_back(0);
            }
        }

        // detachment rates
        for (int i = 0; i < num_crosslinkers; i++) {
            if (state.bottom_is_attached[i] && state.top_is_attached[i]) {
                double horizontal_stretch = state.site_spacing * (state.top_positions[i] - state.bottom_positions[i]) + state.microtubule_offset;
                double E = state.spring_constant/2 * (sqrt(pow(horizontal_stretch, 2) + pow(state.microtubule_seperation, 2)) - state.rest_length);
                rates.push_back(2 * k0 * exp(.5 * E / state.k_B_T)); // factor of 2 since either head can detach
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

// we should really only have to check about 10 sites on either side with 99.9% accuracy
struct SecondHeadAttachmentFunction {
    double k0;
    SecondHeadAttachmentFunction(double k0_) : k0(k0_) {}
    inline static std::mt19937 gen{std::random_device{}()};
    
    void operator()(PRC1System& state, int reaction_agent) {
        // error check
        if (state.top_is_attached[reaction_agent] == state.bottom_is_attached[reaction_agent]) {
            throw std::runtime_error("Invalid state: PRC1 has either both or no heads attached");
        }
        // find closest left and right indices that are attached on both microtubules
        int left_index = -1;
        int right_index = state.num_sites;
        int prc1_index = -1;
        if (state.top_is_attached[reaction_agent]) {
            prc1_index = state.top_positions[reaction_agent];
            for (int i = prc1_index-1; i >= 0; i--) {
                bool found = false;
                if (state.top_sites_are_taken[i]) {
                    for (int j = 0; j < state.num_agents; j++) {
                        if (state.top_positions[j] == i && state.bottom_is_attached[j]) {
                            found = true;
                            left_index = state.bottom_positions[j];
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
            }
            for (int i = prc1_index+1; i < state.num_sites; i++) {
                bool found = false;
                if (state.top_sites_are_taken[i]) {
                    for (int j = 0; j < state.num_agents; j++) {
                        if (state.top_positions[j] == i && state.bottom_is_attached[j]) {
                            found = true;
                            right_index = state.bottom_positions[j];
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
            }
        } else {
            prc1_index = state.bottom_positions[reaction_agent];
            for (int i = prc1_index-1; i >= 0; i--) {
                bool found = false;
                if (state.bottom_sites_are_taken[i]) {
                    for (int j = 0; j < state.num_agents; j++) {
                        if (state.bottom_positions[j] == i && state.top_is_attached[j]) {
                            found = true;
                            left_index = state.top_positions[j];
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
            }
            for (int i = prc1_index+1; i < state.num_sites; i++) {
                bool found = false;
                if (state.bottom_sites_are_taken[i]) {
                    for (int j = 0; j < state.num_agents; j++) {
                        if (state.bottom_positions[j] == i && state.top_is_attached[j]) {
                            found = true;
                            right_index = state.top_positions[j];
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
            }
        }

        // for each available site, determine rate of attachment
        std::vector<double> attachment_probabilities;
        if (state.top_is_attached[reaction_agent]) {
            for (int site = left_index+1; site < right_index; site++) {
                double horizontal_stretch = state.site_spacing * (state.top_positions[reaction_agent] - site) + state.microtubule_offset;
                double E = state.spring_constant/2 * (sqrt(pow(horizontal_stretch, 2) + pow(state.microtubule_seperation, 2)) - state.rest_length);
                attachment_probabilities.push_back(k0 * exp(-.5 * E / state.k_B_T));
            }
        } else {
            for (int site = left_index+1; site < right_index; site++) {
                double horizontal_stretch = state.site_spacing * (site - state.bottom_positions[reaction_agent]) + state.microtubule_offset;
                double E = state.spring_constant/2 * (sqrt(pow(horizontal_stretch, 2) + pow(state.microtubule_seperation, 2)) - state.rest_length);
                attachment_probabilities.push_back(k0 * exp(-.5 * E / state.k_B_T));
            }
        }

        // in case there are no available sites to attach to
        if (attachment_probabilities.size() == 0) {
            return;
        }

        // choose which site to attach to
        std::partial_sum(attachment_probabilities.cbegin(), attachment_probabilities.cend(), attachment_probabilities.begin());
        double totalRate = attachment_probabilities.back();
        std::uniform_real_distribution<double> uniform_distribution(0.0, totalRate);
        double random_num = uniform_distribution(gen);
        int site_index = 0;
        while (attachment_probabilities[site_index] < random_num) {
            site_index++;
        }
        int selected_site = site_index + left_index + 1;

        // attach selected site
        if (state.top_is_attached[reaction_agent]) {
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

struct StatisticFunction {
    int operator()(PRC1System& state) {
        return state.num_agents;
    }
};

extern "C" {
__declspec(dllexport) double* run_prc1_sim(
    double initial_binding_rate,
    double singly_bound_unbinding_rate,
    double k0,
    double* eval_times_array,
    int num_eval_times
){
    // define relevant constants
    double microtubule_length = 5000;
    double site_spacing = 0.2;
    double offset = 2000;
    double spring_constant = 2;
    double rest_length = 32;
    double k_B_T = 4.1;
    double microtubule_seperation = 32;


    // cast eval_times_array to a vector
    std::vector<double> eval_times(eval_times_array, eval_times_array + num_eval_times);

	// random generator initialization
	static std::random_device rd;
	static std::mt19937 gen(rd());
    static std::uniform_int_distribution random_bool(0, 1);
    static std::uniform_real_distribution uniform_distribution(0.0, 1.0);
    
    // define rate function
    RateFunction rate_function(
        initial_binding_rate,
        singly_bound_unbinding_rate,
        k0);

    // define attachment function
    InitialAttachmentFunction initial_attachment_function;
    SecondHeadAttachmentFunction second_head_attachment_function(k0);
    DetachmentFunction detachment_function;
    TimestepFunction timestep_function;
    StatisticFunction statistic_function;

    PRC1System initial_state(microtubule_length, site_spacing, offset, spring_constant, rest_length, k_B_T, microtubule_seperation);

    std::vector<std::function<void(PRC1System&, int)>> reaction_functions = {
        initial_attachment_function,
        second_head_attachment_function,
        detachment_function
    };

    double end_time = eval_times.back();
    
    std::tuple<std::vector<int>, std::vector<double>> history =
        run_gillespie(
            initial_state,
            rate_function,
            reaction_functions,
            timestep_function,
            statistic_function,
            end_time
        );
    
    std::vector<int> state_history = std::get<0>(history);
    std::vector<double> timesteps = std::get<1>(history);

    int timestep_index = 0;
    std::vector<double> answer = {}; // awful name but idk what else to call it
    
    for (const double& eval_time : eval_times) {
        while (timesteps[timestep_index] <= eval_time) {
            timestep_index++;
        }
        // PRC1System cur_state = state_history[timestep_index-1];
        answer.push_back(state_history[timestep_index-1]);
    }
    double* answer_array = new double[answer.size()];
    std::copy(answer.begin(), answer.end(), answer_array);
    return answer_array;
}
__declspec(dllexport) void free_array(double* ptr) {
    delete[] ptr;
}
}

int main() {
    std::vector<double> eval_times = {1, 2, 3, 4, 5, 90};

    auto start = std::chrono::high_resolution_clock::now();
    double* answer_array = run_prc1_sim(2.83723976, 3.19827888, 10, eval_times.data(), 6);
    std::vector<double> answer(answer_array, answer_array + 6);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << "\n";

    for (const double& n : answer) {
        std::cout << n << "\n";
    }
    delete answer_array;
}