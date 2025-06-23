#include "gillespie.h"

#include <vector>
#include <random>
#include <functional>
#include <numeric>
#include <tuple>

// modifies state and time to the next timestep
template <typename State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
void gillespie_step(
	State& state,
	double& time,
	RateFunction rate_function,
	std::vector<ReactionFunction> reaction_function,
	TimestepFunction timestep_function,
	int num_agents
){
	// random generator initialization
	static std::random_device rd;
	static std::mt19937 gen(rd());
	
	// gets cumulative rates
	std::vector<double> cumulative_rates = rate_function(state);
	std::partial_sum(cumulative_rates.cbegin(), cumulative_rates.cend(), cumulative_rates.begin());
	double totalRate = cumulative_rates.back();

	// calculate timestep and increment time
	std::exponential_distribution<double> exponential_distribution(1/totalRate);
	double dt = exponential_distribution(gen);
	timestep_function(state, time, dt);

	// determine index of reaction
	std::uniform_real_distribution<double> uniform_distribution(0.0, totalRate);
	double random_num = uniform_distribution(gen);
	int reaction_index = 0;
	while (cumulative_rates[reaction_index] < random_num) {
		reaction_index++;
	}

	// determine which function and agent the index lines up with
	ReactionFunction reactionFunction = reaction_function[reaction_index / num_agents];
	int reaction_agent = reaction_index % num_agents;

	// do reaction
	reactionFunction(state, reaction_agent);
}

template <class State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
std::tuple<std::vector<State>, std::vector<double>> runGillespie(
	State initial_state,
	RateFunction rate_function,
	std::vector<ReactionFunction> reaction_function,
	TimestepFunction timestep_function,
	double end_time
){
	State state = initial_state;
	double time = 0;
	std::vector<double> times = {0};
	std::vector<State> states = {initial_state};

	while (time < end_time) {
		gillespieStep(state, time, rate_function, reaction_function, timestep_function, state.num_agents);
		times.push_back(time);
		states.push_back(State(state));
	}

	return std::make_tuple(states, times);
}