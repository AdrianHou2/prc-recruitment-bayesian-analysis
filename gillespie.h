#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <vector>
#include <tuple>

template <class State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
std::tuple<std::vector<State>, std::vector<double>> run_gillespie(
	State initial_state,
	RateFunction rate_function,
	std::vector<ReactionFunction> reaction_functions,
	TimestepFunction timestep_function,
	double end_time
);

#endif