#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <vector>
#include <tuple>

template <class State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
std::tuple<std::vector<State>, std::vector<double>> runGillespie(
	State initialState,
	RateFunction rateFunction,
	std::vector<ReactionFunction> reactionFunctions,
	TimestepFunction timestepFunction,
	double endTime,
	int numAgents
);

#endif