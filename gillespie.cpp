#include "gillespie.h"

#include <vector>
#include <random>
#include <functional>
#include <numeric>
#include <tuple>

// modifies state and time to the next timestep
template <typename State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
void gillespieStep(
	State& state,
	double& time,
	RateFunction rateFunction,
	std::vector<ReactionFunction> reactionFunctions,
	TimestepFunction timestepFunction,
	int numAgents
){
	// random generator initialization
	static std::random_device rd;
	static std::mt19937 gen(rd());
	
	// gets cumulative rates
	std::vector<double> cumulativeRates = rateFunction(state);
	std::partial_sum(cumulativeRates.cbegin(), cumulativeRates.cend(), cumulativeRates.begin());
	double totalRate = cumulativeRates.back();

	// calculate timestep and increment time
	std::exponential_distribution<double> exponentialDistribution(1/totalRate);
	double dt = exponentialDistribution(gen);
	timestepFunction(state, time, dt);

	// determine index of reaction
	std::uniform_real_distribution<double> uniformDistribution(0, totalRate);
	double randomNum = uniformDistribution(gen);
	int reactionIndex = 0;
	while (cumulativeRates[reactionIndex] < randomNum) {
		reactionIndex++;
	}

	// determine which function and agent the index lines up with
	ReactionFunction reactionFunction = reactionFunctions[reactionIndex / numAgents];
	int reactionAgent = reactionIndex % numAgents;

	// do reaction
	reactionFunction(state, reactionAgent);
}

template <class State, typename RateFunction, typename ReactionFunction, typename TimestepFunction>
std::tuple<std::vector<State>, std::vector<double>> runGillespie(
	State initialState,
	RateFunction rateFunction,
	std::vector<ReactionFunction> reactionFunctions,
	TimestepFunction timestepFunction,
	double endTime,
	int numAgents
){
	State state = initialState;
	double time = 0;
	std::vector<double> times = {0};
	std::vector<State> states = {initialState};

	while (time < endTime) {
		gillespieStep(state, time, rateFunction, reactionFunctions, timestepFunction, numAgents);
		times.push_back(time);
		states.push_back(State(state));
	}

	return std::make_tuple(states, times);
}