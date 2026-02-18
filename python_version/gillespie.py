import numpy as np

def default_timestep_function(state, time, dt):
    """
    default timestep function which just advances time by dt
    """
    return time + dt

def default_statistic_function(state, statistics):
    """
    default statistic function which records full state at each timepoint (very inefficient for large systems)
    """
    if statistics is None:
        statistics = [state.copy()]
    else:
        statistics.append(state.copy())

def gillespie_step(state, time, rate_function, reaction_functions, 
                   timestep_function=default_timestep_function, 
                   max_timestep=np.inf):
    """
    gillespie step for a system with multiple agents and reactions. modifies state in place, returns time after step

    :param state: state of the system at current time
    :param time: current time of the system
    :param rate_function: function(state) -> list of cumulative rates for each reaction for each agent
        e.g. (r1 a1, r1 a2, ..., r2 a1, r2 a2, ...)
    :param reaction_functions: list of functions(state, agent) which modify the state corresponding to each reaction
    :param timestep_function: function(state, time, dt) -> new time after timestep, also modifies state if needed
    :param num_agents: number of agents in the system
    :return time: new time after gillespie step
    """
    # number of agents
    num_agents = len(state) if len(state) > 0 else 1

    # get cumulative rates and total rate
    cumulative_rates = np.cumsum(rate_function(state))
    total_rate = cumulative_rates[-1]


    # determine next time something happens
    if total_rate <= 0:
        print("Total rate <= 0!")
        print(state)
        print("cum rates:", cumulative_rates)
    dt = np.random.exponential(1 / total_rate)
    
    # do timestep (if dt > max_timestep, nothing happens and we skip to the next timestep)
    if dt > max_timestep:
        return timestep_function(state, time, max_timestep)
    new_time = timestep_function(state, time, dt)

    # determine which reaction occurs and to which agent
    random_value = np.random.uniform(0, 1) * total_rate
    index = np.searchsorted(cumulative_rates, random_value)
    # print("state:", state, "\n\nrates:", cumulative_rates, "\nchoice:", index)
    reaction_index = index // num_agents
    agent = index % num_agents
    reaction_function = reaction_functions[reaction_index]

    # do reaction
    reaction_function(state, agent)
    return new_time

def run_gillespie(initial_state, end_time, rate_function, reaction_functions, 
                  statistic_function=default_statistic_function,
                  timestep_function=default_timestep_function,
                  max_timestep=np.inf):
    """
    run a gillespie simulation until end_time

    :param initial_state: initial state of the system
    :param end_time: time to end the simulation
    :param rate_function: function(state) -> list of cumulative rates for each reaction for each agent
        e.g. (r1 a1, r1 a2, ..., r2 a1, r2 a2, ...)
    :param reaction_function: list of functions(state, agent) which modify the state corresponding to each reaction
    :param statistic_function: function(state, statistic object, time) updates statistic object with current state.
                                calls function(state, None) to initialize
    :param timestep_function: function(state, time, dt) -> new time after timestep, also modifies state if needed
    :param max_timestep: maximum timestep allowed before forcing a timestep
    :return statistics: array of recorded statistics
    :return times: array of times at which statistics were recorded
    """
    state = initial_state
    time = 0.0
    times = [time]
    statistics = [statistic_function(state, None, time)]
    while time < end_time:
        # try:
        time = gillespie_step(state, time, rate_function, reaction_functions, timestep_function, max_timestep)
        # except:
            # return statistics, np.array(times)
        times.append(time)
        statistic_function(state, statistics, time)
    return statistics, np.array(times)