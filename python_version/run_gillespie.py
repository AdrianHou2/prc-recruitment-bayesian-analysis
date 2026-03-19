## file to run gillespie simulation on specifically the prc1_state.py State object
from gillespie import run_gillespie, gillespie_step
from prc1_state import State
from copy import deepcopy
import numpy as np


def run_gillespie_prc1(initial_binding_rate, singly_bound_detachment_rate,
                       base_double_attachment_rate, base_double_detachment_rate, end_time):
    # define rate params (now passed in)
    # initial_binding_rate = 2.8
    # singly_bound_detachment_rate = 3.2
    # k0 = 10

    # define (initial) state params
    microtubule_length = 5000.
    site_spacing = 0.2
    microtubule_offset = 2000
    spring_constant = 2
    rest_length = 32
    k_B_T = 4.1
    microtubule_separation = 32
    params = (microtubule_length, site_spacing, microtubule_offset, spring_constant,
              rest_length, k_B_T, microtubule_separation, singly_bound_detachment_rate,
              base_double_attachment_rate, base_double_detachment_rate)

    # define gillespie params and funcs
    initial_state = State(*params)
    # end_time = 10

    # returns (single attach rates, double attach rates, detach rates)
    def rate_function(state):
        if len(state) == 0:
            return [initial_binding_rate, 0, 0]
        rates = [initial_binding_rate] + [0]*(len(state)-1)
        for prc1 in state:
            rates.append(prc1.double_attachment_rate)
        for prc1 in state:
            rates.append(prc1.detachment_rate)
        return rates

    # define list of reaction functions
    def single_attach_func(state, index):
        # print("single attaching", index)
        state.single_attach_prc1()
        # print(state)

    def double_attach_func(state, index):
        # print("double attaching", index)
        state.double_attach_prc1(index)
        # print(state)

    def detach_func(state, index):
        # print("detaching", index)
        state.detach_prc1(index)
        # print(state)

    reaction_functions = [single_attach_func, double_attach_func, detach_func]

    # define statistic function, returns list of num of prc1 at every timestep
    def statistic_function(state, stat_list, time):
        if stat_list is None:
            stat_list = []
        # stat_list.append(deepcopy(state))
        stat_list.append(len(state))


    # define timestep function
    def timestep_function(state, time, dt):
        return time + dt

    max_timestep = np.inf

    return run_gillespie(initial_state, end_time, rate_function, reaction_functions,
                        statistic_function, timestep_function, max_timestep)


def run_gillespie_prc1_on_grid(initial_binding_rate, singly_bound_detachment_rate, base_double_attachment_rate,
                               base_double_detachment_rate, times_obs,
                               max_steps=200_000, max_timestep=np.inf):
    """
    Runs the *same* spatial PRC1 Gillespie model, but only records the observable y(t)=len(state)
    on the user-supplied observation grid times_obs (previous-value / right-continuous sampling).

    This avoids storing every event time/state and avoids per-path np.searchsorted.

    Returns
    -------
    y : np.ndarray, shape (len(times_obs),)
        y[k] = len(state) at time times_obs[k]
    """
    times_obs = np.asarray(times_obs, dtype=float)
    if times_obs.ndim != 1:
        raise ValueError("times_obs must be a 1D array")
    if len(times_obs) == 0:
        return np.zeros(0, dtype=float)
    if np.any(np.diff(times_obs) < 0):
        raise ValueError("times_obs must be sorted nondecreasing")

    end_time = float(times_obs[-1])

    # define (initial) state params (same as run_gillespie_prc1)
    microtubule_length = 5000.
    site_spacing = 0.2
    microtubule_offset = 2000
    spring_constant = 2
    rest_length = 32
    k_B_T = 4.1
    microtubule_separation = 32
    params = (microtubule_length, site_spacing, microtubule_offset, spring_constant,
              rest_length, k_B_T, microtubule_separation, singly_bound_detachment_rate,
              base_double_attachment_rate, base_double_detachment_rate)

    state = State(*params)

    # returns (single attach rates, double attach rates, detach rates)
    def rate_function(state):
        if len(state) == 0:
            return [initial_binding_rate, 0, 0]
        rates = [initial_binding_rate] + [0]*(len(state)-1)
        for prc1 in state:
            rates.append(prc1.double_attachment_rate)
        for prc1 in state:
            rates.append(prc1.detachment_rate)
        return rates

    # define list of reaction functions
    def single_attach_func(state, index):
        state.single_attach_prc1()

    def double_attach_func(state, index):
        state.double_attach_prc1(index)

    def detach_func(state, index):
        state.detach_prc1(index)

    reaction_functions = [single_attach_func, double_attach_func, detach_func]

    # define timestep function
    def timestep_function(state, time, dt):
        return time + dt

    y = np.zeros(len(times_obs), dtype=float)
    idx = 0
    t = 0.0
    current = float(len(state))

    # fill any observation times at/behind t=0
    while idx < len(times_obs) and times_obs[idx] <= t + 1e-12:
        y[idx] = current
        idx += 1

    steps = 0
    while idx < len(times_obs) and t < end_time and steps < int(max_steps):
        t = gillespie_step(state, t, rate_function, reaction_functions,
                           timestep_function=timestep_function,
                           max_timestep=max_timestep)
        current = float(len(state))

        # fill obs times up to current time (right-continuous sampling, matches searchsorted(..., 'right')-1)
        while idx < len(times_obs) and times_obs[idx] <= t + 1e-12:
            y[idx] = current
            idx += 1

        steps += 1

    # fill remaining with last value
    if idx < len(times_obs):
        y[idx:] = current

    return y