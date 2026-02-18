## file to run gillespie simulation on specifically the prc1_state.py State object
from gillespie import run_gillespie
from prc1_state import State
from copy import deepcopy

def run_gillespie_prc1(initial_binding_rate, singly_bound_detachment_rate, k0, end_time):
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
                    rest_length, k_B_T, microtubule_separation, singly_bound_detachment_rate, k0)

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

    max_timestep = 0.1

    return run_gillespie(initial_state, end_time, rate_function, reaction_functions,
                        statistic_function, timestep_function, max_timestep)