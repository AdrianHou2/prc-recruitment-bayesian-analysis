# prc-recruitment-bayesian-analysis
## Python version:
- prc1_state.py contains the State class, which is a container class that handles all the operation for the prc1-microtubule system
- prc1.py contains the Prc1 class, which represents a single prc1 and interfaces with the State class
- gillespie.py has generic driver code for gillespie algorithm (and defines relevant constants)
- run_gillespie.py calls the generic driver code with the specific functions for the prc1 system

to run the simulation, use:

`from run_gillespie import run_gillespie_prc1`

`run_gillespie_prc1(initial_binding_rate, singly_bound_detachment_rate, k0, end_time)`




## Cpp version:
- prc1System.h and prc1System.cpp define the PRC1System class, which is the state object
- gillespie.h defines a very generic gillespie algorithm
- prc1sim.cpp defines and calls the gillespie algorithm for specifically the prc1 case (this is where the main function is)
