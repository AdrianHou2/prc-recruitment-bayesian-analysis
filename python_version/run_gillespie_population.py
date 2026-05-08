"""
Population-level Gillespie simulator for PRC1.

This module provides a fast population-level model (2 states: n1, n2) that can serve
as Stage 1 in a two-stage ABC scheme. The population model is ~100-1000× faster than
the spatial model, making it suitable for quick parameter exploration.

Model:
    free ⇌ n₁ ⇌ n₂  (two binding states)
    
Reactions:
    r1: ∅ → n₁  with rate a₁(X) = R₀₁(1 + ca·n₁ + 2ca₂·n₂)(1 - (n₁+2n₂)/nₛ)
    r2: n₁ → ∅  with rate a₂(X) = k₁₀·n₁(1 - cd₁·n₁ - 2cd₂·n₂)
    r3: n₁ → n₂ with rate a₃(X) = k₁₂·n₁(1 + ca₁·n₁ + 2ca₂·n₂)(1 - (n₁+2n₂)/nₛ)
    r4: n₂ → n₁ with rate a₄(X) = k₂₁·n₂(1 - cd₁·n₁ - 2cd₂·n₂)

Observable: N(t) = n₁(t) + n₂(t)
"""

import math
import numpy as np
from typing import Tuple


def gillespie_population_one_path(
    R01: float,
    k10: float,
    k12: float,
    k21: float,
    times_obs: np.ndarray,
    rng: np.random.Generator,
    ns: float = 100.0,
    ca1: float = 0.0,
    ca2: float = 0.0,
    cd1: float = 0.0,
    cd2: float = 0.0,
    n1_init: int = 1,
    n2_init: int = 0,
    clip_negative_factors: bool = True,
    max_steps: int = 200_000,
) -> np.ndarray:
    """
    Simulate one path of the population model and record N(t) = n₁ + n₂ at observation times.
    
    Parameters
    ----------
    R01 : float
        Attachment rate for free → n₁ transition (per capita)
    k10 : float
        Detachment rate for n₁ → free transition
    k12 : float
        Transition rate for n₁ → n₂
    k21 : float
        Transition rate for n₂ → n₁
    times_obs : np.ndarray
        1D array of observation times (must be sorted)
    rng : np.random.Generator
        Random number generator (thread-safe)
    ns : float
        Number of binding sites (finite population constraint)
    ca1 : float
        Attachment cooperativity for singly-bound PRC1
    ca2 : float
        Attachment cooperativity for doubly-bound PRC1
    cd1 : float
        Detachment cooperativity for singly-bound PRC1
    cd2 : float
        Detachment cooperativity for doubly-bound PRC1
    n1_init : int
        Initial number of singly-bound PRC1
    n2_init : int
        Initial number of doubly-bound PRC1
    clip_negative_factors : bool
        If True, clip factors to [0, ∞) to prevent negative rates
    max_steps : int
        Maximum number of Gillespie steps to prevent infinite loops
    
    Returns
    -------
    np.ndarray
        N(t) sampled at times_obs
    """
    times_obs = np.asarray(times_obs, dtype=float)
    T = len(times_obs)
    out = np.empty(T, dtype=float)
    
    t = 0.0
    n1 = int(n1_init)
    n2 = int(n2_init)
    i_obs = 0
    step_count = 0
    
    # Output all observations before t=0
    while i_obs < T and times_obs[i_obs] <= t:
        out[i_obs] = float(n1 + n2)
        i_obs += 1
    
    # Main Gillespie loop
    while i_obs < T and step_count < max_steps:
        step_count += 1
        
        # Occupancy
        H = n1 + 2 * n2
        
        # Cooperative factors
        bind_factor = 1.0 + ca1 * n1 + 2 * ca2 * n2
        unbind_factor = 1.0 - cd1 * n1 - 2 * cd2 * n2
        space_factor = 1.0 - (H / ns)
        
        # Clip negative factors if requested
        if clip_negative_factors:
            bind_factor = max(bind_factor, 0.0)
            unbind_factor = max(unbind_factor, 0.0)
            space_factor = max(space_factor, 0.0)
        
        # Reaction rates
        a1 = R01 * bind_factor * space_factor  # free → n₁
        a2 = k10 * n1 * unbind_factor         # n₁ → free
        a3 = k12 * n1 * bind_factor * space_factor  # n₁ → n₂
        a4 = k21 * n2 * unbind_factor         # n₂ → n₁
        
        a0 = a1 + a2 + a3 + a4
        
        # If total rate is zero, system is stuck
        if a0 <= 0.0:
            out[i_obs:] = float(n1 + n2)
            break
        
        # Draw waiting time (exponential)
        u = rng.random()
        dt = -math.log(max(u, 1e-300)) / a0
        t_next = t + dt
        
        # Record all observations in (t, t_next]
        while i_obs < T and times_obs[i_obs] <= t_next:
            out[i_obs] = float(n1 + n2)
            i_obs += 1
        
        # Execute reaction at t_next
        u_rxn = rng.random() * a0
        if u_rxn < a1:
            n1 += 1
        elif u_rxn < a1 + a2:
            n1 = max(0, n1 - 1)
        elif u_rxn < a1 + a2 + a3:
            n1 = max(0, n1 - 1)
            n2 += 1
        else:  # u_rxn < a0
            n2 = max(0, n2 - 1)
            n1 += 1
        
        t = t_next
    
    # Fill remaining observations
    if i_obs < T:
        out[i_obs:] = float(n1 + n2)
    
    return out


def run_gillespie_population_on_grid(
    R01: float,
    k10: float,
    k12: float,
    k21: float,
    times_obs: np.ndarray,
    n_runs: int = 10,
    seed: int = 0,
    ns: float = 100.0,
    ca1: float = 0.0,
    ca2: float = 0.0,
    cd1: float = 0.0,
    cd2: float = 0.0,
    n1_init: int = 1,
    n2_init: int = 0,
    clip_negative_factors: bool = True,
    max_steps: int = 200_000,
) -> np.ndarray:
    """
    Run multiple paths of the population model.
    
    Parameters
    ----------
    R01, k10, k12, k21 : float
        Base rates for population model
    times_obs : np.ndarray
        Observation times
    n_runs : int
        Number of independent paths to simulate
    seed : int
        Random seed
    ns : float
        Number of binding sites
    ca1, ca2, cd1, cd2 : float
        Cooperativity parameters
    n1_init, n2_init : int
        Initial state
    clip_negative_factors : bool
        Whether to clip negative rate factors
    max_steps : int
        Max Gillespie steps per run
    
    Returns
    -------
    np.ndarray
        Shape (n_runs, len(times_obs)) — trajectories sampled on grid
    """
    times_obs = np.asarray(times_obs, dtype=float)
    rng = np.random.default_rng(int(seed))
    
    runs = []
    for _ in range(int(n_runs)):
        path = gillespie_population_one_path(
            R01=float(R01),
            k10=float(k10),
            k12=float(k12),
            k21=float(k21),
            times_obs=times_obs,
            rng=rng,
            ns=float(ns),
            ca1=float(ca1),
            ca2=float(ca2),
            cd1=float(cd1),
            cd2=float(cd2),
            n1_init=int(n1_init),
            n2_init=int(n2_init),
            clip_negative_factors=bool(clip_negative_factors),
            max_steps=int(max_steps),
        )
        runs.append(path)
    
    return np.asarray(runs, dtype=float)


def run_gillespie_population_ensemble(
    R01: float,
    k10: float,
    k12: float,
    k21: float,
    times_obs: np.ndarray,
    n_runs: int = 10,
    seed: int = 0,
    ns: float = 100.0,
    ca1: float = 0.0,
    ca2: float = 0.0,
    cd1: float = 0.0,
    cd2: float = 0.0,
    n1_init: int = 1,
    n2_init: int = 0,
    clip_negative_factors: bool = True,
    max_steps: int = 200_000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Run ensemble and return (mean, std, individual runs).
    Useful for summary statistics.
    
    Returns
    -------
    mean : np.ndarray
        Mean trajectory
    std : np.ndarray
        Standard deviation across runs
    runs : np.ndarray
        Individual trajectories
    """
    runs = run_gillespie_population_on_grid(
        R01=R01,
        k10=k10,
        k12=k12,
        k21=k21,
        times_obs=times_obs,
        n_runs=n_runs,
        seed=seed,
        ns=ns,
        ca1=ca1,
        ca2=ca2,
        cd1=cd1,
        cd2=cd2,
        n1_init=n1_init,
        n2_init=n2_init,
        clip_negative_factors=clip_negative_factors,
        max_steps=max_steps,
    )
    
    mean = runs.mean(axis=0)
    std = runs.std(axis=0, ddof=1) if runs.shape[0] > 1 else np.zeros_like(mean)
    
    return mean, std, runs
