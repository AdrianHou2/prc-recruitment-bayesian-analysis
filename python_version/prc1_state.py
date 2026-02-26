import numpy as np
from sortedcontainers import SortedSet
from prc1 import Prc1

class SortedSetAndComplement:
    def __init__(self, universal_set):
        self.complement = SortedSet(universal_set)
        self.current_set = SortedSet()
    def add(self, v):
        self.current_set.add(v)
        self.complement.remove(v)
    def remove(self, v):
        self.current_set.remove(v)
        self.complement.add(v)
    def discard(self, v):
        if v in self.current_set: self.remove(v)

class State:
    def __init__(self, microtubule_length, site_spacing, microtubule_offset, spring_constant,
                 rest_length, k_B_T, microtubule_separation, singly_bound_detachment_rate, k0):
        self.singly_bound_detachment_rate = singly_bound_detachment_rate
        self.microtubule_length, self.microtubule_separation = microtubule_length, microtubule_separation
        self.microtubule_offset = microtubule_offset
        self.spring_constant, self.rest_length, self.k_B_T, self.k0 = spring_constant, rest_length, k_B_T, k0
        self.site_spacing = site_spacing
        self.num_sites = int(microtubule_length / site_spacing) + 1
        self.__top_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.__bottom_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.doubly_attached_prc1 = SortedSet()
        self.top_attached_prc1 = SortedSet()
        self.bottom_attached_prc1 = SortedSet()
        self.precompute_rates()
        self.last_reaction = self.last_reaction_prc1 = None
    
    @property
    def top_taken_sites(self): return self.__top_taken_sites.current_set
    @property
    def bottom_taken_sites(self): return self.__bottom_taken_sites.current_set

    def get_prc1(self, index):
        if index < len(self.doubly_attached_prc1): return self.doubly_attached_prc1[index]
        index -= len(self.doubly_attached_prc1)
        if index < len(self.top_attached_prc1): return self.top_attached_prc1[index]
        return self.bottom_attached_prc1[index - len(self.top_attached_prc1)]

    def __len__(self): return len(self.doubly_attached_prc1) + len(self.top_attached_prc1) + len(self.bottom_attached_prc1)
        
    def __iter__(self):
        for p in self.doubly_attached_prc1: yield p
        for p in self.top_attached_prc1: yield p
        for p in self.bottom_attached_prc1: yield p

    def single_attach_prc1(self):
        prc1 = Prc1(self)
        while True:
            site, bottom = np.random.randint(0, self.num_sites), np.random.choice([True, False])
            if bottom and site not in self.bottom_taken_sites:
                prc1.binding_site_bottom = site; break
            elif not bottom and site not in self.top_taken_sites:
                prc1.binding_site_top = site; break
        prc1.set_closest_neighbors()
        self.last_reaction, self.last_reaction_prc1 = "single attach", prc1

    def double_attach_prc1(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        rates, (l_bound, _) = prc1.get_rates_and_range()
        if len(rates) == 0: return
        target_site = np.searchsorted(rates, np.random.uniform(0, rates[-1])) + l_bound
        if prc1.bottom_head_is_attached:
            if target_site in self.top_taken_sites: return
            prc1.binding_site_top = target_site
        else:
            if target_site in self.bottom_taken_sites: return
            prc1.binding_site_bottom = target_site

        idx = self.doubly_attached_prc1.index(prc1)
        left = self.doubly_attached_prc1[idx - 1] if idx > 0 else None
        right = self.doubly_attached_prc1[idx + 1] if idx < len(self.doubly_attached_prc1) - 1 else None
        prc1.closest_neighbor_left, prc1.closest_neighbor_right = left, right
        if left: left.closest_neighbor_right = prc1
        if right: right.closest_neighbor_left = prc1
        self.set_neighbors_between_prc1(left, prc1)
        self.set_neighbors_between_prc1(prc1, right)
        self.last_reaction, self.last_reaction_prc1 = "double attach", prc1
            
    def detach_prc1(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        was_doubly, old_l, old_r = prc1.is_doubly_attached, prc1.closest_neighbor_left, prc1.closest_neighbor_right
        if prc1.is_singly_attached:
            if prc1.bottom_head_is_attached: prc1.binding_site_bottom = None
            else: prc1.binding_site_top = None
        else:
            if np.random.choice([True, False]): prc1.binding_site_bottom = None
            else: prc1.binding_site_top = None
        if was_doubly:
            if old_l: old_l.closest_neighbor_right = old_r
            if old_r: old_r.closest_neighbor_left = old_l
            self.set_neighbors_between_prc1(old_l, old_r)
            if prc1.is_singly_attached: prc1.set_closest_neighbors()
        self.last_reaction_prc1 = prc1

    def set_neighbors_between_prc1(self, left_p, right_p):
        def update_set(p_set):
            l_idx = 0 if left_p is None else p_set.bisect_right(left_p)
            r_idx = len(p_set) if right_p is None else p_set.bisect_left(right_p)
            for p in p_set[l_idx:r_idx]:
                p.closest_neighbor_left, p.closest_neighbor_right = left_p, right_p
        update_set(self.bottom_attached_prc1); update_set(self.top_attached_prc1)

    def precompute_rates(self, divisions=1):
        num = self.num_sites * 2
        self.__pre_div_size, self.__pre_div_num = self.site_spacing / divisions, divisions
        off = np.linspace(0, self.site_spacing, divisions, endpoint=False)
        h_dist = (np.arange(-(num//2), (num//2)+1) * self.site_spacing)[:, np.newaxis] + off[np.newaxis, :]
        E = 0.5 * self.spring_constant * np.maximum(np.sqrt(h_dist**2 + self.microtubule_separation**2) - self.rest_length, 0)**2
        self.__pre_cum_rates = np.cumsum(self.k0 * np.exp(-.5 * E / self.k_B_T), axis=0)
    
    @property
    def precomputed_cumulative_rates(self):
        return self.__pre_cum_rates[:, int(self.microtubule_offset // self.__pre_div_size) % self.__pre_div_num]