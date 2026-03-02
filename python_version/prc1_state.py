# TODO: implement double attach function and rate calculations for doubly attaching
# precompute attachment rates for different offsets?
# make sure to check that there is an available spot when attaching

import numpy as np
from sortedcontainers import SortedList, SortedSet
from prc1 import Prc1
from copy import deepcopy

# probably not needed, but I currently use it to define top and bottom taken/untaken sites
class SortedSetAndComplement:
    """
    helper class to maintain a set and its complement
    """
    def __init__(self, universal_set):
        """
        Initialized to empty, with complement as the full universal set

        :param universal_set: iterable of all possible values in the set
        """
        self.complement = SortedSet(universal_set)
        self.current_set = SortedSet()

    def __contains__(self, value):
        return value in self.current_set
    
    def __len__(self):
        return len(self.current_set)
    
    def add(self, value):
        self.current_set.add(value)
        self.complement.remove(value)
    
    def remove(self, value):
        self.current_set.remove(value)
        self.complement.add(value)

class State:
    """
    maintains sorted list of doubly attached PRC1
    handles updating list when prc1 bind/unbind
    updates prc1 rates based on neighbors
    """
    # initialization
    def __init__(self, microtubule_length, site_spacing, microtubule_offset, spring_constant,
                 rest_length, k_B_T, microtubule_separation, singly_bound_detachment_rate, k0):
        # rates
        self.singly_bound_detachment_rate = singly_bound_detachment_rate

        # microtubule and site parameters
        self.microtubule_length = microtubule_length
        self.microtubule_separation = microtubule_separation
        
        # horizontal displacement of the top microtubule's left end from the bottom microtubule's left end
        # (if positive, top is further right than bottom)
        self.microtubule_offset = microtubule_offset

        # constants
        self.spring_constant = spring_constant
        self.rest_length = rest_length
        self.k_B_T = k_B_T
        self.k0 = k0

        # site parameters
        self.site_spacing = site_spacing
        self.num_sites = int(microtubule_length / site_spacing) + 1

        # keep track of taken (and untaken) sites for fast lookup
        self.__top_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.__bottom_taken_sites = SortedSetAndComplement(range(self.num_sites))

        # self.prc1_set = SortedSet()
        self.doubly_attached_prc1 = SortedSet()
        self.top_attached_prc1 = SortedSet()
        self.bottom_attached_prc1 = SortedSet()

        self.precompute_rates()

        # for debugging
        self.last_reaction = None
        self.last_reaction_prc1 = None
    
    # taken and untaken sites properties
    @property
    def top_untaken_sites(self):
        return self.__top_taken_sites.complement
    
    @property
    def bottom_untaken_sites(self):
        return self.__bottom_taken_sites.complement
    
    @property
    def top_taken_sites(self):
        return self.__top_taken_sites.current_set
    
    @property
    def bottom_taken_sites(self):
        return self.__bottom_taken_sites.current_set
    


    # BASIC OPERATIONS

    def clear(self):
        print("here")
        self.__top_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.__bottom_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.top_attached_prc1.clear()
        self.bottom_attached_prc1.clear()
        self.doubly_attached_prc1.clear()
        
    
    @property
    def num_prc1(self):
        return (len(self.doubly_attached_prc1) +
                len(self.top_attached_prc1) +
                len(self.bottom_attached_prc1))

    def get_prc1(self, index):
        if index < len(self.doubly_attached_prc1):
            return self.doubly_attached_prc1[index]
        index -= len(self.doubly_attached_prc1)
        if index < len(self.top_attached_prc1):
            return self.top_attached_prc1[index]
        index -= len(self.top_attached_prc1)
        if index < len(self.bottom_attached_prc1):
            return self.bottom_attached_prc1[index]
        
    def __getitem__(self, index):
        return self.get_prc1(index)

    def __len__(self):
        return (len(self.doubly_attached_prc1) 
                + len(self.top_attached_prc1) 
                + len(self.bottom_attached_prc1))
        
    def __iter__(self):
        for prc1 in self.doubly_attached_prc1:
            yield prc1
        for prc1 in self.top_attached_prc1:
            yield prc1
        for prc1 in self.bottom_attached_prc1:
            yield prc1
    


    # STATE MODIFICATION FUNCTIONS
    # function checklist:
    # - updates closest (doubly attached) neighbors

    def single_attach_prc1(self):
        prc1 = Prc1(self)
        while True:
            site = np.random.randint(0, self.num_sites)
            attachment_head_bottom = np.random.choice([True, False])
            if attachment_head_bottom:
                if site in self.bottom_taken_sites:
                    continue
                prc1.binding_site_bottom = site
                break
            else:
                if site in self.top_taken_sites:
                    continue
                prc1.binding_site_top = site
                break
        
        # update neighbors
        prc1.set_closest_neighbors()
        self.last_reaction = "single attach"
        self.last_reaction_prc1 = str(prc1)


    def double_attach_prc1(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        # closest_index = prc1.closest_index_on_other_side

        precomputed_rates, (left_bound, right_bound) = prc1.get_rates_and_range()
        self.last_reaction = "double attach"
        self.last_reaction_prc1 = str(prc1)
        # if len(self) > 2:
        #     print("attachment range: ", (left_bound, right_bound))
        #     print(self, "\n")
        if len(precomputed_rates) == 0:
            raise RuntimeError("tried to double attach prc1 with 0 attachment rate", "attachment rates:", precomputed_rates)
        total_rate = precomputed_rates[-1]
        random_value = np.random.uniform(0, total_rate)
        relative_attachment_index = np.searchsorted(precomputed_rates, random_value)
        attachment_index = relative_attachment_index + left_bound

        # cur_prc1_index = self.doubly_attached_prc1.index(prc1)
        if attachment_index <= prc1.left_neighbor_opposite_index:
            raise RuntimeError("LEFT BOUND WRONG", self.last_reaction_prc1, self)
        elif right_bound > prc1.right_neighbor_opposite_index:
            raise RuntimeError("RIGHT BOUND WRONG", self.last_reaction_prc1, self)

        if prc1.bottom_head_is_attached:
            if attachment_index in self.top_taken_sites:
                return
            prc1.binding_site_top = attachment_index

        elif prc1.top_head_is_attached:
            if attachment_index in self.bottom_taken_sites:
                return
            prc1.binding_site_bottom = attachment_index

        else:
            print("warning: tried to double attached non-attached prc1")

        # update doubly attached neighbors
        left_neighbor = prc1.closest_neighbor_left
        right_neighbor = prc1.closest_neighbor_right

        if left_neighbor is not None:
            left_neighbor.closest_neighbor_right = prc1
        if right_neighbor is not None:
            right_neighbor.closest_neighbor_left = prc1
        
        # other neighbors
        self.set_neighbors_between_prc1(left_neighbor, prc1)
        self.set_neighbors_between_prc1(prc1, right_neighbor)

            
    def detach_prc1(self, prc1_index):


        prc1 = self.get_prc1(prc1_index)

        # if singly attached, detach the only head
        if prc1.is_singly_attached:
            self.last_reaction = "single detach"
            self.last_reaction_prc1 = str(prc1)
            if prc1.bottom_head_is_attached:
                # detach from bottom
                prc1.binding_site_bottom = None
            else:
                # detach from top
                prc1.binding_site_top = None

        # if doubly attached, randomly detach one head
        else:
            self.last_reaction = "double detach"
            self.last_reaction_prc1 = str(prc1)
            detach_bottom_head = np.random.choice([True, False])
            if detach_bottom_head:
                # detach from bottom
                prc1.binding_site_bottom = None
            else:
                # detach from top
                prc1.binding_site_top = None

            left_neighbor = prc1.closest_neighbor_left
            right_neighbor = prc1.closest_neighbor_right

            # update doubly attached neighbors
            if left_neighbor is not None:
                left_neighbor.closest_neighbor_right = right_neighbor
            if right_neighbor is not None:
                right_neighbor.closest_neighbor_left = left_neighbor
                
            # update singly attached neighbors
            self.set_neighbors_between_prc1(left_neighbor, right_neighbor)

    # useful function for updating neighbors for attachment/detachment functions
    # sets the neighbors of all prc1 between left_prc1 and right_prc1 to left_prc1 and right_prc1
    def set_neighbors_between_prc1(self, left_prc1, right_prc1):
        # helper function to avoid code duplication :)
        if left_prc1 is not None and left_prc1.is_unattached:
            raise RuntimeError("LEFT NEIGHBOR UNATTACHED")
        if right_prc1 is not None and right_prc1.is_unattached:
            raise RuntimeError("RIGHT NEIGHBOR UNATTACHED")
        def set_neighbors_for_set(prc1_set):
            left_index = 0 if left_prc1 is None else prc1_set.bisect_left(left_prc1)
            right_index = len(prc1_set) if right_prc1 is None else prc1_set.bisect_right(right_prc1)
            for prc1 in prc1_set[left_index: right_index]:
                prc1.closest_neighbor_left = left_prc1
                prc1.closest_neighbor_right = right_prc1

        set_neighbors_for_set(self.bottom_attached_prc1)
        set_neighbors_for_set(self.top_attached_prc1)


    # PRECOMPUTED ATTACHMENT RATES

    # definitely don't have to compute for this many sites
    def precompute_rates(self, divisions=1):
        num_sites_to_compute = self.num_sites*2 # should probably be a multiple of 2
        self.__precomputed_division_size = self.site_spacing / divisions
        self.__precomputed_division_num = divisions
        offset = np.linspace(0, self.site_spacing, divisions, endpoint=False)
        base_horizontal_distances = np.arange(-(num_sites_to_compute//2), (num_sites_to_compute//2)+1) * self.site_spacing
        horizontal_distances = base_horizontal_distances[:, np.newaxis] + offset[np.newaxis, :] # index (site num, offset num)
        distances = np.sqrt(horizontal_distances**2 + self.microtubule_separation**2)
        E = 0.5 * self.spring_constant * np.maximum(distances - self.rest_length, 0)**2
        self.__precomputed_rates = self.k0 * np.exp(-.5 * E / self.k_B_T)  # index (site num, offset num)
        self.__precomputed_cumulative_rates = np.cumsum(self.__precomputed_rates, axis=0)
    
    @property
    def precomputed_cumulative_rates(self):
        """returns precomputed cumulative rates from left to right. 0 offset is included in left rates."""
        # offset = self.microtubule_offset % self.site_spacing
        # print(offset, self.site_spacing, self.__precomputed_division_size)        
        index = int(self.microtubule_offset // self.__precomputed_division_size) % self.__precomputed_division_num
        return self.__precomputed_cumulative_rates[:, index]

    # PRINTING
    def __repr__(self):
        return (
                f"State(\n"
                f"doubly_attached_prc1 = {self.doubly_attached_prc1}\n"
                f"top_attached_prc1 = {self.top_attached_prc1}\n"
                f"bottom_attached_prc1 = {self.bottom_attached_prc1}\n"
                f"top_taken_sites = {self.top_taken_sites}\n"
                f"bottom_taken_sites = {self.bottom_taken_sites}\n"
                f")"
                )
