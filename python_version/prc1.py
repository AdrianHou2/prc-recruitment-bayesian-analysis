import numpy as np

class Prc1:
    """""
    Basic data class for PRC1 with some properties for readibility.
    Represents a PRC1 molecule and stores the binding sites.
    Binding sites are the integer indices if bound or None if unbound.
    """
    def __init__(self, state):
        """
        initialize PRC1 molecule
        """
        self.state = state  # so it can access state constants without having to pass them around
        self.__binding_site_top = None
        self.__binding_site_bottom = None

        # closest neighbors that are doubly attached (for fast rate calculations)
        self.closest_neighbor_left = None
        self.closest_neighbor_right = None

        # for avoiding recomputing detachment rates
        self.__prev_binding_sites = (-1, -1)
        self.__prev_detachment_rate = 0

    # GENERAL FUNCTIONS

    def update_prc1_attributes(self, new_bottom_index, new_top_index):
        """
        updates prc1 with new indices,
        re-sorts prc1 in prc1_set and doubly_attached_prc1,
        updates top_taken_sites and bottom_taken_sites
        """
        state = self.state

        # temporarily remove prc1 from state
        if self.is_doubly_attached:
            state.top_taken_sites.discard(self.__binding_site_top)
            state.bottom_taken_sites.discard(self.__binding_site_bottom)
            state.doubly_attached_prc1.remove(self)
        else:
            if self.top_head_is_attached:
                state.top_attached_prc1.remove(self)
                state.top_taken_sites.discard(self.__binding_site_top)
            elif self.bottom_head_is_attached:
                state.bottom_attached_prc1.remove(self)
                state.bottom_taken_sites.discard(self.__binding_site_bottom)

        # update binding sites
        self.__binding_site_top = new_top_index
        self.__binding_site_bottom = new_bottom_index
        if self.top_head_is_attached: state.top_taken_sites.add(self.binding_site_top)
        if self.bottom_head_is_attached: state.bottom_taken_sites.add(self.binding_site_bottom)

        # add prc1 to sorted sets as necessary, also modify taken sites
        if self.is_doubly_attached:
            state.doubly_attached_prc1.add(self)
        elif self.is_singly_attached:
            if self.top_head_is_attached:
                state.top_attached_prc1.add(self)
                state.top_taken_sites.add(self.binding_site_top)
            elif self.bottom_head_is_attached:
                state.bottom_attached_prc1.add(self)
                state.bottom_taken_sites.add(self.binding_site_bottom)

    @property
    def binding_site_top(self):
        return self.__binding_site_top
    
    @binding_site_top.setter
    def binding_site_top(self, value):
        self.update_prc1_attributes(self.__binding_site_bottom, value)

    @property
    def binding_site_bottom(self):
        return self.__binding_site_bottom
    
    @binding_site_bottom.setter
    def binding_site_bottom(self, value):
        self.update_prc1_attributes(value, self.__binding_site_top)
        
    

    def set_closest_neighbors(self):
        doubly_attached_prc1 = self.state.doubly_attached_prc1
        index = doubly_attached_prc1.bisect_left(self)
        if index > 0:
            self.closest_neighbor_left = doubly_attached_prc1[index - 1]
        else:
            self.closest_neighbor_left = None
        if index < len(doubly_attached_prc1):
            self.closest_neighbor_right = doubly_attached_prc1[index]
        else:
            self.closest_neighbor_right = None
    
    # sorts by top or bottom if they are both attached on the same side,
    # otherwise sorts singly attached by their only attachment
    def __lt__ (self, other):
        if self.is_doubly_attached and other.is_doubly_attached:
            return self.binding_site_bottom < other.binding_site_bottom
        elif self.bottom_head_is_attached:
            if other.bottom_head_is_attached:
                return self.binding_site_bottom < other.binding_site_bottom
            else:
                return self.binding_site_bottom < other.binding_site_top
        elif self.top_head_is_attached:
            if other.top_head_is_attached:
                return self.binding_site_top < other.binding_site_top
            else:
                return self.binding_site_top < other.binding_site_bottom
        return True  # self is unbound, arbitrary order (shouldn't happen)
    # ATTACHMENT RELATED PROPERTIES

    @property
    def bottom_head_is_attached(self):
        return self.binding_site_bottom is not None
    
    @property
    def top_head_is_attached(self):
        return self.binding_site_top is not None
    
    @property
    def is_doubly_attached(self):
        return self.top_head_is_attached and self.bottom_head_is_attached
    
    @property
    def is_singly_attached(self):
        return not self.is_doubly_attached


    # DISTANCE PROPERTIES

    @property
    def distance_between_heads(self):
        return np.sqrt(self.horizontal_distance_between_heads**2 + self.vertical_distance_between_heads**2)
    
    @property
    def horizontal_distance_between_heads(self):
        """always positive"""
        top_index = self.binding_site_top
        bottom_index = self.binding_site_bottom
        offset = self.state.microtubule_offset
        site_spacing = self.state.site_spacing
        horizontal_displacement = (top_index * site_spacing + 
                               offset - 
                               bottom_index * site_spacing)
        return np.abs(horizontal_displacement)
    
    @property
    def vertical_distance_between_heads(self):
        return self.state.microtubule_separation
    
    @property
    def closest_index_on_other_side(self):
        """returns closest binding site index (where the bottom site is on the left) on the opposite microtubule"""
        offset = self.state.microtubule_offset
        site_spacing = self.state.site_spacing

        if self.is_doubly_attached:
            print("Warning: closest index on other side requested for doubly attached PRC1")
            return None
        
        elif self.top_head_is_attached:
            # # find closest bottom site
            closest_bottom_index = np.floor(self.binding_site_top + offset/site_spacing).astype(int)
            return closest_bottom_index
        
        elif self.bottom_head_is_attached:
            # find closest top site
            closest_top_index = np.ceil(self.binding_site_bottom - offset/site_spacing).astype(int)
            return closest_top_index
        
        else:
            print("Warning: closest index on other side requested for unbound PRC1")
            return None
        
    # NEIGHBOR RELATED PROPERTIES
    def __neighbor_opposite_index(self, neighbor, out_of_bounds):
        """finds the opposite head of neighbor, returns out_of_bounds if no neighbor"""
        if neighbor is None:
            return out_of_bounds
        elif self.bottom_head_is_attached:
            return neighbor.binding_site_top
        elif self.top_head_is_attached:
            return neighbor.binding_site_bottom
        else:
            raise ValueError("Neighbor index requested for unbound PRC1")
        
    @property
    def left_neighbor_opposite_index(self):
        return self.__neighbor_opposite_index(self.closest_neighbor_left, -1)
    
    @property
    def right_neighbor_opposite_index(self):
        return self.__neighbor_opposite_index(self.closest_neighbor_right, self.state.num_sites)
    

    # ATTACHMENT RATE PROPERTIES
    def get_rates_and_range(self):
        """
        returns an array of cumulative attachment rates, as well as a list of indices it corresponds to\n
        these rates do not take into account taken sites, so are the rates to attempt to attach to each site\n
        range is [left inclusive, right exclusive)\n
        """
        # get precomputed rates and the index that corresponds to current position
        precomputed_rates = self.state.precomputed_cumulative_rates
        zero_index = (len(precomputed_rates) // 2)

        # get the full attachment range to the other side (left inclusive, right exclusive)
        left_index = self.left_neighbor_opposite_index + 1
        right_index = self.right_neighbor_opposite_index
        attachment_range = np.array([left_index, right_index])

        # subtract closest_index to get the attachment_range relative to this prc1's position
        closest_index = self.closest_index_on_other_side
        relative_attachment_range = attachment_range-closest_index

        # make sure left_range, right_range are in bounds to access precomputed_rates
        left_range, right_range = relative_attachment_range + zero_index
        if left_range < 0: left_range = 0
        if right_range > len(precomputed_rates): right_range = len(precomputed_rates)
        actual_range = np.array([left_range, right_range]) - zero_index + closest_index

        # subtract off cumulative rates to the left of the interval
        # to get the actual cumulative rates for the interval
        extra_rate = precomputed_rates[left_range-1] if left_range > 0 else 0
        return precomputed_rates[left_range:right_range] - extra_rate, actual_range
    
    @property
    def double_attachment_rate(self):
        if self.is_doubly_attached:
            return 0
        
        cumulative_rates, _ = self.get_rates_and_range()
        if len(cumulative_rates) == 0:
            return 0
        else:
            return cumulative_rates[-1]
    
    @property
    def attachment_rate(self):
        if self.is_doubly_attached:
            return 0
        elif self.is_singly_attached:
            return self.double_attachment_rate
        else:
            print("Warning: PRC1 attachment rate requested for unbound PRC1")
            return 0  # unbound (shouldn't happen?)
        
    # DETACHMENT RATE PROPERTIES
    @property
    def detachment_rate(self):
        spring_constant = self.state.spring_constant
        k0 = self.state.k0
        k_B_T = self.state.k_B_T
        rest_length = self.state.rest_length

        if self.is_doubly_attached:
            # if already computed detachment rate for these sites, return that
            cur_binding_sites = (self.binding_site_bottom, self.binding_site_top)
            if cur_binding_sites == self.__prev_binding_sites:
                return self.__prev_detachment_rate
            
            # otherwise recompute detachment rate
            self.__prev_binding_sites = cur_binding_sites
            E = 0.5 * spring_constant * np.maximum(self.distance_between_heads - rest_length, 0)**2
            self.__prev_detachment_rate = 2 * k0 * np.exp(.5 * E / k_B_T) # factor of 2 since either head can detach

            return self.__prev_detachment_rate
        elif self.is_singly_attached:
            return self.state.singly_bound_detachment_rate
        else:
            print("Warning: PRC1 detachment rate requested for unbound PRC1")
            return 0  # unbound (shouldn't happen?)
        
    # PRINTING
    def __repr__(self):
        return f"({self.binding_site_bottom}, {self.binding_site_top})"