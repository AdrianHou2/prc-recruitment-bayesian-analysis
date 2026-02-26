import numpy as np

class Prc1:
    def __init__(self, state):
        self.state = state
        self.__binding_site_top = None
        self.__binding_site_bottom = None
        self.closest_neighbor_left = None
        self.closest_neighbor_right = None
        self.__prev_binding_sites = (-1, -1)
        self.__prev_detachment_rate = 0

    def update_prc1_attributes(self, new_bottom_index, new_top_index):
        """Removes, updates, and re-adds the molecule to maintain sorted order."""
        state = self.state
        if self.is_doubly_attached:
            state.doubly_attached_prc1.remove(self)
            state.top_taken_sites.discard(self.__binding_site_top)
            state.bottom_taken_sites.discard(self.__binding_site_bottom)
        elif self.top_head_is_attached:
            state.top_attached_prc1.remove(self)
            state.top_taken_sites.discard(self.__binding_site_top)
        elif self.bottom_head_is_attached:
            state.bottom_attached_prc1.remove(self)
            state.bottom_taken_sites.discard(self.__binding_site_bottom)

        self.__binding_site_top = new_top_index
        self.__binding_site_bottom = new_bottom_index
        
        if self.is_doubly_attached:
            state.top_taken_sites.add(self.binding_site_top)
            state.bottom_taken_sites.add(self.binding_site_bottom)
            state.doubly_attached_prc1.add(self)
        elif self.top_head_is_attached:
            state.top_taken_sites.add(self.binding_site_top)
            state.top_attached_prc1.add(self)
        elif self.bottom_head_is_attached:
            state.bottom_taken_sites.add(self.binding_site_bottom)
            state.bottom_attached_prc1.add(self)

    @property
    def binding_site_top(self): return self.__binding_site_top
    @binding_site_top.setter
    def binding_site_top(self, value): self.update_prc1_attributes(self.__binding_site_bottom, value)
    
    @property
    def binding_site_bottom(self): return self.__binding_site_bottom
    @binding_site_bottom.setter
    def binding_site_bottom(self, value): self.update_prc1_attributes(value, self.__binding_site_top)

    def set_closest_neighbors(self):
        doubly_attached_prc1 = self.state.doubly_attached_prc1
        if self in doubly_attached_prc1:
            idx = doubly_attached_prc1.index(self)
            self.closest_neighbor_left = doubly_attached_prc1[idx - 1] if idx > 0 else None
            self.closest_neighbor_right = doubly_attached_prc1[idx + 1] if idx < len(doubly_attached_prc1) - 1 else None
        else:
            idx = doubly_attached_prc1.bisect_left(self)
            self.closest_neighbor_left = doubly_attached_prc1[idx - 1] if idx > 0 else None
            self.closest_neighbor_right = doubly_attached_prc1[idx] if idx < len(doubly_attached_prc1) else None
    
    def __lt__(self, other):
        """Compare position with object ID as a tie-breaker for SortedSet uniqueness."""
        def get_pos_key(p):
            if p.bottom_head_is_attached:
                return (0, p.binding_site_bottom)
            return (1, p.binding_site_top * p.state.site_spacing + p.state.microtubule_offset)
        
        s_key, o_key = get_pos_key(self), get_pos_key(other)
        if s_key != o_key: return s_key < o_key
        return id(self) < id(other)

    @property
    def bottom_head_is_attached(self): return self.binding_site_bottom is not None
    @property
    def top_head_is_attached(self): return self.binding_site_top is not None
    @property
    def is_doubly_attached(self): return self.top_head_is_attached and self.bottom_head_is_attached
    @property
    def is_singly_attached(self): return self.top_head_is_attached ^ self.bottom_head_is_attached
    
    @property
    def distance_between_heads(self):
        return np.sqrt(self.horizontal_distance_between_heads**2 + self.vertical_distance_between_heads**2)
    
    @property
    def horizontal_distance_between_heads(self):
        top_pos = self.binding_site_top * self.state.site_spacing + self.state.microtubule_offset
        bottom_pos = self.binding_site_bottom * self.state.site_spacing
        return np.abs(top_pos - bottom_pos)
    
    @property
    def vertical_distance_between_heads(self): return self.state.microtubule_separation

    @property
    def closest_index_on_other_side(self):
        offset, spacing = self.state.microtubule_offset, self.state.site_spacing
        if self.top_head_is_attached: return np.floor(self.binding_site_top + offset/spacing).astype(int)
        elif self.bottom_head_is_attached: return np.ceil(self.binding_site_bottom - offset/spacing).astype(int)
        return None

    @property
    def left_neighbor_opposite_index(self):
        if self.closest_neighbor_left is None: return -1
        return self.closest_neighbor_left.binding_site_top if self.bottom_head_is_attached else self.closest_neighbor_left.binding_site_bottom

    @property
    def right_neighbor_opposite_index(self):
        if self.closest_neighbor_right is None: return self.state.num_sites
        return self.closest_neighbor_right.binding_site_top if self.bottom_head_is_attached else self.closest_neighbor_right.binding_site_bottom

    def get_rates_and_range(self):
        precomputed_rates = self.state.precomputed_cumulative_rates
        zero_index = len(precomputed_rates) // 2
        l_idx, r_idx = self.left_neighbor_opposite_index + 1, self.right_neighbor_opposite_index
        c_idx = self.closest_index_on_other_side
        l_range, r_range = np.clip([l_idx - c_idx + zero_index, r_idx - c_idx + zero_index], 0, len(precomputed_rates))
        actual_range = np.array([l_range, r_range]) - zero_index + c_idx
        extra_rate = precomputed_rates[l_range-1] if l_range > 0 else 0
        return precomputed_rates[l_range:r_range] - extra_rate, actual_range

    @property
    def double_attachment_rate(self):
        if self.is_doubly_attached: return 0
        rates, _ = self.get_rates_and_range()
        return rates[-1] if len(rates) > 0 else 0

    @property
    def detachment_rate(self):
        if self.is_doubly_attached:
            cur = (self.binding_site_bottom, self.binding_site_top)
            if cur == self.__prev_binding_sites: return self.__prev_detachment_rate
            self.__prev_binding_sites = cur
            E = 0.5 * self.state.spring_constant * np.maximum(self.distance_between_heads - self.state.rest_length, 0)**2
            self.__prev_detachment_rate = 2 * self.state.k0 * np.exp(.5 * E / self.state.k_B_T)
            return self.__prev_detachment_rate
        return self.state.singly_bound_detachment_rate

    def __repr__(self): return f"({self.binding_site_bottom}, {self.binding_site_top})"