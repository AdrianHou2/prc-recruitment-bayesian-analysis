import numpy as np
from sortedcontainers import SortedSet
from prc1 import Prc1


class SortedSetAndComplement(SortedSet):
    """
    Maintain a sorted set, its complement, and the set of untaken sites
    adjacent to a taken site.
    """
    def __init__(self, universal_set):
        self.universal_set = tuple(universal_set)
        super().__init__()
        self.complement = SortedSet(self.universal_set)
        self.adjacencies = SortedSet()

    def _refresh_local_adjacency(self, value):
        for x in (value - 1, value, value + 1):
            if x not in self.complement:
                self.adjacencies.discard(x)
                continue
            if (x - 1 in self) or (x + 1 in self):
                self.adjacencies.add(x)
            else:
                self.adjacencies.discard(x)

    def add(self, value):
        if value not in self.complement:
            return
        SortedSet.add(self, value)
        self.complement.discard(value)
        self._refresh_local_adjacency(value)

    def discard(self, value):
        if value not in self:
            return
        SortedSet.discard(self, value)
        self.complement.add(value)
        self._refresh_local_adjacency(value)

    def remove(self, value):
        if value not in self:
            raise KeyError(value)
        self.discard(value)


class State:
    """
    Maintains sorted lists of PRC1 objects and occupancy sets.
    Handles updates when PRC1 bind, unbind, or hop.
    """
    def __init__(self, microtubule_length, site_spacing, microtubule_offset, spring_constant,
                 rest_length, k_B_T, microtubule_separation, initial_binding_rate_per_site,
                 singly_bound_detachment_rate, base_double_attachment_rate,
                 base_double_detachment_rate, base_hopping_rate,
                 cooperativity_energy=0,
                 enable_initial_attach_cooperativity=False,
                 enable_second_head_attach_cooperativity=False,
                 enable_second_head_detach_cooperativity=False):

        # fit params
        self.singly_bound_detachment_rate = singly_bound_detachment_rate
        self.base_double_attachment_rate = base_double_attachment_rate
        self.base_double_detachment_rate = base_double_detachment_rate
        self.cooperativity_energy = cooperativity_energy

        # mode flags
        self.enable_initial_attach_cooperativity = enable_initial_attach_cooperativity
        self.enable_second_head_attach_cooperativity = enable_second_head_attach_cooperativity
        self.enable_second_head_detach_cooperativity = enable_second_head_detach_cooperativity

        # known base rates
        self.base_hopping_rate = base_hopping_rate
        self.initial_binding_rate_per_site = initial_binding_rate_per_site

        # microtubule and site parameters
        self.microtubule_length = microtubule_length
        self.microtubule_separation = microtubule_separation
        self.microtubule_offset = microtubule_offset

        # constants
        self.spring_constant = spring_constant
        self.rest_length = rest_length
        self.k_B_T = k_B_T

        # site parameters
        self.site_spacing = site_spacing
        self.num_sites = int(microtubule_length / site_spacing) + 1

        # occupancy sets
        self.top_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.bottom_taken_sites = SortedSetAndComplement(range(self.num_sites))

        # PRC1 collections
        self.doubly_attached_prc1 = SortedSet()
        self.top_attached_prc1 = SortedSet()
        self.bottom_attached_prc1 = SortedSet()

        self.precompute_rates()

        # debugging info
        self.last_reaction = None
        self.last_reaction_prc1 = None

        # turn on only while debugging
        self.debug_validate = False

    # ------------------------------------------------------------------
    # Initial attachment rates
    # ------------------------------------------------------------------

    @property
    def initial_binding_rate(self):
        return float(np.sum(self.cooperative_initial_binding_rates) +
                     np.sum(self.uncooperative_initial_binding_rates))

    @property
    def cooperative_initial_binding_rates(self):
        """
        (bottom binding rate, top binding rate) based on valid adjacent untaken sites.
        """
        if not self.enable_initial_attach_cooperativity:
            return np.array([0.0, 0.0], dtype=float)

        cooperativity_rate = self.initial_binding_rate_per_site * np.exp(
            0.5 * self.cooperativity_energy / self.k_B_T
        )
        return np.array([
            len(self.bottom_adjacent_sites),
            len(self.top_adjacent_sites),
        ], dtype=float) * cooperativity_rate

    @property
    def uncooperative_initial_binding_rates(self):
        """
        (bottom binding rate, top binding rate) based on valid non-adjacent untaken sites.
        """
        n_bottom_untaken = len(self.bottom_untaken_sites)
        n_top_untaken = len(self.top_untaken_sites)

        if not self.enable_initial_attach_cooperativity:
            return np.array([n_bottom_untaken, n_top_untaken], dtype=float) * self.initial_binding_rate_per_site

        n_bottom_uncoop = n_bottom_untaken - len(self.bottom_adjacent_sites)
        n_top_uncoop = n_top_untaken - len(self.top_adjacent_sites)

        return np.array([n_bottom_uncoop, n_top_uncoop], dtype=float) * self.initial_binding_rate_per_site

    # ------------------------------------------------------------------
    # Taken / untaken site helpers
    # ------------------------------------------------------------------

    @property
    def top_untaken_sites(self):
        return self.top_taken_sites.complement

    @property
    def bottom_untaken_sites(self):
        return self.bottom_taken_sites.complement

    @property
    def top_adjacent_sites(self):
        return self.top_taken_sites.adjacencies

    @property
    def bottom_adjacent_sites(self):
        return self.bottom_taken_sites.adjacencies

    def _bottom_nonadjacent_untaken_sites(self):
        return [s for s in self.bottom_untaken_sites if s not in self.bottom_adjacent_sites]

    def _top_nonadjacent_untaken_sites(self):
        return [s for s in self.top_untaken_sites if s not in self.top_adjacent_sites]

    # ------------------------------------------------------------------
    # Basic operations
    # ------------------------------------------------------------------

    def clear(self):
        self.top_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.bottom_taken_sites = SortedSetAndComplement(range(self.num_sites))
        self.top_attached_prc1.clear()
        self.bottom_attached_prc1.clear()
        self.doubly_attached_prc1.clear()

    @property
    def num_prc1(self):
        return (
            len(self.doubly_attached_prc1)
            + len(self.top_attached_prc1)
            + len(self.bottom_attached_prc1)
        )

    def get_prc1(self, index):
        if index < len(self.doubly_attached_prc1):
            return self.doubly_attached_prc1[index]
        index -= len(self.doubly_attached_prc1)

        if index < len(self.top_attached_prc1):
            return self.top_attached_prc1[index]
        index -= len(self.top_attached_prc1)

        if index < len(self.bottom_attached_prc1):
            return self.bottom_attached_prc1[index]

        raise IndexError("PRC1 index out of range")

    def __getitem__(self, index):
        return self.get_prc1(index)

    def __len__(self):
        return (
            len(self.doubly_attached_prc1)
            + len(self.top_attached_prc1)
            + len(self.bottom_attached_prc1)
        )

    def __iter__(self):
        for prc1 in self.doubly_attached_prc1:
            yield prc1
        for prc1 in self.top_attached_prc1:
            yield prc1
        for prc1 in self.bottom_attached_prc1:
            yield prc1

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_state(self):
        top_sites = []
        bottom_sites = []

        for prc1 in self.doubly_attached_prc1:
            if prc1.binding_site_top is None or prc1.binding_site_bottom is None:
                raise RuntimeError(f"Malformed doubly attached PRC1: {prc1}")
            top_sites.append(prc1.binding_site_top)
            bottom_sites.append(prc1.binding_site_bottom)

        for prc1 in self.top_attached_prc1:
            if prc1.binding_site_top is None or prc1.binding_site_bottom is not None:
                raise RuntimeError(f"Malformed top-attached PRC1: {prc1}")
            top_sites.append(prc1.binding_site_top)

        for prc1 in self.bottom_attached_prc1:
            if prc1.binding_site_bottom is None or prc1.binding_site_top is not None:
                raise RuntimeError(f"Malformed bottom-attached PRC1: {prc1}")
            bottom_sites.append(prc1.binding_site_bottom)

        if len(top_sites) != len(set(top_sites)):
            raise RuntimeError(f"Duplicate top occupancy: {top_sites}")
        if len(bottom_sites) != len(set(bottom_sites)):
            raise RuntimeError(f"Duplicate bottom occupancy: {bottom_sites}")

        if set(top_sites) != set(self.top_taken_sites):
            raise RuntimeError("Mismatch between PRC1 top sites and top_taken_sites")
        if set(bottom_sites) != set(self.bottom_taken_sites):
            raise RuntimeError("Mismatch between PRC1 bottom sites and bottom_taken_sites")

    def _maybe_validate(self):
        if self.debug_validate:
            self.validate_state()

    # ------------------------------------------------------------------
    # State modification functions
    # ------------------------------------------------------------------

    def single_attach_prc1(self):
        if not self.enable_initial_attach_cooperativity:
            bottom_candidates = list(self.bottom_untaken_sites)
            top_candidates = list(self.top_untaken_sites)

            bottom_rate = len(bottom_candidates) * self.initial_binding_rate_per_site
            top_rate = len(top_candidates) * self.initial_binding_rate_per_site
            total_rate = bottom_rate + top_rate

            if total_rate <= 0:
                raise RuntimeError("single_attach_prc1 called with no available sites")

            probabilities = np.array([bottom_rate, top_rate], dtype=float) / total_rate
            attachment_head_bottom = bool(np.random.choice([True, False], p=probabilities))

            if attachment_head_bottom:
                site = int(np.random.choice(bottom_candidates))
            else:
                site = int(np.random.choice(top_candidates))

        else:
            bottom_coop_candidates = list(self.bottom_adjacent_sites)
            top_coop_candidates = list(self.top_adjacent_sites)
            bottom_uncoop_candidates = self._bottom_nonadjacent_untaken_sites()
            top_uncoop_candidates = self._top_nonadjacent_untaken_sites()

            coop_rate_per_site = self.initial_binding_rate_per_site * np.exp(
                0.5 * self.cooperativity_energy / self.k_B_T
            )

            bottom_rate_coop = len(bottom_coop_candidates) * coop_rate_per_site
            top_rate_coop = len(top_coop_candidates) * coop_rate_per_site
            bottom_rate_uncoop = len(bottom_uncoop_candidates) * self.initial_binding_rate_per_site
            top_rate_uncoop = len(top_uncoop_candidates) * self.initial_binding_rate_per_site

            total_rate = (
                bottom_rate_coop + top_rate_coop
                + bottom_rate_uncoop + top_rate_uncoop
            )

            if total_rate <= 0:
                raise RuntimeError("single_attach_prc1 called with no valid sites")

            probabilities = np.array(
                [bottom_rate_coop, top_rate_coop, bottom_rate_uncoop, top_rate_uncoop],
                dtype=float
            ) / total_rate

            choice = int(np.random.choice([0, 1, 2, 3], p=probabilities))

            if choice == 0:
                attachment_head_bottom = True
                site = int(np.random.choice(bottom_coop_candidates))
            elif choice == 1:
                attachment_head_bottom = False
                site = int(np.random.choice(top_coop_candidates))
            elif choice == 2:
                attachment_head_bottom = True
                site = int(np.random.choice(bottom_uncoop_candidates))
            else:
                attachment_head_bottom = False
                site = int(np.random.choice(top_uncoop_candidates))

        prc1 = Prc1(self)
        if attachment_head_bottom:
            prc1.binding_site_bottom = site
        else:
            prc1.binding_site_top = site

        prc1.set_closest_neighbors()
        self.last_reaction = "single attach"
        self.last_reaction_prc1 = str(prc1)
        self._maybe_validate()

    def double_attach_prc1(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)

        precomputed_rates, (left_bound, right_bound) = prc1.get_rates_and_range()
        self.last_reaction = "double attach"
        self.last_reaction_prc1 = str(prc1)

        if len(precomputed_rates) == 0:
            raise RuntimeError("tried to double attach prc1 with 0 attachment rate")

        total_rate = precomputed_rates[-1]
        random_value = np.random.uniform(0, total_rate)
        relative_attachment_index = np.searchsorted(precomputed_rates, random_value)
        attachment_index = relative_attachment_index + left_bound

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
            raise RuntimeError("tried to double attach non-attached PRC1")

        left_neighbor = prc1.closest_neighbor_left
        right_neighbor = prc1.closest_neighbor_right

        if left_neighbor is not None:
            left_neighbor.closest_neighbor_right = prc1
        if right_neighbor is not None:
            right_neighbor.closest_neighbor_left = prc1

        self.set_neighbors_between_prc1(left_neighbor, prc1)
        self.set_neighbors_between_prc1(prc1, right_neighbor)
        self._maybe_validate()

    def detach_prc1(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)

        if prc1.is_singly_attached:
            self.last_reaction = "single detach"
            self.last_reaction_prc1 = str(prc1)

            if prc1.bottom_head_is_attached:
                prc1.binding_site_bottom = None
            else:
                prc1.binding_site_top = None

        else:
            self.last_reaction = "double detach"
            self.last_reaction_prc1 = str(prc1)

            bottom_rate = prc1.top_detachment_rate
            top_rate = prc1.bottom_detachment_rate
            probabilities = np.array([bottom_rate, top_rate], dtype=float) / (bottom_rate + top_rate)
            detach_bottom_head = np.random.choice([True, False], p=probabilities)

            if detach_bottom_head:
                prc1.binding_site_bottom = None
            else:
                prc1.binding_site_top = None

            left_neighbor = prc1.closest_neighbor_left
            right_neighbor = prc1.closest_neighbor_right

            if left_neighbor is not None:
                left_neighbor.closest_neighbor_right = right_neighbor
            if right_neighbor is not None:
                right_neighbor.closest_neighbor_left = left_neighbor

            self.set_neighbors_between_prc1(left_neighbor, right_neighbor)

        self._maybe_validate()

    def hop_top(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        if not prc1.top_head_is_attached:
            raise RuntimeError("tried to hop top head of PRC1 with unattached top head")

        self.last_reaction = "hop top"
        self.last_reaction_prc1 = str(prc1)

        hopping_rates = np.asarray(prc1.top_hopping_rates, dtype=float)
        total = hopping_rates.sum()
        if total <= 0:
            raise RuntimeError("hop_top called with zero total hopping rate")

        hopping_probabilities = hopping_rates / total
        step = int(np.random.choice([-1, 1], p=hopping_probabilities))
        new_site = prc1.binding_site_top + step

        if new_site not in self.top_untaken_sites:
            raise RuntimeError("tried to hop top head to invalid or occupied site")

        prc1.binding_site_top = new_site
        self._maybe_validate()

    def hop_bottom(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        if not prc1.bottom_head_is_attached:
            raise RuntimeError("tried to hop bottom head of PRC1 with unattached bottom head")

        self.last_reaction = "hop bottom"
        self.last_reaction_prc1 = str(prc1)

        hopping_rates = np.asarray(prc1.bottom_hopping_rates, dtype=float)
        total = hopping_rates.sum()
        if total <= 0:
            raise RuntimeError("hop_bottom called with zero total hopping rate")

        hopping_probabilities = hopping_rates / total
        step = int(np.random.choice([-1, 1], p=hopping_probabilities))
        new_site = prc1.binding_site_bottom + step

        if new_site not in self.bottom_untaken_sites:
            raise RuntimeError("tried to hop bottom head to invalid or occupied site")

        prc1.binding_site_bottom = new_site
        self._maybe_validate()

    def set_neighbors_between_prc1(self, left_prc1, right_prc1):
        if left_prc1 is not None and left_prc1.is_unattached:
            raise RuntimeError("LEFT NEIGHBOR UNATTACHED")
        if right_prc1 is not None and right_prc1.is_unattached:
            raise RuntimeError("RIGHT NEIGHBOR UNATTACHED")

        def set_neighbors_for_set(prc1_set):
            left_index = 0 if left_prc1 is None else prc1_set.bisect_left(left_prc1)
            right_index = len(prc1_set) if right_prc1 is None else prc1_set.bisect_right(right_prc1)
            for prc1 in prc1_set[left_index:right_index]:
                prc1.closest_neighbor_left = left_prc1
                prc1.closest_neighbor_right = right_prc1

        set_neighbors_for_set(self.bottom_attached_prc1)
        set_neighbors_for_set(self.top_attached_prc1)

    # ------------------------------------------------------------------
    # Precomputed attachment rates
    # ------------------------------------------------------------------

    def precompute_rates(self, divisions=1):
        num_sites_to_compute = self.num_sites * 2
        self.__precomputed_division_size = self.site_spacing / divisions
        self.__precomputed_division_num = divisions

        offset = np.linspace(0, self.site_spacing, divisions, endpoint=False)
        base_horizontal_distances = np.arange(
            -(num_sites_to_compute // 2),
            (num_sites_to_compute // 2) + 1
        ) * self.site_spacing

        horizontal_distances = base_horizontal_distances[:, np.newaxis] + offset[np.newaxis, :]
        distances = np.sqrt(horizontal_distances**2 + self.microtubule_separation**2)
        E = 0.5 * self.spring_constant * np.maximum(distances - self.rest_length, 0)**2

        self.__precomputed_rates = self.base_double_attachment_rate * np.exp(-0.5 * E / self.k_B_T)
        self.__precomputed_cumulative_rates = np.cumsum(self.__precomputed_rates, axis=0)

    @property
    def precomputed_cumulative_rates(self):
        index = int(self.microtubule_offset // self.__precomputed_division_size) % self.__precomputed_division_num
        return self.__precomputed_cumulative_rates[:, index]

    @property
    def precomputed_rates(self):
        index = int(self.microtubule_offset // self.__precomputed_division_size) % self.__precomputed_division_num
        return self.__precomputed_rates[:, index]

    # ------------------------------------------------------------------
    # Printing and misc
    # ------------------------------------------------------------------

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

    def get_distance_between_indices(self, bottom_index, top_index):
        offset = self.microtubule_offset
        site_spacing = self.site_spacing
        horizontal_displacement = (top_index * site_spacing + offset - bottom_index * site_spacing)
        vertical_displacement = self.microtubule_separation
        return np.sqrt(horizontal_displacement**2 + vertical_displacement**2)

    def get_energy_between_indices(self, bottom_index, top_index):
        if self.enable_second_head_detach_cooperativity:
            num_neighbors = [
                (top_index is not None) and (top_index - 1 in self.top_taken_sites),
                (top_index is not None) and (top_index + 1 in self.top_taken_sites),
                (bottom_index is not None) and (bottom_index - 1 in self.bottom_taken_sites),
                (bottom_index is not None) and (bottom_index + 1 in self.bottom_taken_sites),
            ].count(True)
            E = -self.cooperativity_energy * num_neighbors
        else:
            E = 0.0

        if top_index is None or bottom_index is None:
            return E

        distance = self.get_distance_between_indices(bottom_index, top_index)
        E += 0.5 * self.spring_constant * np.maximum(distance - self.rest_length, 0)**2
        return E