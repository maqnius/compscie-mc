#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen,
#                       Simon Treu, Hanna Wulkow, Thomas Hadler
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#

import numpy as np
import itertools as it

class Neighbouring(object):
    r"""
    Parameters
    ----------
    particle_positions : array (float)
        Positions of all particles inside the simulation-box
    radius : float
        Cutoff-radius for short-range coulomb interaction.
    """
    def __init__(self, particle_positions, radius):
        self.particle_positions = particle_positions
        self.n = len(particle_positions)
        self.r = float(radius)
        self._neighbourlist = None

    @property
    def particle_positions(self):
        return self._particle_positions

    @particle_positions.setter
    def particle_positions(self, value):
        # make some asserts
        if not type(value) == np.ndarray:
            raise TypeError('particle_positions must be an ndarray')
        self._particle_positions = value

    # private methods

    # public methods
    def create_neighbourlist(self):  # only internal
        r"""
        Returns
        -------
        list
            empty list
        """
        return []


    def get_particles_within_radius(self, particle_id):
        pass


class NeighbouringPrimitiveLists(Neighbouring):
    r"""
    Parameters
    ----------
    particle_positions : array (float)
        Positions of all particles inside the simulation-box
    radius : float
        Cutoff-radius for short-range coulomb interaction.
    box_size : int or float
        Side-length of the simulation-box.
    """
    def __init__(self, particle_positions, radius, box_size):
        super(NeighbouringPrimitiveLists, self).__init__(particle_positions, radius)
        self.box_size = box_size
        self.create_neighbourlist()

    # private methods

    # public methods
    def create_neighbourlist(self):
        r"""
        Creates a neighbour-list with a naive method.
        """
        n, r, pos, box_size   = self.n, self.r, self.particle_positions, self.box_size
        nlist       = [[] for i in range(n)]
        for i in range(n):
            for j in range(n):
                periodic_distance = np.linalg.norm(0.5 * box_size- (pos[i] - pos[j] + 0.5 * box_size) % box_size)
                if periodic_distance>=r or i==j: continue
                nlist[i].append(j)
        self._neighbourlist = nlist


    def get_particles_within_radius(self, particle_id):
        r"""
        Find all neighbouring particles within the cutoff-radius.

        Parameters
        ----------
        particle_id : int
            Particle-index

        Returns
        -------
        list (int)
            List of indices from neighbouring particles.
        """
        return self._neighbourlist[particle_id]  # returns the indices of the points


class NeighbouringCellLinkedLists(Neighbouring):
    r"""
    Parameters
    ----------
    particle_positions : array (float)
        Positions of all particles inside the simulation-box
    radius : float
        Cutoff-radius for short-range coulomb interaction.
    box_size : int or float
        Side-length of the simulation-box.
    """
    def __init__(self, particle_positions, radius, box_size):
        super(NeighbouringCellLinkedLists, self).__init__(particle_positions, radius)
        self.box_size = float(box_size)
        self._cell_len = -1
        self.nr_cells = int(max(1, self.box_size/self.r))
        self._cell_len =  self.r + (self.box_size % self.r )/self.nr_cells

        self.create_neighbourlist()


    # private methods

    # public methods
    def create_neighbourlist(self): # in O(self.n)
        """
        Creates a cell-linked neighbour-list.
        """
        pos = self.particle_positions
        cell_linked_list = [[[[] for i in range(self.nr_cells)]for j in range(self.nr_cells)]for k in range(self.nr_cells)]
        for i in range(self.n):
            x, y, z = (pos[i]/self._cell_len).astype(int) # // ist ganzzahlige division (ohne rest)
            cell_linked_list[x][y][z].append(i) # we need only indices
        self._neighbourlist = cell_linked_list


    def get_particles_within_radius(self, particle_id):
        r"""
        Find all neighbouring particles within the cutoff-radius and calculate their distances.

        Parameters
        ----------
        particle_id : int
            Particle-index

        Returns
        -------
        ret_idx : list (int)
            Indices of all neighbouring particles within the cutoff-radius.
        ret_dist : list (float)
            Distances from one particle to all neighbouring particles within the cutoff-radius.

        """
        n, r, pos, cell_ll, cell_len = self.n, self.r, self.particle_positions, self._neighbourlist, self._cell_len
        box_size = self.box_size
        nr_cells = len(cell_ll)
        ret_idx, ret_dist = [], []

        p = pos[particle_id]
        cell = (p / cell_len).astype("int")
        cell_dir = np.array(list(it.product([-1,0,1],repeat=3)))
        cells = np.array([(cell + cell_dir[i])%nr_cells for i in range(27)] )

        for cell_idx in cells:
            idx_x, idx_y, idx_z = cell_idx
            cell_list = cell_ll[idx_x][idx_y][idx_z]
            for neigh_idx in cell_list:
                neigh = pos[neigh_idx]
                periodic_distance = np.linalg.norm(0.5 * box_size- (p - neigh + 0.5 * box_size) % box_size)
                if periodic_distance>=r or particle_id==neigh_idx: continue
                ret_idx.append(neigh_idx)
                ret_dist.append(periodic_distance)
        return ret_idx, ret_dist


    def update_cells(self, new_positions):
        r"""
        Updates the cell-linked list.

        Parameters
        ----------
        new_positions : array (float)
            New particle positions.

        """
        cell_linked_list = [[[[] for i in range(self.nr_cells)]for j in range(self.nr_cells)]for k in range(self.nr_cells)]
        for i in range(self.n):
            x, y, z = (new_positions[i]/self._cell_len).astype(int) # // ist ganzzahlige division (ohne rest)
            cell_linked_list[x][y][z].append(i) # we need only indices

        self._neighbourlist = cell_linked_list


class NeighbouringCellLinkedListsArray(Neighbouring):
    r"""
    Parameters
    ----------
    particle_positions : array (float)
        Positions of all particles inside the simulation-box
    box_size : int or float
        Side-length of the simulation-box.
    radius : float
        Cutoff-radius for short-range coulomb interaction.

    """
    def __init__(self, particle_positions, box_size, radius):
        super(NeighbouringCellLinkedListsArray, self).__init__(particle_positions, radius)
        self.box_size = float(box_size)
        self.nr_cells_one_d = int(max(1, self.box_size / self.r))
        self.nr_cells = self.nr_cells_one_d**3
        self._cell_len =  self.r + (self.box_size % self.r )/self.nr_cells
        self.head = np.ones(self.nr_cells, dtype=int) * -1
        self.cell_ll = np.ones(self.n, dtype=int) * -1
        self.update_neighbourlist(self.particle_positions)


    def update_neighbourlist(self, xyz):
        r"""
        Updates the neighbour-list.

        Parameters
        ----------
        xyz : ndarray(n,3), float
            Position of n particles in x,y,z coordinates.

        """
        self.head[:] = -1
        self.cell_ll[:] = -1

        for i in range(self.n):
            cell_index = self._calc_cell_index(i, xyz)
            old_head = self.head[cell_index]
            self.head[cell_index] = i
            self.cell_ll[i] = old_head


    def _get_particles_in_cell(self,cell_index):
        r"""
        Finds all particles inside a cell.

        Parameters
        ----------
        cell_index : list (int)
            x, y, z indices of the cell

        Returns
        -------
        list (int)
            indices of all particles inside a cell

        """
        head_tmp = self.head[cell_index]
        particles_in_cell_idxs = []
        while(head_tmp != -1):
            particles_in_cell_idxs.append(head_tmp)
            head_tmp = self.cell_ll[head_tmp]
        return particles_in_cell_idxs


    def get_particles_within_radius(self, particle_id):
        r"""
        Find all neighbouring particles within the cutoff-radius and calculate their distances.

        Parameters
        ----------
        particle_id : int
            Particle-index

        Returns
        -------
        result_idx : list (int)
            Indices of all neighbouring particles within the cutoff-radius.
        result_dist : list (float)
            Distances from one particle to all neighbouring particles within the cutoff-radius.
        """
        result_idx, result_dist = [], []

        p = self.particle_positions[particle_id]
        cell = (p / self._cell_len).astype("int")
        cell_dir = np.array(list(it.product([-1,0,1],repeat=3)))
        cells = np.array([(cell + cell_dir[i])%self.nr_cells for i in range(27)] )

        for cell_idx in cells:
            idx_x, idx_y, idx_z = cell_idx
            cell_idx_one_dim = self._recalc_cell_index(idx_x, idx_y, idx_z)
            cell_list = self._get_particles_in_cell(cell_idx_one_dim)

            for neigh_idx in cell_list:
                neigh = self.particle_positions[neigh_idx]

                periodic_distance = np.linalg.norm(0.5 * self.box_size- (p - neigh + 0.5 * self.box_size) % self.box_size)
                if periodic_distance>=self.r or particle_id==neigh_idx: continue
                result_idx.append(neigh_idx)
                result_dist.append(periodic_distance)

        return result_idx, result_dist


    def _calc_cell_index(self, i, xyz):
        r"""
        Gets the cell indices of a particle.

        Parameters
        ----------
        i : int
            particle index
        xyz : ndarray(n,3), float
            Position of n particles in x,y,z coordinates.

        Returns
        -------
        int
            cell index of particle i
        """
        x, y, z = (xyz[i] / self._cell_len).astype(int)
        return self._recalc_cell_index(x,y,z)


    def _recalc_cell_index(self, cell_x, cell_y, cell_z):
        r"""
        Calculates the cell index from 3D call indices.

        Parameters
        ----------
        cell_x : int
            cell index in x direction
        cell_y : int
            cell index in y direction
        cell_z : int
            cell index in z direction

        Returns
        -------
        int
            total cell index

        """
        return cell_x*1 + cell_y * self.nr_cells_one_d + cell_z * self.nr_cells_one_d**2
