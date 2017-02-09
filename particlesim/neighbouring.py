#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen, Simon Treu
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
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# imports
import numpy as np
import itertools as it
from particlesim.k_cython import fast_distances

class Neighbouring(object):
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
        return []

    def get_particles_within_radius(self, particle_id):
        pass


class NeighbouringPrimitiveLists(Neighbouring):
    def __init__(self, particle_positions, radius, box_size=5):
        super(NeighbouringPrimitiveLists, self).__init__(particle_positions, radius)
        self.box_size = box_size
        self.create_neighbourlist()

    # private methods

    # public methods
    def create_neighbourlist(self):
        n, r, pos, box_size   = self.n, self.r, self.particle_positions, self.box_size
        nlist       = [[] for i in range(n)]
        for i in range(n):
            for j in range(n):
                periodic_distance = np.linalg.norm(0.5 * box_size- (pos[i] - pos[j] + 0.5 * box_size) % box_size)
                if periodic_distance>=r or i==j: continue
                nlist[i].append(j)
        self._neighbourlist = nlist

    def get_particles_within_radius(self, particle_id):
        return self._neighbourlist[particle_id]  # returns the indices of the points
        # return [self.particle_positions[i] for i in self._neighbourlist[particle_id]] # alternative: return the points themselves

class NeighbouringCellLinkedLists(Neighbouring):
    def __init__(self, particle_positions, radius, box_size=1.0):
        super(NeighbouringCellLinkedLists, self).__init__(particle_positions, radius)
        self.box_size = float(box_size)
        self._cell_len = -1
        self.nr_cells = int(max(1, self.box_size/self.r))
        self._cell_len =  self.r + (self.box_size % self.r )/self.nr_cells

        self.create_neighbourlist()


    # private methods

    # public methods
    def create_neighbourlist(self): # in O(self.n)
        pos = self.particle_positions
        cell_linked_list = [[[[] for i in range(self.nr_cells)]for j in range(self.nr_cells)]for k in range(self.nr_cells)]
        # print("cll shape: ", len(cell_linked_list), len(cell_linked_list[0]), len(cell_linked_list[0][0])) #TODO no print in the end
        for i in range(self.n):
            x, y, z = (pos[i]/self._cell_len).astype(int) # // ist ganzzahlige division (ohne rest)
            #print ("i:", i, ", xyz: ", x,y,z, ", pos[i]:", pos[i]) #TODO no print in the end
            cell_linked_list[x][y][z].append(i) # we need only indices
        self._neighbourlist = cell_linked_list

    def get_particles_within_radius(self, particle_id):
        n, r, pos, cell_ll, cell_len = self.n, self.r, self.particle_positions, self._neighbourlist, self._cell_len
        box_size = self.box_size
        nr_cells = len(cell_ll)
        ret_idx, ret_dist = [], []

        p = pos[particle_id]
        cell = (p / cell_len).astype("int")
        #cell_dir = np.rint((p % r)/r).astype("int")
        #cell_dir[cell_dir == 0] = -1  # setting direction to nearest cell in xyz direction
        cell_dir = np.array(list(it.product([-1,0,1],repeat=3)))
        #cells = np.array([(cell + [x, y, z]) % nr_cells for x in [cell_dir[0], 0] for y in [cell_dir[1], 0] for z in [cell_dir[2], 0]])
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

        cell_linked_list = [[[[] for i in range(self.nr_cells)]for j in range(self.nr_cells)]for k in range(self.nr_cells)]
        for i in range(self.n):
            x, y, z = (new_positions[i]/self._cell_len).astype(int) # // ist ganzzahlige division (ohne rest)
            #print ("i:", i, ", xyz: ", x,y,z, ", pos[i]:", pos[i]) #TODO no print in the end
            cell_linked_list[x][y][z].append(i) # we need only indices

        self._neighbourlist = cell_linked_list


class NeighbouringCellLinkedListsArray(Neighbouring):
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
        self.head[:] = -1
        self.cell_ll[:] = -1

        for i in range(self.n):
            cell_index = self._calc_cell_index(i, xyz)
            old_head = self.head[cell_index]
            self.head[cell_index] = i
            self.cell_ll[i] = old_head

    def _get_particles_in_cell(self,cell_index):
        head_tmp = self.head[cell_index]
        particles_in_cell_idxs = []
        while(head_tmp != -1):
            particles_in_cell_idxs.append(head_tmp)
            head_tmp = self.cell_ll[head_tmp]
        return particles_in_cell_idxs



    def get_particles_within_radius(self, particle_id):

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
        x, y, z = (xyz[i] / self._cell_len).astype(int)
        return self._recalc_cell_index(x,y,z)

    def _recalc_cell_index(self, cell_x, cell_y, cell_z):
        return cell_x*1 + cell_y * self.nr_cells_one_d + cell_z * self.nr_cells_one_d**2

#n = 20000
#distances = np.zeros((n,n), dtype="float")
#box_len = 5.0
#xyz = np.arange(3*n, dtype="float").reshape(n,3) % box_len

#print("start")
#fast_distances(xyz, box_len, distances)
#print("done")
#print(xyz)
#print(distances)