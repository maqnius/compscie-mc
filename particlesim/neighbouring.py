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


class Neighbouring(object):
    def __init__(self, particle_positions, radius=float("inf")):
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
    def create_neighbourlist(self, particle_positions, radius):  # only internal
        return []

    def get_particles_within_radius(self, particle_id):
        pass


class NeighbouringPrimitiveLists(Neighbouring):
    def __init__(self, particle_positions, radius=float("inf"), box_size=5):
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
    def __init__(self, particle_positions, radius=float("inf"), box_size=1.0):
        super(NeighbouringCellLinkedLists, self).__init__(particle_positions, radius)
        self.box_size = float(box_size)
        self._cell_len = -1
        self.create_neighbourlist()


    # private methods

    # public methods
    def create_neighbourlist(self): # in O(self.n)
        n, r, pos, box_size = self.n, self.r, self.particle_positions, self.box_size
        nr_cells = int(box_size / r)
        cell_len = r + (box_size % r )/nr_cells # evenly distributing the overlap
        self._cell_len = cell_len
        nr_cells = max(1, nr_cells) # catches case that nr_cells is 0
        cell_linked_list = [[[[] for i in range(nr_cells)]for j in range(nr_cells)]for k in range(nr_cells)]
        # print("cll shape: ", len(cell_linked_list), len(cell_linked_list[0]), len(cell_linked_list[0][0])) #TODO no print in the end
        for i in range(n):
            x, y, z = (pos[i]/cell_len).astype(int) # // ist ganzzahlige division (ohne rest)
            #print ("i:", i, ", xyz: ", x,y,z, ", pos[i]:", pos[i]) #TODO no print in the end
            cell_linked_list[x][y][z].append(i) # we need only indices
        self._neighbourlist = cell_linked_list

    def get_particles_within_radius(self, particle_id):
        n, r, pos, cell_ll, cell_len = self.n, self.r, self.particle_positions, self._neighbourlist, self._cell_len
        box_size = self.box_size
        nr_cells = len(cell_ll)
        ret = []

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
                ret.append(neigh_idx)
        return ret


if __name__=="__main__":
    box_size = float(7.6)
    nr_particles = 250
    for i in range(100):
        particle_pos = np.random.rand(nr_particles, 3)*box_size
        #print (particle_pos)

        NL_PL = NeighbouringPrimitiveLists(particle_pos, radius=1.2, box_size=box_size)
        NL_CLL = NeighbouringCellLinkedLists(particle_pos, radius=1.2, box_size=box_size)

        #for pid in range(nr_particles):
            #print(sorted(NL_PL.get_particles_within_radius(pid)))
            #print(sorted(NL_CLL.get_particles_within_radius(pid)))

        print(set(NL_CLL.get_particles_within_radius(0)) == (set(NL_PL.get_particles_within_radius(0))))




#
#     print ("particle position at 4: ", particle_pos[4])
#     print ("indices: ", nlist.get_particles_within_radius(4))
#     print ("particles close to 4", particle_pos[nlist.get_particles_within_radius(4)])
#     for i in nlist.get_particles_within_radius(4):
#         periodic_distance = np.linalg.norm(0.5*box_size - (particle_pos[4] - particle_pos[i] + 0.5*box_size)%box_size)
#         print (i, particle_pos[i], particle_pos[4], periodic_distance)
#
