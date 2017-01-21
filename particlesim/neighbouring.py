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


class Neighbouring(object):

    def __init__(self, particle_positions, radius=float("inf")):
        self.particle_positions = particle_positions
        self.n                  = len(particle_positions)
        self.r                  = float(radius)
        self._neighbourlist     = self._create_neighbourlist() # only internal
    
    @property
    def particle_postions(self):
        return self._particle_positions
    
    @particle_postions.setter
    def particle_postions(self, value):
        # make some asserts
        self._particle_positions = value

    # private methods
    def _create_neighbourlist(self, particle_positions, radius): # only internal
        return []

    # public methods
    def get_particles_within_radius(self, particle_id):
        pass




class NeighbouringPrimitiveLists(Neighbouring):
    def __init__(self, particle_positions):
        super(NeighbouringPrimitiveLists, self).__init__(particle_positions, radius=float("inf"))

    # private methods
    def _create_neighbourlist(self):
        n, r, pos   = self.n, self.radius, self.particle_positions
        ret         = [[] for i in xrange(n)]
        for i in xrange(n):
            for j in xrange(n):
                if np.linalg.norm(pos[i]-pos[j])>=r or i==j: continue
                ret[i].append(j)
        return ret

    # public methods
    def get_particles_within_radius(self, particle_id):
        return self._neighbourlist[particle_id] #returns the indices of the points
        #return [self.particle_positions[i] for i in self._neighbourlist[particle_id]] # alternative: return the points themselves



class NeighbouringCellLinkedLists(Neighbouring):
    def __init__(self, particle_positions, box_side_length):
        self.box_side_length = float(box_side_length)
        super(NeighbouringCellLinkedLists, self).__init__(particle_positions)


    # private methods
    def _create_neighbourlist(self): # in O(self.n)
        n, r, pos           = self.n, self.r, self.particle_positions
        nr_cells            = int(self.box_side_length/r + 0.5)
        cell_linked_list    = [[[[] for i in xrange(nr_cells)] for j in xrange(nr_cells)] for k in xrange(nr_cells)]
        for i in xrange(n):
            x, y, z = (pos[i]/r).astype("int")
            cell_linked_list[x][y][z].append(pos[i])
        return cell_linked_list


    # public methods
    def get_particles_within_radius(self, particle_id):
        n, r, pos, cell_ll  = self.n, self.r, self.particle_positions, self._neighbourlist
        box_side_length     = self.box_side_length
        nr_cells            = len(cell_ll)
        ret = []

        p           = pos[particle_id]
        cell        = (p/r).astype("int")
        cell_dir    = np.rint((p%r)).astype("int"); cell_dir[cell_dir==0]=-1 # setting direction to nearest cell in xyz direction
        cells       = np.array([(cell+[x,y,z])%nr_cells for x in [cell_dir[0],0] for y in [cell_dir[1],0] for z in [cell_dir[2],0]])

        for cell_idx in cells:
            idx_x, idx_y, idx_z = cell_idx
            cell_list = cell_ll[idx_x, idx_y, idx_z]
            for neigh in cell_list:
                periodic_distance = np.linalg.norm(0.5*box_side_length - (p-neigh+0.5*box_side_length)%box_side_length)
                if periodic_distance>r or p==neigh: continue
                ret.append(neigh)

        return ret


if __name__=="__main__":
    pass












