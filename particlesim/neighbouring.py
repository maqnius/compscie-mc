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

class Neighbouring(object):

    def __init__(self, particle_positions):
        self.particle_positions = particle_positions
        self._neighbourlist = self._create_neighbourlist(particle_positions) # only internal
    
    @property
    def particle_postions(self):
        return self._particle_positions
    
    @particle_postions.setter
    def particle_postions(self, value):
        # make some asserts
        self._particle_positions = value

    # private methods
    
    def _create_neighbourlist(self,particle_positions): # only internal
        return []

    # public methods

    def get_particles_within_radius(self,particle_id,radius):
        pass



class Neighbouring_Cell_Linked_Lists(Neighbouring):
    def __init__(self, particle_positions):
        super(Neighbouring_Cell_Linked_Lists,self).__init__(particle_positions)