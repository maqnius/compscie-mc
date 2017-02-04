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

import numpy as np
import pytest
from .neighbouring import Neighbouring
from .neighbouring import NeighbouringPrimitiveLists
from .neighbouring import NeighbouringCellLinkedLists

def test_create_Neighbouring_instance():
    with pytest.raises(TypeError):
        neighbouring = Neighbouring([1.,2.3,5.2])


def test_neighbour_structur_is_created():
    array = np.random.rand(5,3)
    neighbouring = NeighbouringCellLinkedLists(array)
    assert neighbouring._neighbourlist != None
    # usually we would not call private (with _name) attributes from outside

def test_create_Neighbouring_Cell_Linked_Lists_instance():

    with pytest.raises(TypeError):
        n = NeighbouringCellLinkedLists("something")


def test_equivalence_CLL_PL():
    box_size = float(7.6)
    nr_particles = 25
    pid = 5 % nr_particles
    particle_pos = np.random.rand(nr_particles, 3)*box_size
    NL_PL = NeighbouringPrimitiveLists(particle_pos, radius=1.2, box_size=box_size)
    NL_CLL = NeighbouringCellLinkedLists(particle_pos, radius=1.2, box_size=box_size)
    assert set(NL_CLL.get_particles_within_radius(pid)) == set(NL_PL.get_particles_within_radius(pid))