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
from .neighbouring import Neighbouring
from .neighbouring import Neighbouring_Cell_Linked_Lists

def test_create_Neighbouring_instance():
    neighbouring = Neighbouring("some input")
    assert neighbouring.particle_positions == "some input"

def test_neighbour_structur_is_created():
    neighbouring = Neighbouring("some input")
    assert neighbouring._neighbourlist == []
    # usually we would not call private (with _name) attributes from outside

def test_create_Neighbouring_Cell_Linked_Lists_instance():
    n = Neighbouring_Cell_Linked_Lists("something")
    assert "something" == n.particle_positions
