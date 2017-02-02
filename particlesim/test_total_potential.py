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
#

import numpy as np
from .helpers_for_tests import create_system_configuration
from .helpers_for_tests import create_positions

from .total_potential import TotalPotential

def test_longrange_potential_is_float():
    n = 100
    system_conf = create_system_configuration(n)
    total = TotalPotential(system_conf, k_cutoff=3)

    new_positions = create_positions(n)

    longrange = total.longrange_energy(new_positions)

#    print("Calculated longrange energy: %s eV" % longrange)
    assert isinstance(longrange, float)

    assert longrange != 0.
