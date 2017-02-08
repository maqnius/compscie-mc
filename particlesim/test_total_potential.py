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
    system_conf = create_system_configuration(n, box_size=10)
    total = TotalPotential(system_conf, k_cutoff=3)

    new_positions = create_positions(n)

    longrange = total.longrange_energy(new_positions)
    # print("Calculated longrange energy: %s " % longrange)
    assert isinstance(longrange, float)

    assert longrange != 0.

def test_shortrange_potential_is_float():
    n = 100
    system_conf = create_system_configuration(n,box_size=10)
    total = TotalPotential(system_conf, k_cutoff=3)

    new_positions = create_positions(n)

    shortrange = total.shortrange_energy(new_positions)
    # print("Calculated shortrange energy: %s " % shortrange)
    assert isinstance(shortrange, float)

    assert shortrange != 0.


def test_parameter_guess():
    n = 100

    system_conf = create_system_configuration(n, box_size=10)
    r_cutoff_prev = 0
    k_cutoff_prev = 0

    for p_i in range(1, 10):
        total = TotalPotential(system_conf, p_error=p_i)

        # Cutoff needs to get bigger when raising the accuracy
        assert(r_cutoff_prev < total.r_cutoff)
        assert(k_cutoff_prev < total.k_cutoff)

        r_cutoff_prev = total.r_cutoff
        k_cutoff_prev = total.k_cutoff
