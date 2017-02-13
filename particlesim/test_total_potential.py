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
from scipy.special import erfc
from particlesim.utils.conversion import prefactor
from .helpers_for_tests import create_system_configuration
from .helpers_for_tests import create_positions
from .api import SystemConfiguration
from .total_potential import TotalPotential

def test_longrange_potential_is_float():
    n = 100
    system_conf = create_system_configuration(n, box_size=10)
    total = TotalPotential(system_conf)

    new_positions = create_positions(n)

    longrange = total.longrange_energy(new_positions)
    # print("Calculated longrange energy: %s " % longrange)
    assert isinstance(longrange, float)

    assert longrange != 0.

def test_shortrange_potential_is_float():
    n = 100
    system_conf = create_system_configuration(n,box_size=10)
    total = TotalPotential(system_conf)

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
        system_conf.p_error = p_i
        total = TotalPotential(system_conf)

        # Cutoff should get not or not much smaller when raising the accuracy
        if not r_cutoff_prev == 0:
            assert(total.r_cutoff/r_cutoff_prev > 0.9)
        if not k_cutoff_prev == 0:
            assert(total.k_cutoff/(k_cutoff_prev) > 0.9)

        r_cutoff_prev = total.r_cutoff
        k_cutoff_prev = total.k_cutoff




def test_shortrange_coulomb_with_4_charges():
    """
    Compare the calculated shortrange potential to a simple charge distribution of 4 point charges.
    The reference potential is calculated manually.
    """
    boxsize = 10.
    xyz = np.array([[0., 0., 0.], [0., 0., 1.], [0., 1., 0.], [0., 1., 1.]])
    charges = np.array([1, -1, -1, 1])
    n = len(charges)
    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize)
    potential = TotalPotential(system_conf)
    sigma = potential.sigma_c
    shortrange_pot = 0.0
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            dist = np.linalg.norm(xyz[i]-xyz[j])
            shortrange_pot += 0.5*charges[i]*charges[j]/dist * erfc(dist/(np.sqrt(2)*sigma))
    shortrange_pot *= 1/(4*np.pi) * prefactor
    shortrange_pot_sim = potential.shortrange_energy(xyz, lennard_jones=False)
    np.testing.assert_almost_equal(actual=shortrange_pot_sim, desired=shortrange_pot, decimal=5)


