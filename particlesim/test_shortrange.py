#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen, Simon Treu, Thomas Hadler, Hanna Wulkow
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
from .shortrange import Shortrange
from .api import SystemConfiguration
from scipy.special import erfc


def test_shortrange_ewald():
    """
    Test Coulomb energy by comparing it to a simple test-system.
    """
    ewald_sigma = 1.0
    q = 2.0
    boxsize = 10.
    charges = np.array([q, -1 * q])

    xyz = np.array([[0., 0., 0.], [0., 0., 1.]])
    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize)
    shortrange = Shortrange(system_conf,sigma_c = ewald_sigma, r_cutoff=2.)
    theoretical_shortrange_energy = -4 * erfc(1/(np.sqrt(2)*ewald_sigma))
    shortrange_energy = shortrange.shortrange(xyz, lj=False)

    # Test the pbc
    xyz = np.array([[0., 0., 0.], [0., 0., 9.]])
    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize)
    shortrange = Shortrange(system_conf, sigma_c=ewald_sigma, r_cutoff=2.)
    shortrange_energy_periodic = shortrange.shortrange(xyz, lj=False)


    np.testing.assert_almost_equal(actual=shortrange_energy, desired=theoretical_shortrange_energy, decimal=5)
    np.testing.assert_almost_equal(actual=shortrange_energy_periodic, desired=theoretical_shortrange_energy, decimal=5)


