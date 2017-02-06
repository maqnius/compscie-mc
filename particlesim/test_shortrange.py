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
from .ewald_summation import EwaldSummation
from .total_potential import TotalPotential
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



def test_lj_potential():
    """
    Test Lennard-Jones energy by comparing it to a simple test-system.
    """
    ewald_sigma = 1.
    epsilon_Na = 0.0469
    epsilon_Cl = 0.15
    epsilon_NaCl = np.sqrt(epsilon_Cl*epsilon_Na)
    epsilons = np.array([epsilon_Na, epsilon_Cl])

    sigma_Na = 2*1.21496
    sigma_Cl = 2*2.02234
    sigma_NaCl = 0.5*(sigma_Na + sigma_Cl)
    sigmas = np.array([sigma_Na, sigma_Cl])

    xyz = np.array([[0., 0., 0.], [0., 0., 1.]])
    charges = np.array([1., -1.])
    boxsize = 10.

    thoeretical_lj_energy = 8 * epsilon_NaCl * (sigma_NaCl**12 - sigma_NaCl**6)

    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize, sigmas=sigmas, epsilons=epsilons)
    shortrange = Shortrange(system_conf, sigma_c=ewald_sigma, r_cutoff=2.)
    shortrange_energy = shortrange.shortrange(xyz, coulomb=False)

    np.testing.assert_almost_equal(actual=shortrange_energy, desired=thoeretical_lj_energy, decimal=5)


'''
def test_longrange_potential():
    """
    Test longrange Energy by comparing it to a simple test-system.
    """
    xyz = np.array([[0., 0., 0.], [0., 0., 1.]])
    charges = np.array([1., -1.])
    boxsize = 10.
    sigma = 1.
    k_cutoff = 1.

    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize)
    ewald_summation = EwaldSummation(system_conf, sigma=sigma, k_cutoff=k_cutoff)

    theoretical_longrange_energy = 1/(boxsize*2*np.pi) * np.exp(-sigma**2*(2*np.pi/boxsize)**2 / 2) *\
                                   (np.abs(1-np.exp(1j*2*np.pi/boxsize))**2 + np.abs(1-np.exp(-1j*2*np.pi/boxsize))**2) \
                                    - 2/(np.sqrt(2*np.pi*sigma))

    longrange_energy = ewald_summation.longrange_energy(xyz)

    np.testing.assert_almost_equal(actual=longrange_energy, desired=theoretical_longrange_energy, decimal=5)
'''


def test_total_potential():
    """

    """
    ewald_sigma = 1.
    epsilon_Na = 0.0469
    epsilon_Cl = 0.15
    epsilon_NaCl = np.sqrt(epsilon_Cl * epsilon_Na)
    epsilons = np.array([epsilon_Na, epsilon_Cl])

    sigma_Na = 2 * 1.21496
    sigma_Cl = 2 * 2.02234
    sigma_NaCl = 0.5 * (sigma_Na + sigma_Cl)
    sigmas = np.array([sigma_Na, sigma_Cl])

    xyz = np.array([[0., 0., 0.], [0., 0., 1.]])
    charges = np.array([1., -1.])
    boxsize = 10.

    system_conf = SystemConfiguration(xyz=xyz, charges=charges, box_size=boxsize, sigmas=sigmas, epsilons=epsilons)
    total_pot = TotalPotential(system_conf, sigma_c = 1., k_cutoff = 10, r_cutoff = 4)

    theoretical_potential = - 2 + 8 * epsilon_NaCl * (sigma_NaCl**12 - sigma_NaCl**6)

    potential = total_pot.potential(xyz)

    np.testing.assert_allclose(actual=potential, desired=theoretical_potential, rtol=0.01)
