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



def test_coulomb_random():
    """
    Do 10 repetitions of test-function for coulomb potential.
    """
    test_repetitions = 10
    for i in range(test_repetitions):
        coulomb_random()


def test_shortrange_with_different_neighbouring():
    """
    Do 10 repetitions of test-function for the comparison of the two neighbouring methods.
    """
    test_repetitions = 10
    for i in range(test_repetitions):
        shortrange_with_different_neighbouring()


def test_lennard_jones_rondom():
    """
    Do 10 repetitions of test-function for the lennard jones potential.
    """
    test_repetitions = 10
    for i in range(test_repetitions):
        lennard_jones_rondom()



def create_test_system():
    """
    Create a system configuration with a number of particles distributed in a 3x3x3 box inside of a 120x120x120
    box. The boxsize is so big, that particles from neighbouring boxes are outside the cutoff radius.
    Calculate the shortrange energies for the test system with a naive method.

    Returns
    -------
    system_conf : object
        System configuration for testing

    test_potential : floats
        coulomb and lennard jones energy
    """
    n = 4
    boxsize = 120
    particle_box = 3
    xyz = np.random.uniform(0, particle_box, (n, 3))
    ones = np.ones(2)
    charges = np.append(ones, -1 * ones)
    # Lennard-Jones parameters
    epsilon_Na = 0.0469
    epsilon_Cl = 0.15
    epsilon_NaCl = np.sqrt(epsilon_Cl * epsilon_Na)
    epsilons = np.append(epsilon_Na * ones, epsilon_Cl * ones)
    sigma_Na = 1.21496
    sigma_Cl = 2.02234
    sigma_NaCl = 0.5 * (sigma_Na + sigma_Cl)
    sigmas = np.append(sigma_Na * ones, sigma_Cl * ones)

    system_conf = SystemConfiguration(xyz, box_size=boxsize, charges=charges, sigmas=sigmas, epsilons=epsilons,
                                      r_cutoff=8.)
    total_potential = TotalPotential(system_conf)
    # Calculate Test_potential
    coulomb = 0.
    lj = 0.
    sigma_lj = 0.
    epsilon_lj = 0.
    for i in range(n):
        for j in range(i):
            if charges[i] == 1 and charges[j] == 1:
                epsilon_lj = epsilon_Na
                sigma_lj = sigma_Na
            if charges[i] == -1 and charges[j] == -1:
                epsilon_lj = epsilon_Cl
                sigma_lj = sigma_Cl
            if charges[i]*charges[j] == -1:
                epsilon_lj = epsilon_NaCl
                sigma_lj = sigma_NaCl
            dist = np.linalg.norm(xyz[i] - xyz[j])

            coulomb += charges[i] * charges[j] / dist * erfc(dist / (np.sqrt(2) * total_potential.sigma_c))
            lj += 4 * epsilon_lj * ((sigma_lj / dist) ** 12 - (sigma_lj / dist) ** 6)

    coulomb *= (1 / (4 * np.pi)) * prefactor
    test_potential = coulomb, lj

    return system_conf, test_potential


def shortrange_with_different_neighbouring():
    """
    Test if the different neighbouring methods for the shortrange potential deliver the same result.
    """
    system_conf, test_potential = create_test_system()
    total_potential = TotalPotential(system_conf)

    potential_neigh_true = total_potential.shortrange_energy(system_conf.xyz)

    system_conf.neighbouring = False
    total_potential = TotalPotential(system_conf)
    potential_neigh_false = total_potential.shortrange_energy(system_conf.xyz)

    np.testing.assert_almost_equal(actual=potential_neigh_true, desired=potential_neigh_false, decimal=5)


def coulomb_random():
    """
    Test the shortrange coulomb energy with a number of particles distributed in a 3x3x3 box inside of a 120x120x120
    box. The boxsize is so big, that particles from neighbouring boxes are outside the cutoff radius.
    """
    system_conf, test_potential = create_test_system()
    total_potential = TotalPotential(system_conf)

    sim_coulomb = total_potential.shortrange_energy(system_conf.xyz, lennard_jones=False)
    test_coulomb = test_potential[0]

    np.testing.assert_almost_equal(actual=sim_coulomb, desired=test_coulomb, decimal=5)


def lennard_jones_rondom():
    """
    Test the lennard jones energy with a number of particles distributed in a 3x3x3 box inside of a 120x120x120
    box. The boxsize is so big, that particles from neighbouring boxes are outside the cutoff radius.
    """
    system_conf, test_potential = create_test_system()
    total_potential = TotalPotential(system_conf)

    sim_lennard_jones = total_potential.shortrange_energy(system_conf.xyz, coulomb=False)
    test_lennard_jones = test_potential[1]

    np.testing.assert_allclose(actual=sim_lennard_jones, desired=test_lennard_jones, rtol=0.001)
