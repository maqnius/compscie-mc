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
import numpy.testing as testing
import pyximport; pyximport.install()
import particlesim.k_cython as k_cython
import time

'''
def test_energy_is_positive():
    """
    Not sure if necessary
    """
    sigma, k_cutoff = calc_sigma_and_k_cutoff()
    shape = np.array([100, 100, 10])
    test_config = create_test_system(100, shape, 10)
    energy = longrange_energy(test_config, shape, sigma, k_cutoff)
    assert energy >= 0.
'''
def calc_sigma_and_k_cutoff():
    p = -np.log(1e-3)
    k_cutoff = 4
    sigma = np.sqrt(2*p/k_cutoff**2)

    return sigma, k_cutoff

def create_test_system(N, box_size, max_charge):
    """
    Creates a random particle configuration for testing other functions.

    Parameters
    ----------

    N : int
        Number of particles

    box_size : float
        Boxsize for particle positions

    max_charge : int
        Maximum for absolute values of particle charges


    Returns
    -------

    np.array
        Array with particle positions (3D) and charges
    """

    # Charges:
    charges_half = np.random.randint(-max_charge, max_charge, N // 2)
    charges = np.append(charges_half, -1*charges_half)

    # Positions:
    positions = np.zeros((N,3))
    for i in range(shape.shape[0]):
        positions[:, i] = np.random.randint(-box_size, box_size, N)

    test_config = [positions, charges]

    return test_config

# def test_k_vectors():
#     """
#     Speed test of calculation
#     """
#     K = 100
#
#     timestamp_start = time.time()
#     new = k_cython.calc_k_vectors(K)
#     timestamp_stop = time.time()
#     print("New Method took %s seconds." %(-timestamp_start + timestamp_stop,))
#
#     timestamp_start = time.time()
#     test = k_cython.calc_k_vectors_test(K)
#     timestamp_stop = time.time()
#     print("Test Method took %s seconds." %(-timestamp_start + timestamp_stop,))
#
#     timestamp_start = time.time()
#     old = calc_k_vectors_old(K)
#     timestamp_stop = time.time()
#     print("Old Method took %s seconds." %(-timestamp_start + timestamp_stop,))
#
#     testing.assert_array_equal(old, new)