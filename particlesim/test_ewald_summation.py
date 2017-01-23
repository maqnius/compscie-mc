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
import scipy.constants
from .ewald_summation import ewald_summation

def test_energy_is_float():
    test_config = create_test_system(200, 100, 10)
    energy = ewald_summation(test_config, 100, 1, 4)
    assert isinstance(energy, float)


def test_energy_is_positive():
    """
    Not sure if necessary
    """
    test_config = create_test_system(200, 100, 10)
    energy = ewald_summation(test_config, 100, 1, 4)
    assert energy >= 0.

def create_test_system(N, size, max_charge):
    """
    Creates a random particle configuration for testing other functions.

    Parameters
    ----------

    N : int
        Number of particles

    size : int
        Boxsize for particle positions

    max_charge : int
        Maximum for absolute values of particle charges


    Returns
    -------

    np.array
        Array with particle positions (3D) and charges
    """

    # Charges:
    charges_half = np.random.randint(-max_charge, max_charge, N / 2)
    charges = np.append(charges_half, -1*charges_half) * scipy.constants.e

    # Positions:
    positions = np.random.randint(-size, size, (N, 3))

    test_config = [positions, charges]

    return test_config