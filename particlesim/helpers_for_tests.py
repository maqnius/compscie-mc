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
from .api import *
import numpy as np


def create_system_configuration(number_of_particles, box_size, max_charge=1.0, r_cutoff=None, k_cutoff=None):
    r"""
    adds particles with random position, charge = 1, sigma = 1 and epsilon = 1 to the system_configuration

    Parameters
    ----------
    box_size : float
        Side-length of a cubic box.
    max_charge : float or int
        Max. positive or negative charge of test configuration.

    Returns
    -------
    :obj:
        System configuration for testing.
    """

    # Charges:
    charges_half = np.random.randint(-max_charge, max_charge, (number_of_particles + 1) // 2 )
    charges = np.append(charges_half, -1 * charges_half)
    charges = charges[:number_of_particles] #check that nr of charges is correct
    # Positions:
    positions = create_positions(number_of_particles, box_size)
    # Sigmas for Lennard-Jones potential
    sigma = [1.] * number_of_particles

    return SystemConfiguration(xyz=positions,
                                       charges=charges,
                                       sigmas=sigma,
                                       epsilons=sigma,box_size=box_size, r_cutoff=r_cutoff, k_cutoff=k_cutoff)


def create_positions(number_of_particles, box_size = 1.0):
    return np.random.rand(number_of_particles,3) * box_size


def create_sampler(number_of_particles, box_size):
    r"""
    Creates a basic system_configuration with random positions, charge = 1, sigma = 1 and epsilon = 1
    and a sampler object

    Parameters
    ----------
    number_of_particles : int
        Number of particles inside one box.
    box_size : float or int
        Side-length of a cubic box.

    Returns
    -------
    :obj:
        Sampler object for testing.
    :obj:
        System configuration for testing.
    """
    system_configuration = create_system_configuration(number_of_particles, box_size=box_size)
    sampler = Sampler(system_configuration)

    return sampler, system_configuration


def periodic_distance(pos1, pos2, box_size):
    r"""
    Calculate the distance between two particles with periodic boundary conditions.

    Parameters
    ----------
    pos1 : array-like of floats
        Position of particle 1.
    pos2 : array-like of floats
        Position of particle 2.
    box_size : int or float
        Side-length of the particle box.

    Returns
    -------
    float
        Distance between two particles.

    """
    return np.linalg.norm(0.5 * box_size - (pos1 - pos2 + 0.5 * box_size) % box_size)