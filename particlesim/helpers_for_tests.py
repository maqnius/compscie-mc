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


def create_system_configuration(number_of_particles, box_size=1.0, max_charge=1.0):
    r"""
    adds particles with random position, charge = 1, sigma = 1 and epsilon = 1 to the system_configuration
    :param system_configuration: object, instance of SystemConfiguration
    :param number_of_particles: int, positive int
    :return:

    Parameters
    ----------
    box_size
    max_charge
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
                                       epsilons=sigma)


def create_positions(number_of_particles, box_size = 1.0):
    return np.random.rand(number_of_particles,3) * box_size

def create_sampler(number_of_particles):
    r"""
    creates a basic system_configuration with random positions, charge = 1, sigma = 1 and epsilon = 1
    and a sampler object
    :param number_of_particles: int, positive integer
    :return: sampler, system_configuration
    """
    system_configuration = create_system_configuration(number_of_particles)
    sampler = Sampler(system_configuration)

    return sampler, system_configuration