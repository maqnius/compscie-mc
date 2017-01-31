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


def add_basic_random_particle_group_to_system_config(system_configuration, number_of_particles):
    r"""
    adds particles with random position, charge = 1, sigma = 1 and epsilon = 1 to the system_configuration
    :param system_configuration: object, instance of SystemConfiguration
    :param number_of_particles: int, positive int
    :return:
    """
    particle_positions = np.random.rand(number_of_particles,3)
    charge = [1.] * number_of_particles
    sigma = [1.] * number_of_particles
    system_configuration.add_particles(xyz=particle_positions,
                                       charges=charge,
                                       sigmas=sigma,
                                       epsilons=sigma)

def create_sampler(number_of_particles):
    r"""
    creates a basic system_configuration with random positions, charge = 1, sigma = 1 and epsilon = 1
    and a sampler object
    :param number_of_particles: int, positive integer
    :return: sampler, system_configuration
    """
    system_configuration = SystemConfiguration()
    sampler = Sampler()
    add_basic_random_particle_group_to_system_config(system_configuration,
                                                     number_of_particles)
    return sampler, system_configuration