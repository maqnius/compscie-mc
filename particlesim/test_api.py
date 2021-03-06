#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen,
#                       Simon Treu, Hanna Wulkow, Thomas Hadler
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
import pytest
from .api import *
from .helpers_for_tests import *

def test_system_configuration_potential_value():
    n_particles = 4
    system_configuration = create_system_configuration(n_particles, box_size=10)
    potential = system_configuration.potential(system_configuration.xyz, lennard_jones=True, coulomb=True)
    assert isinstance(potential, (float, int))


def test_add_particles_not_matching_input():
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charges = [1] * 10
    sigmas = [1.]*50 + [1.5] * 50
    with pytest.raises(TypeError):
        SystemConfiguration(xyz=particle_positions,charges=charges,
                                           sigmas=sigmas,
                                           epsilons=sigmas)
def test_system_configuration_sigmas_array_not_same_as_length():
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charges = [1.] * number_of_particles
    epsilons = [1.] * number_of_particles
    sigmas = [1.]*50
    with pytest.raises(TypeError):
        SystemConfiguration(xyz=particle_positions,charges=charges,
                                           sigmas=sigmas,
                                           epsilons=epsilons)

def test_system_configuration_sigmas_array_not_same_as_length():
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charges = [1.] * number_of_particles
    epsilons = [1.]*50
    sigmas = [1.] * number_of_particles
    with pytest.raises(TypeError):
        SystemConfiguration(xyz=particle_positions,charges=charges,
                                           sigmas=sigmas,
                                           epsilons=epsilons)

def test_add_particles_outside_box():
    particle_positions = np.array([[10,1,1]])
    charges = [1]
    sigmas = [1.]
    with pytest.raises(ValueError):
        SystemConfiguration(xyz=particle_positions,charges=charges,
                                           sigmas=sigmas,
                                           epsilons=sigmas,box_size=5)

def test_cutoff_too_large():
    number_of_particles = 10
    particle_positions = np.random.rand(number_of_particles,3)
    box_size = 5.0

    with pytest.raises(ValueError):
        SystemConfiguration(xyz = particle_positions, box_size = box_size)



def test_create_lennard_jones_sigmas():
    number_of_particles = 3
    particle_positions = np.random.rand(number_of_particles, 3)
    epsilons = [1] * 3
    sigmas = [2, 1.5, 1]
    charges = [1] * 3
    system_configuration = SystemConfiguration(xyz=particle_positions, charges=charges, epsilons=epsilons, sigmas=sigmas)
    np.testing.assert_array_equal(sigmas, system_configuration.lj_sigma_matrix.diagonal())

def test_only_input_positions_within_box_are_excepted():
    pass