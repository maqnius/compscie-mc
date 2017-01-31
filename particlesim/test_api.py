import numpy as np
import pytest
from .api import *
from .helpers_for_tests import *

def test_system_configuration_potential_value():
    # assert systemconfiguration.potential returns value type float
    n_particles = 4
    system_configuration = create_system_configuration(n_particles)
    potential = system_configuration.potential(system_configuration.xyz)
    assert isinstance(potential, (float, int))

def test_add_particles():
    system_configuration = create_system_configuration(100)
    assert 100 == system_configuration.number_of_particle_types()

def test_add_particles_not_matching_input():
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charges = [1] * 10
    sigmas = [1.]*50 + [1.5] * 50
    with pytest.raises(TypeError):
        SystemConfiguration(xyz=particle_positions,charges=charges,
                                           sigmas=sigmas,
                                           epsilons=sigmas)
def test_create_lennard_jones_epsilons():
    number_of_particles = 3
    particle_positions = np.random.rand(number_of_particles, 3)
    sigmas = [1]*3
    epsilons = [2,3,4]
    charges =  [1]*3
    system_configuration = SystemConfiguration(xyz=particle_positions,charges=charges,sigmas=sigmas,epsilons=epsilons)
    np.testing.assert_array_equal(epsilons, system_configuration.lj_epsilon_matrix.diagonal())

def test_create_lennard_jones_sigmas():
    number_of_particles = 3
    particle_positions = np.random.rand(number_of_particles, 3)
    epsilons = [1] * 3
    sigmas = [2, 3, 4]
    charges = [1] * 3
    system_configuration = SystemConfiguration(xyz=particle_positions, charges=charges, epsilons=epsilons, sigmas=sigmas)
    np.testing.assert_array_equal(sigmas, system_configuration.lj_sigma_matrix.diagonal())

def test_only_input_positions_within_box_are_excepted():
    pass