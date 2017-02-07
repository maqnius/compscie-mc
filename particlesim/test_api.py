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
def test_create_systemconfig_with_single_sigma_epsilon_charge():
    number_of_particles = 5
    particle_positions = np.random.rand(number_of_particles, 3)
    assert number_of_particles == len(SystemConfiguration(particle_positions,sigmas=2.0,
                                                               epsilons=2.0, charges = 1).xyz)



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