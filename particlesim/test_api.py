import numpy as np
import pytest
from .api import *
from .helpers_for_tests import *

def test_system_configuration_potential_value():
    # assert systemconfiguration.potential returns value type float
    system_configuration = SystemConfiguration()
    potential = system_configuration.potential()
    assert isinstance(potential, (float, int))

def test_add_particles():
    system_configuration = SystemConfiguration()
    add_basic_random_particle_group_to_system_config(system_configuration, 100)
    add_basic_random_particle_group_to_system_config(system_configuration, 100)
    assert 200 == system_configuration.number_of_particle_types()

def test_add_particles_of_same_type():
    system_configuration = SystemConfiguration()
    number_of_particles = 100
    epsilon=5.0
    particle_positions = np.random.rand(number_of_particles,3)
    system_configuration.add_particles_same_type(particle_positions, epsilon=epsilon)
    assert number_of_particles == system_configuration.number_of_particle_types()
    assert np.all(system_configuration.epsilons == np.asarray([epsilon]*number_of_particles))

def test_add_particles_not_matching_input():
    system_configuration = SystemConfiguration()
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charges = [1] * 10
    sigmas = [1.]*50 + [1.5] * 50
    with pytest.raises(TypeError):
        system_configuration.add_particles(xyz=particle_positions,
                                           charges=charges,
                                           sigmas=sigmas,
                                           epsilons=sigmas)

def test_sampler_negative_iteration_number():
    # might fail after potential function is implemented
    number_of_particles = 3
    sampler, system_configuration = create_sampler(number_of_particles)
    with pytest.raises(ValueError):
        sampler.metropolis(system_configuration, iteration_number=-1)
    with pytest.raises(ValueError):
        sampler.metropolis(system_configuration,iteration_number=1.5)

def test_sampler_trajectory():
    # might fail after potential function is implemented
    number_of_particles = 5
    iteration_number = 3
    sampler, system_configuration = create_sampler(number_of_particles)
    traj, pot = sampler.metropolis(system_configuration, iteration_number)
    assert np.any(np.not_equal(traj[0],traj[1]))
    assert len(traj) == iteration_number + 1

def test_sampler_no_particles_in_system():
    # might fail after potential function is implemented
    number_of_particles = 0
    sampler, system_configuration = create_sampler(number_of_particles)
    with pytest.raises(ValueError):
        sampler.metropolis(system_configuration,iteration_number=2)

def test_create_lennard_jones_epsilons():
    system_configuration = SystemConfiguration()
    number_of_particles = 3
    particle_positions = np.random.rand(number_of_particles, 3)
    sigmas = [1]*3
    epsilons = [2,3,4]
    charges =  [1]*3
    system_configuration.add_particles(xyz=particle_positions,charges=charges,sigmas=sigmas,epsilons=epsilons)
    np.testing.assert_array_equal(epsilons, system_configuration.lj_epsilon_matrix.diagonal())

    system_configuration.add_particles(particle_positions,charges,sigmas,epsilons)
    np.testing.assert_array_equal(np.append(epsilons, epsilons), system_configuration.lj_epsilon_matrix.diagonal())

def test_create_lennard_jones_sigmas():
    system_configuration = SystemConfiguration()
    number_of_particles = 3
    particle_positions = np.random.rand(number_of_particles, 3)
    epsilons = [1] * 3
    sigmas = [2, 3, 4]
    charges = [1] * 3
    system_configuration.add_particles(xyz=particle_positions, charges=charges, epsilons=epsilons, sigmas=sigmas)
    np.testing.assert_array_equal(sigmas, system_configuration.lj_sigma_matrix.diagonal())

    system_configuration.add_particles(particle_positions, charges,  sigmas, epsilons)
    np.testing.assert_array_equal(np.append(sigmas, sigmas), system_configuration.lj_sigma_matrix.diagonal())

def test_only_input_positions_within_box_are_excepted():
    pass