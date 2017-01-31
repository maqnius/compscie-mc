import numpy as np
import pytest
from .api import *

#
# def test_sampler_mmc_returns_trajectory():
#     # assert sampler.mmc returns an instance of trajectory
#     sampler = Sampler()
#     iteration_number = 1000
#     system_configuration = SystemConfiguration() # kann man mock für Klasse definieren
#     mmc_result = sampler.marcov_mc(iteration_number, system_configuration)
#     assert mmc_result.isinstance(Trajectory)


def test_system_configuration_potential_value():
    # assert systemconfiguration.potential returns value type float
    system_configuration = SystemConfiguration()
    potential = system_configuration.potential()
    assert isinstance(potential, (float, int))


def add_basic_random_particle_group_to_system_config(system_configuration):
    number_of_particles = 100
    particle_positions = np.random.rand(number_of_particles,3)
    charge = [1.] * number_of_particles
    sigma = [1.]*50 + [1.5] * 50
    system_configuration.add_particles(xyz=particle_positions, charges=charge, sigmas=sigma, epsilons=sigma)


def test_add_particles():
    system_configuration = SystemConfiguration()
    add_basic_random_particle_group_to_system_config(system_configuration)
    add_basic_random_particle_group_to_system_config(system_configuration)
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
    charge = [1] * 10
    sigma = [1.]*50 + [1.5] * 50
    with pytest.raises(TypeError):
        system_configuration.add_particles(xyz=particle_positions, charges=charge, sigmas=sigma, epsilons=sigma)

