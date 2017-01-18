import numpy as np
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
    particle_positions = np.random.rand(number_of_particles, 3)
    charge = np.array([1] * number_of_particles)
    system_configuration.add_particles_same_type(particle_positions, charge, sigma=1, epsilon=1, lj_cutoff=3)


def test_add_particles_same_type():
    system_configuration = SystemConfiguration()
    add_basic_random_particle_group_to_system_config(system_configuration)
    assert 0 < system_configuration.number_of_particle_types()

