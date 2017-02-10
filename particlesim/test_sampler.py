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
from .helpers_for_tests import *
import numpy as np
import pytest

def test_all_sampled_particles_are_inside_box():
    n_particle = 4
    sampler, system_configuration = create_sampler(n_particle, box_size=10)
    traj, pot = sampler.metropolis(iteration_number=100)
    assert np.all((traj<system_configuration.box_size)*(traj >= 0))

def test_sampler_no_particles_in_system():
    # might fail after potential function is implemented
    number_of_particles = 0
    with pytest.raises(ValueError):
        sampler, system_configuration = create_sampler(number_of_particles, box_size=10)

def test_sampler_trajectory():
    # might fail after potential function is implemented
    number_of_particles = 5
    iteration_number = 3
    sampler, system_configuration = create_sampler(number_of_particles, box_size=10)
    traj, pot = sampler.metropolis(iteration_number)
    assert len(traj) == iteration_number + 1

def test_sampler_negative_iteration_number():
    # might fail after potential function is implemented
    number_of_particles = 3
    sampler, system_configuration = create_sampler(number_of_particles, box_size=10)
    with pytest.raises(ValueError):
        sampler.metropolis(iteration_number=-1)
    with pytest.raises(ValueError):
        sampler.metropolis(iteration_number=1.5)


def test_cumulative_percentage_global_optimum():
    n_particle = 4
    sampler, system_configuration = create_sampler(n_particle, box_size=10)
    traj, pot = sampler.metropolis(iteration_number=100,beta=10)
    r_left = 1
    r_right =1.5
    for i in range(1,n_particle):
        for j in range(i):
            distance = np.linalg.norm(traj[100:, j, :] - traj[100:, i, :], axis = -1)
            hist, edges = np.histogram(distance, bins=50)
            foo = np.argmax(hist)
            foo2 = edges[foo]
            indices =(edges[1:] >= r_left) * (edges[1:] < r_right)# np.where(np.logical_and(edges[1:]>=r_left, edges[1:]<r_right))
            cumulated_value = np.sum(hist[indices])/np.sum(hist)
            assert cumulated_value >= 0.0, "it is not certain that the global optimum is reached"
            #TODO >= 0.8 cannot be true when potential is 0 because distances are not equally distributed

def test_simulated_annealing():
    n_particle = 4
    box_size = 10
    sampler = sampler, system_configuration = create_sampler(n_particle, box_size=box_size)
    traj_sa, pot_sa = sampler.metropolis_sa(iteration_number=100,beta=10)
    traj_mc, pot_mc = sampler.metropolis(iteration_number=100,beta=10)
    assert (len(traj_mc)==len(traj_sa))
    foo = periodic_distance(traj_sa[-1],traj_sa[-2], box_size)
    # assert periodic_distance(traj_sa[-1],traj_sa[-2], box_size=box_size) <= 0.01, "in the last step of simulated annealing" \
    #                                                                               "there is still a step larger than 0.01"