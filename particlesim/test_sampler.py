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

# def test_all_sampled_particles_are_inside_box():
#     pass
#
# def test_that_sampler_returns_maximum_around_lj_induced_value():
#     sampler, system_configuration = create_sampler(4)
#     traj, pot = sampler.metropolis(system_configuration=system_configuration,iteration_number= 10000)
#     for i in range(1, 4):
#         for j in range(i):
#             hist, edges = np.histogram(np.linalg.norm(traj[1000:, j, :] - traj[1000:, i, :], axis=-1), bins=50)
#             np.testing.assert_almost_equal(actual=edges[hist.argmax()],desired=lj_r0, decimal=1,
#                                            err_msg="It might also be that there is only a local "
#                                                    "and no global max around lj_r0")
def test_cumulative_percentage_global_optimum():
    n_particle = 4
    sampler, system_configuration = create_sampler(n_particle)
    traj, pot = sampler.metropolis(system_configuration=system_configuration, iteration_number=100000)
    r_right = 0.8
    r_left = 0
    for i in range(1,n_particle):
        for j in range(i):
            distance = np.linalg.norm(traj[1000:, j, :] - traj[1000:, i, :], axis = -1)
            hist, edges = np.histogram(distance, bins=50)
            indices = np.where(np.logical_and(edges[1:]>=r_left, edges[1:]<r_right))
            cumulated_value = np.sum(hist[indices])/np.sum(hist)
            assert cumulated_value >= 0.5, "it is not certain that the global optimum is reached"