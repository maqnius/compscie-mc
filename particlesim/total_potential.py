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
from particlesim.ewald_summation import EwaldSummation
from particlesim.shortrange import Shortrange
import numpy as np
import time


class TotalPotential(object):
    r"""
    This class is initialized when a system_configuration is created. All calculations that are independent
    of the particle positions should be done at initialization
    """

    def __init__(self, system_configuration, sigma_c = None, k_cutoff = None, r_cutoff = None, p_error = 10):
        self.sigma_c = sigma_c
        self.p_error = p_error
        self.sigmas_lj = system_configuration.sigmas
        self.system_configuration = system_configuration
        self.k_cutoff = k_cutoff
        self.r_cutoff = r_cutoff

        # Sigma shouldn't be set manually usually but for the sake of testing we leave it
        # as an optional parameter
        if sigma_c is None:
            if p_error <= 0:
                raise ValueError('p_error must be bigger than zero')

            # Set parameters according to the wanted accuracy p_error
            if all(x is None for x in [k_cutoff, r_cutoff]):
                # No parameters set at all - > estimate by time duration a
                # calculation in real and reziprocal space takes

                sigma_c_tmp = 1.0
                r_cutoff_tmp = 1.0
                k_cutoff_tmp = 3.0

                # Create instance for long range coulomb energy
                self.longrange = EwaldSummation(system_configuration, sigma_c_tmp, k_cutoff_tmp)
                # Create instance for calculation of shortrange energy
                self.shortrange = Shortrange(system_configuration, sigma_c_tmp, r_cutoff_tmp)

                self._estimate_all_parameters()

                # Update parameters
                self.longrange.sigma = self.sigma_c
                self.longrange.k_cutoff = self.k_cutoff

                self.shortrange.sigma_c = self.sigma_c
                self.shortrange.recreate_neighbourlist(self.r_cutoff)

            else:
                # At least one parameter set - > we can calculate the optimal value for the
                # missing parameter

                self._estimate_parameters()

                # Create instance for long range coulomb energy
                self.longrange = EwaldSummation(system_configuration, self.sigma_c, self.k_cutoff)
                # Create instance for calculation of shortrange energy
                self.shortrange = Shortrange(system_configuration, self.sigma_c, self.r_cutoff)
        else:
            if k_cutoff is None or r_cutoff is None:
                raise ValueError('k_cutoff and r_cutoff need to be set')

            # Create instance for long range coulomb energy
            self.longrange = EwaldSummation(system_configuration, self.sigma_c, self.k_cutoff)
            # Create instance for calculation of shortrange energy
            self.shortrange = Shortrange(system_configuration, self.sigma_c, self.r_cutoff)



    def longrange_energy(self, positions):
        return self.longrange.longrange_energy(positions)

    def shortrange_energy(self, positions, lennard_jones = True, coulomb = True):
        return self.shortrange.shortrange(positions, lj=lennard_jones, coulomb=coulomb)

    def potential(self, xyz_trial, lennard_jones = True, coulomb = True):
        pot = 0.
        if(coulomb):
            pot += self.longrange_energy(xyz_trial)
        pot += self.shortrange_energy(xyz_trial, lennard_jones, coulomb)
        return pot


    def _estimate_parameters(self):
        '''
        Estimates one missing cutoff parameter
        Returns
        -------
        '''
        if self.k_cutoff is None:
            self.k_cutoff = 2 * self.p_error / self.r_cutoff
        if self.r_cutoff is None:
            self.r_cutoff = 2 * self.p_error / self.k_cutoff

        self.sigma_c = np.sqrt(2 * self.p_error) / self.k_cutoff

    def _estimate_all_parameters(self):
        """

        Returns
        -------

        """
        it_longrange = self.longrange.get_iterations()
        it_shortrange = self.shortrange.get_iterations()

        start = time.time()
        self.longrange.longrange_energy(self.system_configuration.xyz)
        time_long = time.time() - start
        time_long *= 1 / it_longrange

        start = time.time()
        self.shortrange.shortrange(self.system_configuration.xyz, lj=False)
        time_short = time.time() - start
        time_short *= 1 / it_shortrange

        self.r_cutoff = np.sqrt(self.p_error/np.pi) * (time_long/time_short) ** (1/6) * \
                            self.system_configuration.box_size / len(self.system_configuration.xyz) ** (1/6)

        self._estimate_parameters()