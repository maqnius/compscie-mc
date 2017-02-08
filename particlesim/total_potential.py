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

class TotalPotential(object):
    r"""
    This class is initialized when a system_configuration is created. All calculations that are independent
    of the particle positions should be done at initialization
    """

    def __init__(self, system_configuration, sigma_c = 1., k_cutoff = 3, r_cutoff = 3, p_error = 10):
        self.sigma_coulomb = sigma_c
        self.p_error = p_error
        self.sigmas_lj = system_configuration.sigmas
        self.system_configuration = system_configuration
        self.longrange = EwaldSummation(system_configuration, self.sigma_coulomb, k_cutoff)

        # Create instance for calculation of shortrange energy
        self.shortrange = Shortrange(system_configuration, self.sigma_coulomb, r_cutoff)

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


    def _estimate_parameters(self, self.p_error):
        k_cutoff = 4
        sigma = np.sqrt(2 * p_error / k_cutoff)
        r_cutoff =

        return sigma, k_cutoff
        pass



    def _r_cutoff_optimal(self):
        pass