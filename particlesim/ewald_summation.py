#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen, Simon Treu
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

import numpy as np
from particlesim.lib.converstion import prefactor
from .k_cython import calc_k_vectors

class EwaldSummation(object):

    def __init__(self, system_conf, sigma, k_cutoff):
        self.positions = system_conf.xyz # Original Positions
        self.charges = system_conf.charges
        self.volume = system_conf.volume
        self.sigma = sigma
        self.sigma_sq = sigma * sigma

        # Assig cutoff k and calculate vectors in k-space for longrange interaction energy
        self.k_cutoff = k_cutoff # multiple of 2*pi/L
        self.k_vectors = np.multiply(calc_k_vectors(k_cutoff), 2*np.pi/system_conf.box_size)


    def longrange_energy(self, positions):
        """
        Calculates the longrange potential and the self interaction potential
        of a given particle distribution using Ewald Summation.

        Atomic Units are used.

        Parameters
        ----------

        system_conf : np.array
            Array with particle positions and charges

        shape : int or float array
            Dimensions of the supercell

        sigma : float
            Standard deviation of gaussian distribution

        K : int
            Cutoff parameter in reciprocal space


        Returns
        -------

        float
            longrange and selfinteraction potential
        """
        self.positions  = positions

        N = len(self.positions)

        # Calculate longrange potential
        longrange_potential = 0.

        for k_i in self.k_vectors:
            structure_factor = 0.
            k_sq = np.linalg.norm(k_i) ** 2

            for a in range(N):
                structure_factor += self.charges[a] * np.e ** (1j * np.dot(k_i, self.positions[a]))
            structure_factor_squared = np.absolute(structure_factor) ** 2
            longrange_potential += structure_factor_squared * np.e ** (-self.sigma_sq * k_sq / 2) / k_sq

        longrange_potential *= 4 * np.pi / (2 * self.volume )

        # Calculate self-interaction potential
        self_interaction_potential = 0.

        for charge_i in self.charges:
            self_interaction_potential += charge_i**2
        self_interaction_potential *= 1 / (np.sqrt(2*np.pi) * self.sigma)

        # Calculate total potential
        longrange_and_self_potential = longrange_potential - self_interaction_potential

        return longrange_and_self_potential * prefactor

    def get_iterations(self):
        '''

        Returns
        -------
        it : int
            Number of steps for calculating the potential
        '''

        return len(self.k_vectors) * len(self.positions)