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
import scipy.constants as constants
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)
import particlesim.k_cython as k_cython

class EwaldSummation(object):

    def __init__(self, system_conf, sigma, k_cutoff):
        self.positions = system_conf.xyz # Original Positions
        self.charges = system_conf.charges
        self.volume = system_conf.box_size * system_conf.box_size
        self.sigma = sigma
        self.sigma_sq = sigma * sigma

        # Assig cutoff k and calculate vectors in k-space for longrange interaction energy
        self.k_cutoff = k_cutoff
        self.k_vectors = np.multiply(self.calc_k_vectors(k_cutoff), 2*np.pi/self.volume)


    def longrange_energy(self, positions):
        """
        Calculates the longrange potential and the self interaction potential
        of a given particle distribution using Ewald Summation.

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

        longrange_potential *= 1 / (2 * self.volume )

        # Calculate self-interaction potential
        self_interaction_potential = 0.

        for charge_i in self.charges:
            self_interaction_potential += charge_i**2
        self_interaction_potential *= 1 / (4 * np.pi * np.sqrt(2*np.pi) * self.sigma)

        #print("longrange_potential", longrange_potential)
        #print("self_interaction_potential", self_interaction_potential)

        # Calculate total potential
        longrange_and_self_potential = longrange_potential - self_interaction_potential

        return longrange_and_self_potential

    def calc_k_vectors_old(self, k_cutoff):
        """
        Old naive implementation.

        Calculates the k vectors of our lattice until a cutoff Value K.

        Parameters
        ----------

        K : int
            Cutoff value for the absolute value of the k-vectors


        Returns
        -------

        k_vectors : ndarray
            Array of k vectors that have an absolute value below cutoff K
        """
        k_vectors = []

        # Create all k-vectors with absolute value <= K

        for a in range(-k_cutoff, -k_cutoff + 1):
            b_limit = int(np.sqrt(-k_cutoff ** 2 - a ** 2))
            for b in range(-b_limit, b_limit + 1):
                c_limit = int(np.sqrt(-k_cutoff ** 2 - a ** 2 - b ** 2))
                for c in range(-c_limit, c_limit + 1):
                    k_vectors.append([a, b, c])

        # Remove k = [0, 0, 0]
        k_vectors.remove([0, 0, 0])

        return k_vectors

    def calc_k_vectors(self, K):
        return k_cython.calc_k_vectors(K)
