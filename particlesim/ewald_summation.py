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
from particlesim.utils.conversion import prefactor
from .k_cython import calc_k_vectors


class EwaldSummation(object):
    """
    Parameters
    ----------
    system_conf : :obj:
        System configuration containing all parameters for the simulation.
    sigma : float
        Standard deviation of Gaussian charge distribution.
    k_cutoff : float
        Cutoff-radius in reciprocal space.

    """
    def __init__(self, system_conf, sigma, k_cutoff):
        self.system_conf = system_conf
        self.positions = system_conf.xyz # Original Positions
        self.charges = system_conf.charges
        self.volume = system_conf.volume
        self.sigma = sigma
        self.sigma_sq = sigma * sigma

        # Assig cutoff k and calculate vectors in k-space for longrange interaction energy
        self._k_cutoff = k_cutoff # multiple of 2*pi/L
        self.k_vectors = np.multiply(calc_k_vectors(self._k_cutoff), 2*np.pi/system_conf.box_size)

    @property
    def k_cutoff(self):
        return self._k_cutoff

    @k_cutoff.setter
    def k_cutoff(self, value):
        self._k_cutoff = value
        self.k_vectors = np.multiply(calc_k_vectors(self._k_cutoff), 2*np.pi/self.system_conf.box_size)

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

        # Calculate longrange potential vecotrized

        k_sq = np.einsum('ij,ij -> i', self.k_vectors, self.k_vectors)

        structure_factor = np.einsum('ik, k', np.exp(1j * np.einsum('ji, ki', self.k_vectors, self.positions)), self.charges )
        structure_factor_squared = structure_factor.real**2 + structure_factor.imag**2

        longrange_potential = np.dot(structure_factor_squared.T, np.exp(-self.sigma_sq * k_sq / 2) / k_sq) * \
                              prefactor / (2 * self.volume )

        # Calculate self-interaction potential

        self_interaction_potential = np.dot(self.charges.T, self.charges) * \
                                     1 / (np.sqrt(2*np.pi) * self.sigma) * 1/( 4 * np.pi) * prefactor


        # Calculate total potential
        longrange_and_self_potential = longrange_potential - self_interaction_potential

        return longrange_and_self_potential

    def get_iterations(self):
        '''
        Returns
        -------
        it : int
            Number of steps for calculating the potential

        '''

        return len(self.k_vectors) * len(self.positions)