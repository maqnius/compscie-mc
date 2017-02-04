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

import numpy as np
from particlesim.neighbouring import NeighbouringCellLinkedLists
from scipy.special import erfc
import scipy.constants as constants


class LennardJones(object):
    def __init__(self, system_conf, sigma_c, r_cutoff, infinite_wall=False):
        """

        Parameters
        ----------
        system_conf : SystemConfiguration
            Instance of an SystemConfiguration Object that holds essential parameters
            previously set by the user

        sigma_c : Float
            Spread of the Gaussian function used for the Ewald Summation part

        """
        self.infinite_wall = infinite_wall
        self.epsilon_r = system_conf.epsilon_r
        self.box_length = system_conf.box_size
        self.charges = system_conf.charges
        self.epsilons = system_conf.lj_epsilon_matrix
        self.sigmas = system_conf.lj_sigma_matrix
        self.sigma_c = sigma_c

        # Create instance of neighbouring list
        self.nlist = NeighbouringCellLinkedLists(system_conf.xyz, r_cutoff,
                                                 self.box_length)

    def lj_potential(self, r, sigma=1.0, epsilon=1.0):
        r"""
        Compute the Lennard-Jones potential.

        Parameters
        ----------
        r : float or array-like of float
            Euklidean particle-particle distance(s).
        sigma : float, optional, default=1.0
            Zero crossing distance.
        epsilon : float, optional, default=1.0
            Depth of the potential well.
            Epsilon needs to be in eV

        Returns
        -------
        float or array-like of float
            Lennard-Jones energy value(s).

        """
        q = (sigma / r) ** 6
        return 4.0 * (epsilon * (q * (q - 1.0)))

    def shortrange(self, positions):
        r"""
        Compute the interaction potential for a pair of particles as a sum of the Lennard Jones Potential
        and the short coulomb interaction part of the ewald summation

        Parameters
        ----------
        xyz : numpy.ndarray(shape=(n, d))
            d-dimensional coordinates of n particles.

        li_sigma : numpy.ndarray(shape=(n, 1))
            List of zero crossing distances for the Lennard-Jones contribution.

        ewald_sigma: float
            Spread of the gaussian function in the Ewald Summation that makes a cutoff possible

        epsilon : numpy.ndarray(shape=(n, 1))
            List of depths of the potential well for the Lennard-Jones contribution.

        charges: numpy.ndarray(shape=(n,1))
            List of charges


        Returns
        -------
        float
            Total interaction potential in eV

        """

        [n, m] = positions.shape

        lj_interaction = 0
        coulomb_interaction = 0

        self.nlist.update_neighbourlist(positions)

        for i in range(0, n):
            particle1 = i

            neighbors = self.nlist.get_particles_within_radius(particle1)

            lj_interaction_tmp = 0
            coulomb_tmp = 0
            for j in range(0, len(neighbors) - 1):
                particle2 = neighbors[j]

                sigma = self.sigmas[i, j]
                epsilon = self.epsilons[i, j]

                r = np.linalg.norm(
                    0.5 * self.nlist.box_side_length - (positions[particle1] - positions[particle2] +
                                                        0.5 * self.nlist.box_side_length) % self.nlist.box_side_length)

                # Lennard Jones Potential
                lj_interaction_tmp += self.lj_potential(r, sigma=sigma, epsilon=epsilon)

                # Shortrange Coulomb Energy
                coulomb_tmp = self.charges[particle1] * self.charges[particle2] / r * erfc(
                    r / (np.sqrt(2) * self.sigma_c))

            lj_interaction += lj_interaction_tmp
            coulomb_interaction += coulomb_tmp

        return lj_interaction + coulomb_interaction * constants.eV  # Energy Unit is eV

    def external_potential(self, positions):
        r"""
        Compute the external potential for a set of particles.

        Parameters
        ----------
        xyz : numpy.ndarray(shape=(n, d))
            d-dimensional coordinates of n particles.
        box_length : float, optional, default=None
            If not None, the area outside [0, box_length]^d
            is forbidden for each particle.

        Returns
        -------
        float
            Total external potential.

        """
        if self.box_length is None:
            return 0.0
        if np.all(positions >= 0.0) and np.all(positions <= self.box_length):
            return 0.0
        return np.inf

    def phi(self, positions):
        r"""
        Compute the interaction and external potential for a set of particles.

        Parameters
        ----------
        ewald_sigma
        xyz : numpy.ndarray(shape=(n, d))
            d-dimensional coordinates of n particles.
        sigma : numpy.ndarray(shape=(n, 1))
            List of zero crossing distances for the Lennard-Jones contribution.
        epsilon : numpy.ndarray(shape=(n, 1))
            List of depths of the potential well for the Lennard-Jones contribution.
        box_length : float, optional, default=None
            If not None, the area outside [0, box_length]^d
            is forbidden for each particle.

        Returns
        -------
        float
            Total interaction and external potential.

        """
        if (self.infinite_wall):
            # With infinite Wall
            return self.shortrange(positions) + self.external_potential(positions)
        else:
            return self.shortrange(positions)
