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


class Shortrange(object):
    def __init__(self, system_conf, sigma_c, r_cutoff):
        """

        Parameters
        ----------
        system_conf : SystemConfiguration
            Instance of an SystemConfiguration Object that holds essential parameters
            previously set by the user

        sigma_c : Float
            Spread of the Gaussian function used for the Ewald Summation part

        r_cutoff : int or Float
            Lenght of the realspace cutoff for calculating neighbours

        """
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

        Returns
        -------
        float or array-like of float
            Lennard-Jones energy value(s).

        """
        q = (sigma / r) ** 6
        return 4.0 * (epsilon * (q * (q - 1.0)))

    def shortrange(self, positions, coulomb=True, lj=True):
        r"""
        Compute the interaction potential for a pair of particles as a sum of the Lennard Jones Potential
        and the short coulomb interaction part of the ewald summation

        Parameters
        ----------
        xyz : numpy.ndarray(shape=(n, d))
            d-dimensional coordinates of n particles.

        Returns
        -------
        float
            Total interaction potential in Hartree-Energy

        """

        [n, m] = positions.shape

        lj_interaction = 0
        coulomb_interaction = 0

        #self.nlist.update_cells(positions)
        self.nlist.particle_positions = positions
        self.nlist.create_neighbourlist()

        for particle1 in range(0, n):

            neighbors, neigh_dists = self.nlist.get_particles_within_radius(particle1)
            neighbors = np.array(neighbors)
            neigh_dists = np.array(neigh_dists)
            sigma = np.array(self.sigmas)[[particle1], [neighbors]]
            charges = np.array(self.charges)
            if coulomb:
                coulomb_interaction += charges[particle1] * np.sum(charges[neighbors]/neigh_dists * erfc(neigh_dists/(np.sqrt(2) * sigma)))

            for j in range(len(neighbors)):
                particle2 = neighbors[j]
                r = neigh_dists[j]
                sigma = self.sigmas[particle1, particle2]
                epsilon = self.epsilons[particle1, particle2]



                # Lennard Jones Potential
                if lj:
                    lj_interaction += self.lj_potential(r, sigma=sigma, epsilon=epsilon)

        return 0.5 *(lj_interaction + coulomb_interaction)

