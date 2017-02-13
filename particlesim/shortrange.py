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
from particlesim.k_cython import fast_distances
from scipy.special import erfc
from particlesim.neighbouring import NeighbouringCellLinkedLists
from particlesim.utils.conversion import prefactor


class Shortrange(object):
    def __init__(self, system_conf, sigma_c, r_cutoff):
        """
        Parameters
        ----------
        system_conf : SystemConfiguration
            Instance of an SystemConfiguration Object that holds essential parameters
            previously set by the user.
        sigma_c : Float
            Spread of the Gaussian function used for the Ewald Summation part.
        r_cutoff : int or Float
            Length of the real-space cutoff for calculating neighbours.
        """
        self.system_conf = system_conf
        self.epsilon_r = system_conf.epsilon_r
        self.box_length = system_conf.box_size
        self.charges = system_conf.charges
        self.epsilons = system_conf.lj_epsilon_matrix
        self.sigmas = system_conf.lj_sigma_matrix
        self.sigma_c = sigma_c
        self.distances = np.zeros((system_conf.xyz.shape[0],system_conf.xyz.shape[0]))
        self.r_cutoff = r_cutoff
        self.neighbouring = system_conf.neighbouring


        # Create instance of neighbouring list
        if(self.neighbouring):
            self.nlist = NeighbouringCellLinkedLists(system_conf.xyz, r_cutoff,
                                                    self.box_length)

    def lj_potential(self, r, sigma=1.0, epsilon=1.0):
        r"""
        Compute the Lennard-Jones potential.

        Parameters
        ----------
        r : float or array-like of float
            Euklidean particle-particle distance(s).
        sigma : float or array-like of float, optional, default=1.0
            Zero crossing distance(s).
        epsilon : float or array-like of float, optional, default=1.0
            Depth of the potential well.

        Returns
        -------
        float or array-like of float
            Lennard-Jones energy value(s).
        """
        q = (sigma / r) ** 6
        return np.sum(4.0 * (epsilon * (q * (q - 1.0))))


    def shortrange(self, positions, coulomb=True, lj=True):
        r"""
        Compute the interaction potential for a pair of particles as a sum of the Lennard Jones Potential
        and the short coulomb interaction part of the ewald summation

        Parameters
        ----------
        positions : numpy.ndarray(shape=(n, d))
            d-dimensional coordinates of n particles.
        coulomb : bool
            If true calculate coulomb potential.
        lj : bool
            If true calculate lennard jones potential.
        neighbouring : bool
            If false calculate all distances.

        Returns
        -------
        float
            Total interaction potential in Hartree-Energy.

        """

        [n, m] = positions.shape


        lj_interaction = 0
        coulomb_interaction = 0

        #self.nlist.update_cells(positions)

        if self.neighbouring:
            self.nlist.particle_positions = positions
            self.nlist.create_neighbourlist()


            for particle1 in range(0, n):

                neighbors, neigh_dists = self.nlist.get_particles_within_radius(particle1)
                if len(neighbors) == 0:
                    continue
                neighbors = np.array(neighbors)

                neigh_dists = np.array(neigh_dists)
                sigma = self.sigmas[[particle1], [neighbors]]
                epsilon = self.epsilons[[particle1], [neighbors]]
                charges = self.charges
                if coulomb:
                    coulomb_interaction += charges[particle1] * np.sum(charges[neighbors]/neigh_dists * erfc(neigh_dists/(np.sqrt(2) * self.sigma_c)))

                if lj:
                    lj_interaction += self.lj_potential(neigh_dists, sigma=sigma, epsilon=epsilon)

            return 0.5 * (lj_interaction + coulomb_interaction * 1/(4*np.pi) * prefactor)

        else:
            fast_distances(positions, box_len=self.box_length, distances=self.distances)
            for particle1 in range(0,n):
                lj_interaction_tmp, coulomb_interaction_tmp = 0,0
                for particle2 in range(particle1+1, n):
                    if lj:
                        sigma = self.sigmas[particle1, particle2]
                        epsilon = self.epsilons[particle1, particle2]
                        cutoff = self.system_conf.lj_cutoff_matrix[particle1, particle2]
                        distance = self.distances[particle1,particle2]

                        if distance < cutoff:
                            lj_interaction_tmp += self.lj_potential(distance,sigma,epsilon)
                    if coulomb:
                        cutoff = self.r_cutoff
                        distance = self.distances[particle1,particle2]

                        if distance < cutoff:
                            coulomb_interaction_tmp += self.charges[particle1] * self.charges[particle2] / distance * erfc(
                                distance / (np.sqrt(2) * self.sigma_c))

                lj_interaction +=  lj_interaction_tmp
                coulomb_interaction += coulomb_interaction_tmp
            return lj_interaction + coulomb_interaction * 1/(4*np.pi) * prefactor


    def get_iterations(self):
        """
        Number of itarations to estimate cutoff parameters.
        Returns
        -------
        int
            Total number of iterations to calculate the short-range potential.
        """
        if self.neighbouring:
            it = 0
            for particle1 in range(0, len(self.charges)):
                neighbors, neigh_dists = self.nlist.get_particles_within_radius(particle1)
                it += len(neighbors)
            return it
        else:
            return len(self.charges)*(len(self.charges)-1)/2


    def recreate_neighbourlist(self):
        '''
        Necessary for time measurement to estimate cutoff parameters.
        '''
        if self.neighbouring:
            self.nlist = NeighbouringCellLinkedLists(self.system_conf.xyz, self.r_cutoff,
                                                 self.box_length)