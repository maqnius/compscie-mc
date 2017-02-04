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

def lj_potential(r, sigma=1.0, epsilon=1.0):
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
    q = (sigma / r)**6
    return 4.0 * (epsilon * (q * (q - 1.0)))

def shortrange(xyz, lj_sigma, ewald_sigma, epsilon, charges):
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
    sigmalist = lj_sigma
    epslist = epsilon

    nlist = NeighbouringCellLinkedLists(xyz, radius=1.2, box_side_length=2)

    [n,m] = xyz.shape

    lj_interaction = 0
    coulomb_interaction = 0

    for i in range(0, n):
        particle1 = i
        sigma1 = sigmalist[particle1]
        epsilon1 = epslist[particle1]

        neighbors = nlist.get_particles_within_radius(particle1)

        lj_interaction_tmp = 0
        coulomb_tmp = 0
        for j in range(0,len(neighbors)-1):
            particle2 = neighbors[j]
            sigma2 = sigmalist[particle2]
            epsilon2 = epslist[particle2]
            sigma = sigma1
            epsilon = epsilon1

            if (sigma2 != sigma1):
                sigma = (sigma1 + sigma2) / 2

            if (epsilon2 != epsilon1):
                epsilon = (epsilon1 + epsilon2) / 2

            # r = np.linalg.norm(xyz[particle1, :] - xyz[particle2, :]) Periodischer Abstand muss berechnet werden
            # Consider periodic boundary conditions when calculating the distances
            r = np.linalg.norm(0.5 * nlist.box_side_length - (xyz[particle1] - xyz[particle2] +
                                                              0.5 * nlist.box_side_length) % nlist.box_side_length)

            # Lennard Jones Potential
            lj_interaction_tmp += lj_potential(r, sigma=sigma, epsilon=epsilon)

            # Shortrange Coulomb Energy
            coulomb_tmp = charges[particle1] * charges[particle2] / r * erfc(r / (np.sqrt(2) * ewald_sigma))

        lj_interaction += lj_interaction_tmp
        coulomb_interaction += coulomb_tmp

    return lj_interaction + coulomb_interaction * constants.eV # Energy Unit is eV


def external_potential(xyz, box_length=None):
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
    if box_length is None:
        return 0.0
    if np.all(xyz >= 0.0) and np.all(xyz <= box_length):
        return 0.0
    return np.inf

def phi(xyz, lj_sigma, ewald_sigma, epsilon, charges, box_length=None):
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
    return shortrange(xyz, lj_sigma, ewald_sigma, epsilon, charges) + external_potential(xyz, box_length=box_length)

#testing the functions

xyz = np.random.rand(10, 3) * 2.0
sigma = [1,1,1,1,1,1,1,1,1,1]
charges = 10*[1]
ewald_sigma = 0.01
epsilon = [5,1,2,3,4,6,7,8,9,5,5]
box_length = 2

#print(phi(xyz, sigma, ewald_sigma, epsilon, charges, box_length))

