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
import scipy.constants

def ewald_summation(system_conf, volume, sigma, K):
    """
    Calculates the longrange potential and the self interaction potential
    of a given particle distribution using Ewald Summation.

    Parameters
    ----------

    system_conf : np.array
        Array with particle positions and charges

    volume : int or float
        Volume of a supercell

    sigma : float
        Standard deviation of gaussian distribution

    K : int
        Cutoff parameter in reciprocal space


    Returns
    -------

    float
        longrange and selfinteraction potential
    """


    positions = system_conf[0]
    charges = system_conf[1]
    N = len(positions)
    sigma_sq = sigma**2
    epsilon_0 = scipy.constants.epsilon_0
    k_vectors = []

    # Create all k-vectors with absolute value <= K
    for a in range(0, K+1):
        for b in range(0, int(np.sqrt(K**2 - a**2)) + 1):
            for c in range(0, int(np.sqrt(K**2 - a**2 - b**2)) + 1):
                k_vectors.append([a,b,c])

    k_vectors.pop(0)

    # Calculate longrange potential
    longrange_potential = 0.

    for k_i in k_vectors:
        structure_factor = 0.
        k_sq = np.linalg.norm(k_i)**2

        for a in range(N):
            structure_factor += charges[a] * np.e**(1j*np.dot(k_i, positions[a]))
        structure_factor_squared = np.absolute(structure_factor)**2
        longrange_potential += structure_factor_squared * np.e**(-sigma_sq * k_sq / 2) / k_sq

    longrange_potential *= 1/(2*volume*epsilon_0)

    # Calculate self-interaction potential
    self_interaction_potential = 0.

    for charge_i in charges:
        self_interaction_potential += charge_i
    self_interaction_potential *= 1/(2*epsilon_0*sigma*(2*np.pi)**(3/2))

    # Calculate total potential
    longrange_and_self_potential = longrange_potential - self_interaction_potential

    return longrange_and_self_potential