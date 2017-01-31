import numpy as np
import copy
from .total_potential import *


class SystemConfiguration(object):
    r"""
    A class to save the system configuration parameters

    Parameter:
        epsilon_r   :   optional, float, default=1.0
                        Permittivity in medium, default in vacuum
        box_size    :   optional, float
                        shape of box, boundaries
    arithmetic mean for sigma and geometric mean for epsilon
    arithmetic = (a+b)/2; geometric : sqrt(a*a)
    Lorentz Berthelot Rule
    lj_cutoff = 2.5 * sigma
    """
    def __init__(self, box_size=1.0, epsilon_r=1.0):
        self.box_size = box_size
        self.epsilon_r = epsilon_r
        self.xyz = np.ndarray(shape=(0,3),dtype=float)
        self.charges = np.asarray([],dtype=float)
        self.sigmas = np.asarray([],dtype=float)
        self.epsilons = np.asarray([],dtype=float)
        self._total_potential = TotalPotential(self)

    def add_particles(self, xyz, charges, sigmas, epsilons):
        r"""
        adds the particles to the SystemConfiguration

        :param xyz: array(n,3), float, positions of n-particles
        :param charges: array(n,1), float, charges of the particles
        :param sigmas: array(n,1), float, Lennard-Jones params
        :param epsilons: array(n,1), float, Lennard-Jones params
        :return: nil

        """
        # check if number of particles match in all params, raise error if not
        if not len(xyz) == len(charges):
            raise TypeError('charges must have the same length as particle numbers')
        if not len(xyz) == len(sigmas):
            raise TypeError('sigmas must have the same length as particle numbers')
        if not len(xyz) == len(epsilons):
            raise TypeError('epsilons must have the same length as particle numbers')

        # append new particles configuration to existing configuration
        self.xyz = np.concatenate((self.xyz, xyz), axis=0)
        self.charges = np.append(self.charges,charges)
        self.sigmas = np.append(self.sigmas, sigmas)
        self.epsilons = np.append(self.epsilons, epsilons)
        self.create_lj_mean_parameters()

    def add_particles_same_type(self, xyz, charge = 0., sigma = 1.0, epsilon = 1.0):
        r"""
        Add particles with same values for charge, sigma and epsilon to the system configuration
        :param xyz: np.ndarray(n,3)
        :param charge: float, Default = 0
        :param sigma: float, Default = 0
        :param epsilon: float, Default = 0
        :return:

        """

        # append new particles configuration to existing configuration
        number_of_particles = len(xyz)
        self.xyz = np.concatenate((self.xyz, xyz), axis=0)
        self.charges = np.append(self.charges, np.asarray([charge]*number_of_particles))
        self.sigmas = np.append(self.sigmas,np.asarray([sigma]*number_of_particles))
        self.epsilons = np.append(self.epsilons,np.asarray([epsilon]*number_of_particles))
        self.create_lj_mean_parameters()

    def potential(self,xyz_trial):
        # TODO only stub
        return self._total_potential.potential(xyz_trial)

    def number_of_particle_types(self):
        return len(self.xyz)

    def create_lj_mean_parameters(self):
        self.create_lennard_jones_epsilons()
        self.create_lennard_jones_sigmas()

    def create_lennard_jones_epsilons(self):
        self.lj_epsilon_matrix = np.sqrt(np.array([self.epsilons]).transpose()*np.array([self.epsilons]))

    def create_lennard_jones_sigmas(self):
        self.lj_sigma_matrix = (np.array([self.sigmas]).transpose() + np.array([self.sigmas]))/2

class Sampler(object):
    r"""A sampler class for hamiltonian objects."""
    def __init__(self, system_configuration):
        if len(system_configuration.xyz) == 0:
            raise ValueError("no particle in system configuration")
        self.system_configuration = system_configuration

    def _update(self, xyz, pot, step, beta):
        xyz_trial = (xyz + 2.0 * self.system_configuration.box_size * step
                     * (np.random.rand(*xyz.shape)- 0.5))%self.system_configuration.box_size
        pot_trial = self.system_configuration.potential(xyz_trial)
        if pot_trial <= pot or np.random.rand() < np.exp(beta * (pot - pot_trial)):
            return xyz_trial, pot_trial
        return xyz, pot

    def metropolis(self, iteration_number, step=0.1, beta=1.0):
        r"""
        Perform a Metropolis MC sampling procedure.

        Parameters
        ----------
        system_configuration : object
            Encapsulates the system's positions and params, box_size and
            a function to compute potential energies.
        iteration_number : int
            Number of Metropolis update steps.
        step : float, optional, default=0.1
            Maximal size of an update move in each coordinate.
        beta : float, optional, default=1.0
            Inverse temperature factor (1/kT).

        Returns
        -------
        numpy.ndarray of float
            Configuration trajectory.
        numpy.ndarray of float
            Total interaction and external potential trajectory.

        """

    #   check input data
        if not isinstance(iteration_number,int) or iteration_number <= 0:
            raise ValueError("To sample you need at least one iteration step...\n"
                             "iteration_numer has to be a positive integer")
        if not isinstance(step,(float,int)) or step <= 0:
            raise ValueError("stepsize has to be a postive number")
        if not isinstance(beta,(float,int)) or beta <= 0:
            raise ValueError("beta has to be a postive number")

    #   create copy of instance and work with copy, so initial configuration is unchanged
        xyz_traj = [self.system_configuration.xyz]
        pot_traj = [self.system_configuration.potential(self.system_configuration.xyz)]

    #   perform metropolis
        for i in range(iteration_number):
            xyz, pot = self._update(
                xyz_traj[-1]
                , pot_traj[-1],
                step=step, beta=beta)
            xyz_traj.append(xyz)
            pot_traj.append(pot)
        return np.asarray(xyz_traj, dtype=np.float64), np.asarray(pot_traj, dtype=np.float64)

    # TODO
    def metropolis_sa(self, hamiltonian, size, step=0.1, beta=1.0):
        r"""
        Perform a Metropolis-based simulated annealing procedure.

        Parameters
        ----------
        hamiltonian : object
            Encapsulates the system's degrees of freedom
            and interactions.
        size : int
            Number of Metropolis update steps.
        step : float, optional, default=0.1
            Maximal size of an update move in each coordinate.
        beta : float, optional, default=1.0
            Initial inverse temperature factor (1/kT).

        Returns
        -------
        numpy.ndarray of float
            Configuration trajectory.
        numpy.ndarray of float
            Total interaction and external potential trajectory.

        """
        beta_values = 1.0 / np.linspace(1.0E-15, 1.0 / beta, size)[::-1]
        xyz_traj = [np.asarray(hamiltonian.xyz, dtype=np.float64)]
        pot_traj = [hamiltonian.potential()]
        for i in range(size):
            xyz, pot = self._update(
                hamiltonian,
                xyz_traj[-1], pot_traj[-1],
                step=step, beta=beta_values[i])
            xyz_traj.append(xyz)
            pot_traj.append(pot)
        hamiltonian.xyz[:] = xyz_traj[-1]
        return np.asarray(xyz_traj, dtype=np.float64), np.asarray(pot_traj, dtype=np.float64)
