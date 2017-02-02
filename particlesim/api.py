import numpy as np
from .total_potential import *


class SystemConfiguration(object):
    r"""

    Parameters
    ----------
    xyz : ndarray(n,3), float
          position of n particles in x,y,z coordinates
    sigmas : ndarray(n) or float value
            sigma coefficient of lennard jones potential for each particle
            if not array but float value, assigned to all particles
            Default = 1.0 --> assigned to all particles
    epsilons : ndarray(n) or float value
                epsilon coefficient of lennard jones potential for each particle
                if not array but float value, assigned to all particles
                Default = 1.0 --> assigned to all particles
    charges : ndarray(n) or f<loat value
              charges coefficient of lennard jones potential for each particle
              if not array but float value, assigned to all particles
                Default = 0.0 --> assigned to all particles
    box_size : float,
                box_size for cubic simulation box, positive number
                Default = 1.0
    epsilon_r : float,
                relative permittivity constant of system
                Default = 1.0 --> for vacuum by definition

               arithmetic mean for sigma and geometric mean for epsilon
             arithmetic = (a+b)/2; geometric : sqrt(a*a)
             Lorentz Berthelot Rule
             lj_cutoff = 2.5 * sigma
    """

    def __init__(self, xyz, sigmas= 1.0, epsilons = 1.0, charges=0.0, box_size=5.0, epsilon_r=1.0):

        if not np.all((xyz>=0)*(xyz<box_size)):
            raise ValueError("xyz must be in range of zero to %d" %box_size)

        if isinstance(sigmas, (float,int)):
            sigmas = np.asarray([float(sigmas)] * len(xyz))
        elif not len(xyz) == len(sigmas):
            raise TypeError('sigmas must have the same length as particle numbers')

        if isinstance(epsilons, (float,int)):
            epsilons = np.asarray([float(epsilons)] * len(xyz))
        elif not len(xyz) == len(epsilons):
            raise TypeError('epsilons must have the same length as particle numbers')

        if isinstance(charges, (float,int)):
            charges = np.asarray([float(charges)] * len(xyz))
        elif not len(xyz) == len(charges):
                raise TypeError('charges must have the same length as particle numbers')

        self.box_size = box_size
        self.volume = box_size * box_size
        self.epsilon_r = epsilon_r
        self.xyz = xyz
        self.charges = charges
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.create_lj_mean_parameters()
        self._total_potential = TotalPotential(self)

        @property
        def xyz(self):
            return self._xyz

        @xyz.setter
        def xyz(self, value):
            xyz = np.asarray(value)
            if not issubclass(xyz.dtype.type, np.float):
                raise TypeError("values in xyz must be of type float")
            if xyz.ndim != 2 or xyz.shape[0] < 2 or xyz.shape[1] != 3:
                raise ValueError("xyz must be of shape=(n_particles, dim) with n_particles > 1 and dim = 3")
            self._xyz = xyz

        @property
        def box_size(self):
            return self._box_size

        @property
        def volume(self):
            return self._volume

        @box_size.setter
        def box_size(self, value):
            if not isinstance(value, (float, int)) or value <= 0.0:
                raise ValueError("box_size must be a positive number or None")
            self._box_size = float(value)

        @property
        def charges(self):
            return self._charges

        @charges.setter
        def charges(self, value):
            charges = np.asarray(value)
            if not issubclass(charges.dtype.type, np.float):
                raise TypeError("values of charges must be of type float")
            if charges.ndim != 2 or charges.shape[1] != 1:
                raise ValueError("charge must be of shape=(n_charges, dim) dim = 1")
            charges = np.asarray(value)

        @property
        def sigmas(self):
            return self._sigmas

        @sigmas.setter
        def sigmas(self, value):
            if not np.all(value >= 0):
                raise ValueError("sigmas must be positive float")

        @property
        def epsilons(self):
            return

        @epsilons.setter
        def epsilons(self, value):
            if not np.all(value >= 0):
                raise ValueError("epsilons must be positive float")






    # def add_particles_same_type(self, xyz, charge = 0., sigma = 1.0, epsilon = 1.0):
    #     r"""
    #     Add particles with same values for charge, sigma and epsilon to the system configuration
    #     :param xyz: np.ndarray(n,3)
    #     :param charge: float, Default = 0
    #     :param sigma: float, Default = 0
    #     :param epsilon: float, Default = 0
    #     :return:
    #
    #     """
    #
    #     # append new particles configuration to existing configuration
    #     number_of_particles = len(xyz)
    #     self.xyz = np.concatenate((self.xyz, xyz), axis=0)
    #     self.charges = np.append(self.charges, np.asarray([charge]*number_of_particles))
    #     self.sigmas = np.append(self.sigmas,np.asarray([sigma]*number_of_particles))
    #     self.epsilons = np.append(self.epsilons,np.asarray([epsilon]*number_of_particles))
    #     self.create_lj_mean_parameters()

    def potential(self,xyz_trial):
        # TODO only stub
        foo = self._total_potential.potential(xyz_trial)
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
