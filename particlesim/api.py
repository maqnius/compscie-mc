import numpy as np
import scipy.constants as constants

class TotalPotential(object):
    # dummy Klasse fuer totales Potential, später ersetzen...

    # Needs to be estimated
    t_k = 1 # Runtime of one fourierspace interaction of Ewald Simmulation
    t_r = 1 # Runtime of one realspace interaction of Ewald Simmulation

    def __init__(self, system_configuration, p = 1e-5, k_cutoff = None):
        self.system_configuration = system_configuration
        self.p = p
        if(k_cutoff):
            self.k_cutoff = k_cutoff



    def __estimate_parameters(self):
        """

        Returns
        -------

        """
        pass

    def __r_cutoff_optimal(self):
        return np.sqrt(self.p/np.pi) * (t_k/t_r)**(1/6) L/(N)


    def potential(self):
        # TODO only stub
        return 0


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

    def potential(self):
        # TODO only stub
        total_potential = TotalPotential(self)
        return total_potential.potential()

    def number_of_particle_types(self):
        return len(self.xyz)


class Sampler(object):

    def _update(self, system_configuration, stepsize):
        xyz_trial = system_configuration.xyz + 2.0* stepsize * (np.random.rand(*system_configuration.xyz.shape) - 0.5)
        potential_trial = system_configuration.potential(xyz_trial)

    def markov_mc(self, iteration_number, stepsize, system_configuration):
        pass
