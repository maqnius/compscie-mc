import numpy as np


class TotalPotential(object):
    # dummy Klasse fuer totales Potential, später ersetzen...
    def __init__(self, system_configuration):
        self.system_configuration = system_configuration

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
        self.charges = []
        self.sigmas = []
        self.epsilons = []

    def add_particles(self, xyz, charges, sigmas, epsilons):
        r"""
        :param xyz: array(n,3) positions of the particles
        :param charges: array(n,1)
        :param sigmas: array(n,1)
        :param epsilons: array(n,1)
        :return: nil

        adds the particles to the systemConfiguration
        """
        if not len(xyz) == len(charges):
            raise TypeError('charges must have the same length as particle numbers')
        if not len(xyz) == len(sigmas):
            raise TypeError('sigmas must have the same length as particle numbers')
        if not len(xyz) == len(epsilons):
            raise TypeError('epsilons must have the same length as particle numbers')
        self.xyz = np.concatenate((self.xyz, xyz), axis=0)
        self.charges += charges
        self.sigmas += sigmas
        self.epsilons += epsilons

    def add_particles_same_type(self, xyz, charge = 0., sigma = 1.0, epsilon = 1.0):
        r"""
        Add particles with same values for charge, sigma and epsilon to the system configuration
        :param xyz: np.ndarray(n,3)
        :param charge: float, Default = 0
        :param sigma: float, Default = 0
        :param epsilon: float, Default = 0
        :return:

        """
        number_of_particles = len(xyz)
        self.xyz = np.concatenate((self.xyz, xyz), axis=0)
        self.charges += [charge]*number_of_particles
        self.sigmas += [sigma]*number_of_particles
        self.epsilons += [epsilon]*number_of_particles

    def potential(self):
        # TODO only stub
        total_potential = TotalPotential(self)
        return total_potential.potential()

    def number_of_particle_types(self):
        return len(self.xyz)


class Sampler(object):
    def markov_mc(self, iteration_number, system_configuration):
        pass
