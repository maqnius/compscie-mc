import numpy as np
class TotalPotential(object):
    # dummy Klasse fuer totales Potential, später ersetzen...
    def __init__(self, system_configuration):
        pass
    pass


class SystemConfiguration(object):
    r"""
    A class to save the system configuration parameters

    Parameter:
        epsilon_r   :   optional, float, default=1.0
                        Permittivity in medium, default in vacuum
        box_size    :   optional, float
                        shape of box, boundaries

    """

    def __init__(self, box_size,epsilon_r=1.0):
        self.box_size = box_size
        self.epsilon_r = epsilon_r
        pass
    def add_particles_same_type(self, xyz, charge, sigma, epsilon, lj_cutoff):
        pass


    def potential(self):
        total_potential=TotalPotential(self)
        return 0
    pass


class Sampler(object):
    def marcov_mc(self, system_configuration):
        pass
    pass

