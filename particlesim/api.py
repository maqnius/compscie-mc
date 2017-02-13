from .total_potential import *

class SystemConfiguration(object):
    r"""
    Parameters
    ----------
    xyz : ndarray(n,3), float
        Position of n particles in x,y,z coordinates.
    sigmas : ndarray(n) or float value
        Sigma coefficient of lennard jones potential for each particle;
        if not array but float value, assigned to all particles.
        Default = 1.0 --> assigned to all particles
    epsilons : ndarray(n) or float value
        Epsilon coefficient of lennard jones potential for each particle;
        if not array but float value, assigned to all particles.
        Default = 1.0 --> assigned to all particles
    charges : ndarray(n) or float value
        Charges coefficient of lennard jones potential for each particle;
        if not array but float value, assigned to all particles.
        Default = 0.0 --> assigned to all particles
    box_size : float
        Boxsize for cubic simulation box; positive number.
        Default = 1.0
    epsilon_r : float,
        Relative permittivity constant of system.
        Default = 1.0 --> for vacuum by definition
    labels : array-like of string
        Additional information about the particles.
    p_error : int
        Max. error for the total ewald summation.
        Error = e^-p_error
        Default = 10
    r_cutoff : float
        Cutoff-radius for shortrange Ewald summation
        Default = None --> Optimal cutoff is calculated automatically
    k_cutoff : float
        Cutoff-radius for longrange Ewald summation in reciprocal space.
        Dafault = None --> Optimal cutoff is calculated automatically
    neighbouring : bool
        True: Use neighbouring list for calculation of shortrange energies.
        False: Calculate neighbouring with fast_distances function in cython.


    Notes
    -----
        arithmetic mean for sigma and geometric mean for epsilon
        arithmetic = (a+b)/2; geometric : sqrt(a*b)
        Lorentz Berthelot Rule
        lj_cutoff = 2.5 * sigma

    """

    def __init__(self, xyz, sigmas= 1.0, epsilons = 1.0, charges=0.0, box_size=12.0, epsilon_r=1.0, labels = [],
                    p_error=10, r_cutoff = None, k_cutoff = None, neighbouring = False):

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
        self._volume = box_size ** 3
        self.epsilon_r = epsilon_r
        self.xyz = xyz
        self.charges = charges
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.labels = labels
        self.r_cutoff = r_cutoff
        self.k_cutoff = k_cutoff
        self.create_lj_mean_parameters()
        self.create_lennard_jones_cutoff()
        self._neighbouring = neighbouring
        self.p_error = p_error
        self._total_potential = TotalPotential(self)

        if self.box_size <= 2 * max(self.lj_cutoff_matrix.max(),self._total_potential.r_cutoff):
            raise ValueError('Box_size to small. Box_size has to be twice the cutoff radius '
                             'of the Lennard Jones potential.\n'
                             'box_size = %f\n lj_max = %f, coulomb_cutoff(r_cutoff) = %f \n'
                             'set box_size to be larger than %f \n '
                             % (self.box_size, self.lj_cutoff_matrix.max(), self._total_potential.r_cutoff, 2 * max(self.lj_cutoff_matrix.max(),self._total_potential.r_cutoff))
                             )


    @property
    def p_error(self):
        return self._p_error

    @p_error.setter
    def p_error(self, value):
        if value <= 0:
            raise ValueError('p_error must be bigger than zero')
        self._p_error = value

    @property
    def neighbouring(self):
        return self._neighbouring

    @neighbouring.setter
    def neighbouring(self, value):
        if not isinstance(value,bool):
            raise TypeError
        self._total_potential.shortrange.neighbouring = value
        self._neighbouring = value
        pass
    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, value):
        xyz = np.asarray(value,dtype=float)
        if not (issubclass(xyz.dtype.type, np.float) or issubclass(xyz.dtype.type, np.integer)):
            raise TypeError("values in xyz must be of type float or int")
        if xyz.ndim != 2 or xyz.shape[0] < 2 or xyz.shape[1] != 3:
            raise ValueError("xyz must be of shape=(n_particles, dim) with n_particles > 1 and dim = 3")
        self._xyz = xyz

    @property
    def volume(self):
        return self._volume

    @property
    def box_size(self):
        return self._box_size

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
        if not (issubclass(charges.dtype.type, np.float) or issubclass(charges.dtype.type, np.integer)):
            raise TypeError("values of charges must be of type float or int")
        if charges.ndim != 1:
            raise ValueError("charges must be a 1 dim array")
        self._charges = np.asarray(value,dtype=np.float)

    @property
    def sigmas(self):
        return self._sigmas

    @sigmas.setter
    def sigmas(self, value):
        sigmas = np.asarray(value, dtype=np.float)
        if not np.all(sigmas >= 0):
            raise ValueError("sigmas must be positive float")
        self._sigmas = sigmas

    @property
    def epsilons(self):
        return self._epsilons

    @epsilons.setter
    def epsilons(self, value):
        epsilons = np.asarray(value,dtype=np.float)
        if not np.all(epsilons >= 0):
            raise ValueError("epsilons must be positive float")
        self._epsilons = epsilons

    def potential(self,xyz_trial, lennard_jones = True, coulomb = True):
        if not (type(lennard_jones) == bool and type(coulomb == bool)):
            raise TypeError('lennard_jones and coulomb must be booleans')

        xyz_trial = np.asarray(xyz_trial,dtype=float)
        if not (issubclass(xyz_trial.dtype.type, np.float) or issubclass(xyz_trial.dtype.type, np.integer)):
            raise TypeError("values in xyz must be of type float or int")
        if xyz_trial.ndim != 2 or xyz_trial.shape[0] < 2 or xyz_trial.shape[1] != 3:
            raise ValueError("xyz must be of shape=(n_particles, dim) with n_particles > 1 and dim = 3")

        return self._total_potential.potential(xyz_trial, lennard_jones, coulomb)

    def create_lj_mean_parameters(self):
        self.create_lennard_jones_epsilons()
        self.create_lennard_jones_sigmas()

    def create_lennard_jones_epsilons(self):
        self.lj_epsilon_matrix = np.sqrt(np.array([self.epsilons]).transpose()*np.array([self.epsilons]))

    def create_lennard_jones_sigmas(self):
        self.lj_sigma_matrix = (np.array([self.sigmas]).transpose() + np.array([self.sigmas]))/2

    def create_lennard_jones_cutoff(self):
        self.lj_cutoff_matrix = 2.5 * self.lj_sigma_matrix


class Sampler(object):
    r"""A sampler class for hamiltonian objects.

    Parameters
    ----------
    system_configuration : :obj:
        Instance of an SystemConfiguration Object that holds essential parameters
        previously set by the user.

    """
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

    def metropolis_sa(self, iteration_number, step=0.1, beta=1.0):
        r"""
        Perform a Metropolis-based simulated annealing procedure.

        Parameters
        ----------
        self : object,
            Encapsulates the system's positions and params, box_size and
            a function to compute potential energies.
        iteration_number : int
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
        if isinstance(beta, (float, int)):
            # beta determines maximum
            beta_values = 1.0 / np.linspace(1.0E-15, 1.0 / beta, iteration_number)[::-1]
        else:
            try:
                if len(beta) == iteration_number:
                    # Accept beta values
                    beta_values = beta
                elif len(beta) == 2:
                    # beta contains min and max value for beta
                    beta_values = 1.0 / np.linspace(1.0 / beta[1], 1.0 / beta[0], iteration_number)[::-1]
                else:
                    raise ValueError(
                        "beta must be float|int, touple with len 2 or touple with len equal to iteration number")
            except TypeError:
                print("beta must be float|int, touple with len 2 or touple with len equal to iteration number")
                exit(1)

        xyz_traj = [self.system_configuration.xyz]
        pot_traj = [self.system_configuration.potential(self.system_configuration.xyz)]

        for i in range(iteration_number):
            xyz, pot = self._update(xyz_traj[-1], pot_traj[-1],
                step=step, beta=beta_values[i])
            xyz_traj.append(xyz)
            pot_traj.append(pot)
        return np.asarray(xyz_traj, dtype=np.float64), np.asarray(pot_traj, dtype=np.float64)


