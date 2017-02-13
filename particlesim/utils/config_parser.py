import pkg_resources
import glob
import configparser
import re
import numpy as np
from particlesim.utils import distributions
from os.path import join, abspath, isdir
from os import makedirs
from numpy import genfromtxt
from particlesim.api import SystemConfiguration


class ProblemCreator(object):

    default_config_path = pkg_resources.resource_filename('particlesim', 'lib/default_config.py')


    def __init__(self, user_config):
        '''
        Parameters
        ----------
        user_config : String
            Path to the .cfg config file
        '''
        # First read in the Librarys with atomic constants
        self._lib = self._read_lib()

        # Read configurations
        self.config = configparser.ConfigParser()
        self.config.read(self.default_config_path)  # first default configuration
        self.config.read(user_config)               # then user config

        # Set modes for this simulation
        self.use_ewald = self.config['ewald_summation'].getboolean('use_ewald')
        self.use_lj = self.config['ewald_summation'].getboolean('use_lennard_jones')
        self.manual = 'manual' in self.config

        # Initialize empty lists that are later used to create the problem
        self.positions = []
        self.sigmas = []
        self.epsilons = []
        self.charges = []
        self. labels = []
        self.sigma_ewald = None

        # Check necessary parameters

        # Box-Lenght
        if self.config.has_option('general', 'box-size'):
            self.box_size = self.config['general'].getfloat('box-size')
        else:
            raise ValueError("No box-size given")

        # Epsilon_r
        if self.config.has_option('general', 'epsilon_r'):
            self.epsilon_r = self.config['general'].getfloat('epsilon_r')
        else:
            raise ValueError("No epsilon_r given")

        # Sigma for Ewald Summation if it is activated
        if self.use_ewald:
            self.sigma_ewald = self.config['ewald_summation']['sigma']

        # Optional Parameters
        # r_cutoff
        if self.config.has_option('general', 'r_cutoff'):
            self.r_cutoff = self.config['general'].getfloat('r_cutoff')

        # k_cutoff
        if self.config.has_option('general', 'k_cutoff'):
            self.k_cutoff = self.config['general'].getfloat('k_cutoff')

        # sigma_c
        if self.config.has_option('general', 'sigma_c'):
            self.sigma_c = self.config['general'].getfloat('sigma_c')


        # Create particles as intended in the config file
        if self.manual:
            path = self.config['manual']['csv_path']
            self.initial_configuration = genfromtxt(path, delimiter = ",", skip_header = 1, unpack = True,
                dtype=[('label', '|S11'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('charge', '<f8'),
                ('epsilon', '<f8'),('sigma', '<f8')])
            self.n = len(self.initial_configuration[0])
        else:
            # Look for defined particle classes
            sections = self.config.sections()
            r = re.compile("particle_class_*")  # Regular expresssion
            particle_classes = filter(r.match, sections)

            if particle_classes == []:
                raise ValueError('No particles defined')

            for particle_class in particle_classes:
                self.add_particles(particle_class)



    def _read_lib(self):
        '''
        Reads in every cfg file found in /lib to create a dictionary for the atoms
        Keep in mind, that matching fields get overwritten by following reads.

        Returns
        -------
        config : configparser.ConfigParser
            Dictionary with atom information
        '''

        lib_dir = pkg_resources.resource_filename('particlesim', 'lib/')+'particle_types.py'
        log_files = glob.glob(lib_dir)

        config = configparser.ConfigParser()
        for log_file in log_files:
            config.read(log_file)

        return config


    def generate_problem(self):
        if self.manual:
            # Using CSV FILE
            # Unfiddle the input csv file
            try:
                readout = self.initial_configuration
                self.labels = self.initial_configuration['label'].astype(str, copy=False)
                self.positions = readout[['x', 'y', 'z']].view((float, len(readout[['x', 'y', 'z']].dtype.names)))
                self.charges = readout[['charge']].view((float, len(readout[['charge']].dtype.names)))
                self.epsilons = readout[['epsilon']].view((float, len(readout[['epsilon']].dtype.names)))
                self.sigmas = readout[['sigma']].view((float, len(readout[['sigma']].dtype.names)))

            except IndexError:
                print("The csv-file could not be read: The format of your csv-file does not have the required form. \
                      Even when LJ is deactivated, you have express some values for epsilons or sigmas")
                exit(1)
            arguments = dict(xyz=self.positions, sigmas=self.sigmas, epsilons=self.epsilons, charges=self.charges,
                             box_size=self.box_size, epsilon_r=self.epsilon_r , labels=self.labels)
        else:
            # Use local lists
            arguments = dict(xyz=np.array(self.positions), sigmas=np.array(self.sigmas),
                            epsilons=np.array(self.epsilons), charges=np.array(self.charges),
                            box_size=self.box_size, epsilon_r=self.epsilon_r, labels=self.labels)

        if hasattr(self, 'sigma_c'):
            arguments['sigma_c'] = self.sigma_c

        if hasattr(self, 'r_cutoff'):
            arguments['r_cutoff'] = self.r_cutoff

        if hasattr(self, 'k_cutoff'):
            arguments['k_cutoff'] = self.k_cutoff

        self.system_conf = SystemConfiguration(**arguments)

        return self.system_conf

    def add_particles(self, particle_class):
        '''
        Adding a group of particles read from the config list

        Parameters
        ----------
        particle_class : Name of the section in the particle class
        '''
        # Get parameters from dict
        # epsilons, sigmas, number of particles, labels
        type = self.config[particle_class]['type']
        epsilon, sigma = self.get_atom_params(type)

        n = self.config[particle_class].getint('number')
        distr = self.config[particle_class]['distribution']
        charge = self.config[particle_class].getfloat('charge')
        label = self.config[particle_class]['label']

        # Create positions as whished
        if(distr == 'uniform'):
            xyz = distributions.create_uniform_distribution(n, box_size = self.box_size)

        self.charges.extend([charge] * n)
        self.labels.extend([label] * n)
        self.epsilons.extend([epsilon] * n)
        self.sigmas.extend([sigma] * n)
        self.positions.extend(xyz.tolist())

    def get_atom_params(self, type):
        epsilon = self._lib[type].getfloat('epsilon')
        sigma = self._lib[type].getfloat('sigma')

        return epsilon, sigma

    def export_config(self):
        '''
        Exports the used config files. There might be differences with the
        user specified config file because of default parameters.

        The config file will be created in the log folder

        '''
        log_dir = pkg_resources.resource_filename('particlesim', 'logs/')

        # Check if logs exists
        if not isdir(log_dir):
            makedirs(log_dir)

        # Write to .cfg file
        config_path = abspath( join(log_dir, 'config.cfg'))

        with open(config_path, 'w') as configfile:
            self.config.write(configfile)

    def export_csv(self, positions, path):
        '''
        Generates an output csv file for the current configuration.

        Returns
        -------
        '''
        n = len(self.labels)
        labels = np.array(self.labels).reshape(n, 1)

        with open(path, 'wb') as file:
            file.write(str.encode("label, x, y, z, Charge, LJ-Epsilon, LJ-Sigma\n"))
            with_labels = np.concatenate((labels, positions, self.charges.reshape(n, 1),
                                          self.epsilons.reshape(n, 1), self.sigmas.reshape(n, 1)), axis=1)
            np.savetxt(file, with_labels, delimiter=",", fmt='%s')

    @staticmethod
    def convert_charmm_parrams(from_file, to_file):
        '''
        Converts atom parameters of Charmm Param File to our INI Config Format

        Parameters
        ----------
        from_file : path
            Path to the Charm param file

        Returns
        -------
        config : configparser.ConfigParser
            Config File generated from the Charmm Param Fila
        '''

        param_names = genfromtxt(from_file,
                                 comments='!',
                                 usecols= 0,
                                 dtype='|S5'
                                 ).astype(str, copy=False)
        epsilons_charmm, sigma_charmm = genfromtxt(from_file,
                                                   comments='!',
                                                   usecols= (2, 3),
                                                   unpack=True
                                                   )

        # Now lets convert the units from the units used in the charmm param file to atomic units
        # Epsilon from kcal/mol to Hartree Energy
        # ... TODO

        # Sigma from sigma/2 in Angstrom to 1/a_0
        # ... TODO

        config = configparser.ConfigParser()
        config['DEFAULT'] = {'label' : 'None'}

        for i in range(len(param_names)):
            config.add_section(param_names[i])
            config.set(param_names[i], 'epsilon', str(abs(epsilons_charmm[i])))
            config.set(param_names[i], 'sigma', str(sigma_charmm[i]))

        # Write to .cfg file
        with open(to_file, 'w') as configfile:
            config.write(configfile)
