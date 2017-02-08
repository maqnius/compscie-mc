import pkg_resources
import glob
import configparser
import re
from os.path import join, abspath, isdir
from os import makedirs
from numpy import genfromtxt
from particlesim.api import SystemConfiguration
from scipy import constants

class ProblemCreator(object):
    def __init__(self, config_file_path):
        '''
        Parameters
        ----------
        config_file_path : String
            Path to the .cfg config file
        '''
        # First read in the Librarys with atomic constants

        self._lib = self._read_lib()
        self.config = configparser.ConfigParser()
        self.config.read(config_file_path)

        # Check necessary parameters
        # CSV file with particle configuration

        # Box-Lenght
        if(self.config['general']['box-size']):
            self.box_size = self.config['general'].getfloat('box-size')
        else:
            raise ValueError("No box-size given")

        # Sigma for Ewald Summation if it is activated
        ewald = self.config['ewald_summation'].getboolean('use_ewald')
        if(ewald):
            try:
                self.sigma_ewald = self.config['ewald_summation']['sigma']
            except KeyError:
                ewald = False
                self.config['ewald_summation']['use_ewald'] = "no"
                print("No sigma for Ewald Summation found therefore Ewald summation got deactivated.")

        lennard_jones = self.config['lennard_jones'].getboolean('use_lennard_jones')

        # Create particles as intended in the config file
        if('manual' in self.config):
            path = self.config['manual']['csv_path']
            self.initial_configuration = genfromtxt(path, delimiter = ",", skip_header = 1, unpack = True)
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

        lib_dir = pkg_resources.resource_filename('particlesim', 'lib/')+'*.cfg'
        log_files = glob.glob(lib_dir)

        config = configparser.ConfigParser()
        for log_file in log_files:
            config.read(log_file)

        return config


    def generate_problem(self):
        # Unfiddle the input csv file
        try:
            positions = self.initial_configuration[:3].transpose()
            charges = self.initial_configuration[3].transpose()
            epsilons = self.initial_configuration[4].transpose()
            sigmas = self.initial_configuration[5].transpose()
        except IndexError:
            print("The csv-file could not be read: The format of your csv-file does not have the required form.")
            exit(1)

        self.system_conf = SystemConfiguration(positions, sigmas, epsilons, charges, self.box_size)

        return self.system_conf

    def add_particles(self, particle_class):
        '''
        Adding a group of particles read from the config list

        Parameters
        ----------
        particle_class : Name of the section in the particle class
        '''
        type = self.config[particle_class]['type']
        epsilon, sigma = self.get_atom_params(type)

        n = self.config[particle_class]['number']
        distr = self.config[particle_class]['distribution']
        label = self.config[particle_class]['label']


    def generate_positions(self, n, ):

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

    def export_csv(self, positions):
        '''
        Generates an output csv file for the current configuration.

        Returns
        -------
        '''
        pass

    def export_trajectory_to_vtk(self, trajectory):
        '''
        Exports the trajectory file to a vtk file for simulation

        Parameters
        ----------
        trajectory

        Returns
        -------

        '''
        pass

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
        config['DEFAULT'] = {'label' : ''}

        for i in range(len(param_names)):
            config.add_section(param_names[i])
            config.set(param_names[i], 'epsilon', str(epsilons_charmm[i]))
            config.set(param_names[i], 'sigma', str(sigma_charmm[i]))

        # Write to .cfg file
        with open(to_file, 'w') as configfile:
            config.write(configfile)
