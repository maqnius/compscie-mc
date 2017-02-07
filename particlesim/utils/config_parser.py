import configparser
from os.path import dirname, join, abspath, isdir
from numpy import genfromtxt
from particlesim.api import SystemConfiguration
import pkg_resources

class ProblemCreator(object):
    def __init__(self, config_file_path):
        '''
        Parameters
        ----------
        config_file_path : String
            Path to the .cfg config file
        '''
        self.config = configparser.ConfigParser()
        self.config.read(config_file_path)
        self.manual = True
        self.ewald = True
        self.lennard_jones = True

        # Check necessary parameters
        # CSV file with particle configuration
        if(self.config['manual']['csv_path']):
            path = self.config['manual']['csv_path']
            self.initial_configuration = genfromtxt(path, delimiter = ",", skip_header = 1, unpack = True)
            self.n = len(self.initial_configuration[0])
        else:
            self.manual = False

        # Check for particle information
        # ...

        # Bos-Lenght
        if(self.config['general']['box-size']):
            self.box_size = self.config['general'].getfloat('box-size')
        else:
            raise ValueError("No box-size given")

        # Sigma for Ewald Summation if it is activated
        self.ewald = self.config['ewald_summation'].getboolean('use_ewald')
        if(self.ewald):
            try:
                self.sigma_ewald = self.config['ewald_summation']['sigma']
            except KeyError:
                self.ewald= False
                self.config['ewald_summation']['use_ewald'] = "no"
                print("No sigma for Ewald Summation found therefore Ewald summation got deactivated.")

        self.lennard_jones = self.config['lennard_jones'].getboolean('use_lennard_jones')

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

    def export_config(self):
        '''
        Exports the used config files. There might be differences with the
        user specified config file because of default parameters.

        The config file will be created in the log folder

        '''
        log_dir = pkg_resources.resource_filename('particlesim', 'lib/')

        # Check if logs exists
        if not isdir(log_dir):
            os.makedirs(log_dir)

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

#    def import_charmm_parrams(self):
