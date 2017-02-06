import configparser
import os
from numpy import genfromtxt
from particlesim.api import SystemConfiguration

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
        self.ewald = True
        self.lennard_jones = True

        # Check necessary parameters
        if(self.config['general']['csv_path']):
            path = self.config['general']['csv_path']
            self.initial_configuration = genfromtxt(path, delimiter = ",", skip_header = 2, unpack = True)
            self.n = len(self.initial_configuration[0])
        else:
            raise ValueError("No link to initial configuration given")

        if(self.config['general']['box-size']):
            self.box_size = self.config['general'].getfloat('box-size')
        else:
            raise ValueError("No box-size given")

        if(self.config['ewald_summation']['use_ewald']):
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
        positions = self.initial_configuration[:3].transpose()
        charges = self.initial_configuration[3].transpose()
        epsilons = self.initial_configuration[4].transpose()
        sigmas = self.initial_configuration[5].transpose()

        self.system_conf = SystemConfiguration(positions, sigmas, epsilons, charges, self.box_size)

        return self.system_conf

    def export_config(self):
        '''
        Exports the used config files. There might be differences with the
        user specified config file because of default parameters.

        The config file will be created in the log folder

        '''
        parent = (os.path.dirname(os.path.dirname(__file__)))
        log_dir = os.path.abspath(os.path.join(parent, 'logs/'))

        # Check if logs exists
        if not os.path.isdir(log_dir):
            os.makedirs(log_dir)

        # Write to .cfg file
        config_path = os.path.abspath(os.path.join(log_dir, 'config.cfg'))
        with open(config_path, 'w') as configfile:
            self.config.write(configfile)

    def export_output_csv(self):
        '''
        Generates an output csv file for the current configuration

        Returns
        -------

        '''
        pass
