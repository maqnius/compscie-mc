import configparser
from numpy import genfromtxt
import csv

class ProblemCreator(object):
    def __init__(self, config_file_path):
        self.parser = configparser.ConfigParser()
        self.config = self.parser.read(config_file_path)
        self.ewald = True

        # Check necessary parameters
        if(self.config['init_configuration']['csv_path']):
            path = self.config['init_configuration']['csv_path']
            self.initial_configuration = genfromtxt(path, delimiter = "\t", skip_header = 2)
        else:
            raise ValueError("No link to initial configuration given")

        if(self.config['ewald_summation']['use_ewald']):
            self.ewald = self.config['ewald_summation'].getboolean('use_ewald')
            if(self.ewald):
                self.sigma_ewald = self.config['ewald_summation']['sigma']
