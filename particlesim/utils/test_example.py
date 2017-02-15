#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen,
#                       Simon Treu, Hanna Wulkow, Thomas Hadler
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#

from particlesim.utils.config_parser import ProblemCreator
from particlesim.api import SystemConfiguration
from os.path import dirname, join, abspath
import configparser
import pkg_resources


def example_manual_configuration():
    # Create a test config file
    test_config_path = create_test_config()

    # Get a ProblemCreator file from the config file
    creator = ProblemCreator(test_config_path)

    # Initialize a SystemConfiguration Object from the config file
    system_config = creator.generate_problem()

    # Export the actual configuration that is used for the
    # SystemConfiguration Object in case some parameters couldn't
    # be set as the user whished
    creator.export_config()

    # Export csv file
    parent = dirname((dirname(dirname(__file__))))
    csv_path = abspath(join(parent, 'example/example_out_positions.csv'))
    creator.export_csv(system_config.xyz, csv_path)

    assert(isinstance(system_config, SystemConfiguration))

def create_lib():
    charmm_path = pkg_resources.resource_filename('particlesim', 'lib/particle_types.txt')
    to_file = pkg_resources.resource_filename('particlesim', 'lib/particle_types.py')
    lib = ProblemCreator.convert_charmm_parrams(charmm_path, to_file)

def example_type_configuration():
    # Create a test config file
    test_config_path = create_test_config_with_type_declaration()

    # Get a ProblemCreator file from the config file
    creator = ProblemCreator(test_config_path)

    # Initialize a SystemConfiguration Object from the config file
    system_config = creator.generate_problem()

    # Export the actual configuration that is used for the
    # SystemConfiguration Object in case some parameters couldn't
    # be set as the user whished
    creator.export_config()
    assert(len(system_config.labels) == len(system_config.charges))
    assert(system_config.epsilon_r == creator.epsilon_r)
    assert(isinstance(system_config, SystemConfiguration))


def create_test_config():
    # Find path to example_init_positions.csv
    # in ../../example/
    parent = dirname( (dirname ( dirname(__file__) ) ) )
    csv_path = abspath( join (parent, 'example/example_init_positions.csv'))

    # Create example configuration fields
    config = configparser.ConfigParser()

    config['general'] = {'box-size': 20.0}
    config['manual'] = {'csv_path': csv_path}

    config['ewald_summation'] = {'use_ewald': 'yes', 'sigma': 1.0}
    config['lennard_jones'] = {'use_lennard_jones': 'yes'}


    # Write to .cfg file
    config_path = abspath( join (parent, 'example/example_config.cfg'))
    with open(config_path, 'w') as configfile:
        config.write(configfile)

    return config_path


def create_test_config_with_type_declaration():
    # Find path to example_init_positions.csv
    # in ../../example/
    parent = dirname((dirname(dirname(__file__))))
    csv_path = abspath(join(parent, 'example/example_init_positions.csv'))

    # Create example configuration fields
    config = configparser.ConfigParser()

    config['general'] = {'box-size': 20.0, 'sigma_ewald': 1.0, 'k_cutoff': 1.0, 'r_cutoff': 3.0}

    # Create two types of atoms and their distribution
    config['particle_class_1'] = {'type': 'N', 'label': 'Natrium',
                                  'number': 10, 'charge': 1, 'distribution': 'uniform'}

    config['particle_class_2'] = {'type': 'C', 'label': 'Chlor',
                                  'number': 20, 'charge': -1, 'distribution': 'uniform'}

    # Write to .cfg file
    config_path = abspath(join(parent, 'example/example_config_types.cfg'))
    with open(config_path, 'w') as configfile:
        config.write(configfile)

    return config_path