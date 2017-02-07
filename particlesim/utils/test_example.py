from particlesim.utils.config_parser import ProblemCreator
from particlesim.api import SystemConfiguration
from os.path import dirname, join, abspath
import configparser

def test_example_manual_configuration():
    test_config_path = create_test_config()
    creator = ProblemCreator(test_config_path)
    system_config = creator.generate_problem()
    creator.export_config()

    assert(isinstance(system_config, SystemConfiguration))

def example_type_configuration():
    test_config_path = create_test_config_with_type_declaration()
    creator = ProblemCreator(test_config_path)
    system_config = creator.generate_problem()
    creator.export_config()

    assert(isinstance(system_config, SystemConfiguration))

def create_test_config():
    # Find path to example_init_positions.csv
    # in ../../example/
    parent = dirname( (dirname ( dirname(__file__) ) ) )
    csv_path = abspath( join (parent, 'example/example_init_positions.csv'))

    # Create example configuration fields
    config = configparser.ConfigParser()

    config['DEFAULT'] = {'use_ewald': 'yes', 'use_lennard_jones': 'yes'}

    config['general'] = {'box-size': 5.0}
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

    config['DEFAULT'] = {'use_ewald': 'yes', 'use_lennard_jones': 'yes'}

    config['general'] = {'box-size': 5.0}
    config['ewald_summation'] = {'use_ewald': 'yes', 'sigma': 1.0}
    config['lennard_jones'] = {'use_lennard_jones': 'yes'}


    # Create two types of atoms and their distribution
    config['particle_class_1'] = {'type': 'na', 'label': 'Natrium',
                                  'number': 10, 'distribution': 'normal'}

    config['particle_class_2'] = {'type': 'cl', 'label': 'Chlor',
                                  'number': 20, 'distribution': 'normal'}

    # Write to .cfg file
    config_path = abspath(join(parent, 'example/example_config_types.cfg'))
    with open(config_path, 'w') as configfile:
        config.write(configfile)

    return config_path