from particlesim.utils.config_parser import ProblemCreator
from particlesim.api import SystemConfiguration
from os.path import dirname, join, abspath
import configparser
import pkg_resources


def test_example_manual_configuration():
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

    assert(isinstance(system_config, SystemConfiguration))

def test_create_lib():
    charmm_path = pkg_resources.resource_filename('particlesim', 'lib/particle_types.txt')
    to_file = pkg_resources.resource_filename('particlesim', 'lib/particle_types.cfg')
    lib = ProblemCreator.convert_charmm_parrams(charmm_path, to_file)

def test_example_type_configuration():
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
    assert(isinstance(system_config, SystemConfiguration))


def create_test_config():
    # Find path to example_init_positions.csv
    # in ../../example/
    parent = dirname( (dirname ( dirname(__file__) ) ) )
    csv_path = abspath( join (parent, 'example/example_init_positions.csv'))

    # Create example configuration fields
    config = configparser.ConfigParser()

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

    config['general'] = {'box-size': 5.0}

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