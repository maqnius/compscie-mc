from particlesim.utils.config_parser import ProblemCreator
import os
import configparser

def test_example_configuration():
    test_config_path = create_test_config()
    creator = ProblemCreator(test_config_path)
    system_config = creator.generate_problem()
    creator.export_config()

def create_test_config():
    # Find path to example_init_positions.csv
    # in ../../example/
    parent = os.path.dirname((os.path.dirname(os.path.dirname(__file__))))
    csv_path = os.path.abspath(os.path.join(parent, 'example/example_init_positions.csv'))

    # Create example configuration fields
    config = configparser.ConfigParser()

    config['DEFAULT'] = {'use_ewald': 'yes', 'use_lennard_jones': 'yes'}

    config['general'] = {'csv_path': csv_path, 'box-size': 5.0}
    config['ewald_summation'] = {'use_ewald': 'yes', 'sigma': 1.0}
    config['lennard_jones'] = {'use_lennard_jones': 'yes'}

    # Write to .cfg file
    config_path = os.path.abspath(os.path.join(parent, 'example/example_config.cfg'))
    with open(config_path, 'w') as configfile:
        config.write(configfile)

    return config_path