from particlesim.utils import xyz
from os.path import dirname, abspath, join, isfile
import numpy as np

from particlesim.utils.config_parser import ProblemCreator


def export_trajectory():
    parent = dirname((dirname(dirname(__file__))))
    config_path = abspath(join(parent, 'example/example_config.cfg'))
    creator = ProblemCreator(config_path)

    system_conf = creator.generate_problem()

    fake_traj = [system_conf.xyz, system_conf.xyz] # Similar format as used in the simulation
    traj_path = abspath(join(parent, 'example/example_traj.xyz'))
    xyz.export_trajectory(system_conf.labels, fake_traj, traj_path)

    assert(isfile(traj_path))