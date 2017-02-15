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