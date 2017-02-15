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

import numpy as np

def export_trajectory(labels, trajectory, path):
    """
    Exports a trajectory in xyz format

    Parameters
    ----------
    labels : list
        A list containing the labels of the atom

    trajectory : list[np.ndarray]
        A list containing the xyz np.ndarrays

    path : str
        The path to the export file
    """

    n = len(labels)
    labels = np.array(labels).reshape((n,1))

    with open(path, 'wb') as file:
        for xyz in trajectory:
            file.write(str.encode("%s\n" %n))
            file.write(str.encode("label x y z\n"))

            with_labels = np.concatenate((labels, xyz), axis=1)
            np.savetxt(file, with_labels, delimiter=" ", fmt = '%s')
