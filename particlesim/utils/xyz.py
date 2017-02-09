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
    labels = np.array(labels).reshape((30,1))

    with open(path, 'wb') as file:
        for xyz in trajectory:
            file.write(str.encode("%s\n" %n))
            file.write(str.encode("label x y z\n"))

            with_labels = np.concatenate((labels, xyz), axis=1)
            np.savetxt(file, with_labels, delimiter=" ", fmt = '%s')

            file.write(str.encode("\n"))