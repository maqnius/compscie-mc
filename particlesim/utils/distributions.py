import numpy as np


def create_uniform_distribution(number_of_particles, box_size, dim = 3):
    return np.random.rand(number_of_particles, dim) * box_size