import numpy as np
from .api import *


def test_sampler_mmc_returns_trajectory():
    # assert sampler.mmc returns an instance of trajectory
    sampler = Sampler()
    iteration_number = 1000
    system_configuration = SystemConfiguration(box_shape, particles) # kann man mock für Klasse definieren
    mmc_result = sampler.marcov_mc(iteration_number, system_configuration)
    assert mmc_result.isinstance(Trajectory)