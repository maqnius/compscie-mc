#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen, Simon Treu
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
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from .ewald_summation import Ewald_Summation


def test_create_Ewald_Summation():
    system_conf = conf_example()
    ewald = Ewald_Summation(system_conf)
    assert Ewald_Summation.system_conf == conf_example()


def test_energy_is_float():
    system_conf = conf_example()
    ewald = Ewald_Summation(system_conf)
    assert isinstance(Ewald_Summation.long_rangepotential(), float)


def test_energy_is_positive():
    system_conf = conf_example()
    ewald = Ewald_Summation(system_conf)
    assert Ewald_Summation.long_rangepotential() >= 0.


def conf_example():
    return "some conf"