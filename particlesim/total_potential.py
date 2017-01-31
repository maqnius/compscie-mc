#   particlesim
#   Copyright (C) 2017 Mark Niehues, Stefaan Hessmann, Jaap Pedersen, Simon Treu, Thomas Hadler, Hanna Wulkow
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
class TotalPotential(object):
    r"""
    This class is initialized when a system_configuration is created. All calculations that are independent
    of the particle positions should be done at initialization
    """
    def __init__(self, system_configuration):
        self.system_configuration = system_configuration
        # Todo ewald_long_range = EwaldLongRange()
        # TODO shortrange = Shortrange
        # #initialize the neighbouring datastructur. When there is a position update,
        #  just update the existing instance.
    def potential(self,xyz_trial):
        return 0

