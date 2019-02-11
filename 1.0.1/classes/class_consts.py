""" Copyright 2017 Tilo Zienert

    This file is part of cp-tools.

    cp-tools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cp-tools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cp-tools.  If not, see <http://www.gnu.org/licenses/>.

"""

class consts:
    def __doc__(self):
        """
            a central space for constants

            based on the 2014 CODATA recommended values
            see https://physics.nist.gov/cuu/index.html
        """
    def __init__(self):
        self._const_h=6.626070040e-34 # Planck constant in J*s
        self._const_k=1.38064852e-23 # Boltzmann constant in K/K
        self._const_N=6.022140857e23 # Avogadro constant in 1/mol
        self._const_R=8.3144598 # molar gas constant in J/K*mol
