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

    ###

    DESCRIPTION: calculates the debye temperature from elastic constants
    USAGE: python debye-temperature-elast-const.py , answering questions
"""

import math as m
import sys
sys.path.append('classes') # path to class defininitions
# own
from class_debye import debye
debye=debye(1)
#
# collecting input values
print('Give values at T = 0 K\n')
M=float(input('Molar mass in g/mol: '))
q=int(input('Number of atoms per formular unit? : '))
d=float(input('Density in g/cm3: '))
nu=float(input('Poisson-number? : '))
K=float(input('Bulk modulus in GPa: '))
#
# calculate the Debye temperature
TD1=debye.TD_Debye1912(M,d,K,nu)*m.pow(q,1./3.)
TD2,vm=debye.TD_Anderson1963(M,d,K,nu,q)
# print result
print('\nT_D= {0} K [Debye1912]\nT_D= {1} K [Anderson1963]\n\nspecific energy E_s= {2} KJ/g\n'.format(TD1,TD2*1e5,K/d))
print('Have a nice day old fellow!')
