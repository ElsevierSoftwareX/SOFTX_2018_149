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

    DESCRIPTION: calculates cv for a given Debye temperature and writes the results to savefile
    USAGE: python debye-int.py TD
"""

from scipy.integrate import quad
from decimal import *
import math as m
import sys

##########
# config #
##########

T = 1 # start temperature in K for cv calculation
T_end = 2000 # maximum temperature in K for cv calculation
TD = 400 # standard Debye temperature if none is given
savefile = "savefile" # name of the savefile, data will be saved in '{savefile}-{TD}'

##########

getcontext().prec=10 # necessary for very small Debye temperatures

# for simplicity, class_debye is not used
def debye_integrand(x):
    return (Decimal('{0}'.format(m.exp(Decimal(str(x)))))*Decimal(str(x))**4)/((Decimal('{0}'.format(m.exp(Decimal(str(x)))))-1)**2)

def cv(TD,T):
    N=Decimal('6.0221353e23')
    k=Decimal('1.3806503e-23')
    return 9*N*k*((T/TD)**3)*Decimal('{0}'.format(quad(debye_integrand,0,TD/T)[0]))

# get Debye temperature from standard input
arg_list=[]
for arg in sys.argv:
    arg_list.append(arg)

if len(arg_list)>1:
    TD=Decimal(arg_list[1])
else:
    TD=Decimal(str(TD))

# open savefile and write cv to it
fp=open('savefile-{0}'.format(TD),'w')
fp.write('# T[K]\tcv[J/K*mol]\n')

while T != T_end+0.5:
    fp.write('{0}\t{1}\n'.format(T,cv(TD,Decimal(str(T)))))
    if T<0.5:
        T=0.5
    else:
        T+=0.5

fp.close()

print('All data written to savefile-{0}\nHave a nice day old fellow!'.format(TD))
