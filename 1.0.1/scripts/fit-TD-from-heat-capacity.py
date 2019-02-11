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

    DESCRIPTION: get the Debye temperature from fitting the Debye equation to T-cp data
    USAGE: python fit-TD-from-heat-capacity.py DATAFILE
    DATAFILE format: T  cp
"""

import sys
sys.path.append('classes') # path to class defininitions
import numpy as np
from scipy.optimize import leastsq
#
from class_debye import debye
# vars
debye=debye(1)
data_x=[]
data_y=[]
phases=[]

# read data from DATAFILE
fp=open(sys.argv[1],'r')

T_max=float(input('T_max for Debye-fitting? '))
for line in fp:
    if line[0]!='#':
        tmp=line.split()
        if float(tmp[0])<=T_max:
            data_x.append(float(tmp[0]))
            data_y.append(float(tmp[1]))
fp.close()
#
data_x=np.array(data_x)
data_y=np.array(data_y)
# in case of multiphase, give the mole fraction and the Debye temperature of the known phases
# the contribution of the unknown phase is then fitted
while(True):
    tmp=input('If more than one phase ..\nGive mole fraction of fix phase: ')
    if not tmp:
        break
    else:
        T_D=float(input('Give Debye temperature: '))
        phases.append([float(tmp),T_D])
if not len(phases):
    phases=0
# start fitting
p=[100] # T_D = 100 K as first guess
fit_res=leastsq(debye.debye_residuals,p,args=(phases,data_x,data_y))
T_D=fit_res[0][0]
#
print('\nDebye-temperature T_D= {0} K for cp until T= {1} K\n'.format(T_D,T_max))
