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

    DESCRIPTION: calculates m from a fit of B(T) = B_T0 + m*T for T >= T_B0
    
    USAGE: python fit-m-from-B-T.py datafile

    datafile: tab-separated T[K]-B[GPa] values

"""

# imports
import sys
sys.path.append('classes')
#
from class_debye import debye
from class_datafile import datafile
from class_math_helpers import math_helpers
#
debye=debye(1)
df=datafile()
mh=math_helpers()
# vars
data_x=[]
data_y=[]
#
# ask for input paramters for cp-calculation
TD=float(input('Give Debye temperature [K]: '))
alpha_TD=float(input('Give alpha_TD [1e6/K]: '))*1e-6
V0=float(input('Give V0 [cm^3/mol]: '))*1e-6
B0=float(input('Give B0 [GPa]: '))*1e9
#
# set start parameters
debye.set_property('debye_temperature-only',TD)
debye.set_property('alpha_TD',alpha_TD)
debye.set_property('V0',V0)
debye.set_property('B0',B0)
#
# calculate B_T0 and T_B0
B_T0,T_B0=debye.get_T_B0()
print('calculated B_T0 = {0}, T_B0 = {1}'.format(B_T0,T_B0))
#
# ok, read experimental data
df.set_property('datafile',sys.argv[1])
tmp=df.get_data_values('x:y','ask','')
#
# get all B values at T >= T_B0
for i in range(0,len(tmp[0])):
    if tmp[0][i] >= T_B0:
        data_x.append(tmp[0][i])
        data_y.append(tmp[1][i]*1e9)
#
# call fitting routine
if len(data_x)>1:
    fit_res=mh.data_fit('linear',data_x,data_y)
    print('\nFIT: B0 = {0} GPa, m = {1} GPa/K\n'.format(fit_res[0]/1e9,fit_res[1]/1e9))
elif len(data_x)==1: # we do not need to fit anything
    m=mh.get_slope(T_B0,data_x[0],B_T0,data_y[0])
    print('\nslope calc: m = {0} GPa/K\n'.format(m/1e9))
#
print('All done, have a nice day old fellow!')
