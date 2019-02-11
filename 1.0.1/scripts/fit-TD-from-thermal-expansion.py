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

    DESCRIPTION: get the Debye temperature using alpha_V-cv relation [Garai2006] to fit T-alpha_V data
    USAGE: python fit-TD-from-thermal-expansion.py DATAFILE
    DATAFILE format: T  alpha_V[1e6/K]
"""

# imports
import sys
sys.path.append('classes') # path to class defininitions

import numpy as np
from scipy.optimize import leastsq
# own
from class_debye import debye
from class_math_helpers import math_helpers

arg_list=[]
for arg in sys.argv:
    arg_list.append(arg)

if len(arg_list)>1:
    #some vars
    debye=debye(1)
    T=[]
    alpha=[]
    T_max=float(input('Give maximum fit temperature in K: '))
    #
    # read data
    fp=open(arg_list[1],'r')
    for line in fp:
        if line[0]!='#':
            tmp=line.split()
            if float(tmp[0])<T_max:
                T.append(float(tmp[0]))
                alpha.append(float(tmp[1]))
            else:
                break
    fp.close()
    #
    # ok, now fit
    T_max=T[len(T)-1] # highest exp. temperature
    T_min=T[0]
    p=[150] # first guess of the Debye temperature
    #
    fit_res=leastsq(debye.Garai2006_residuals,p,args=(np.array(T),np.array(alpha),T_min,T_max))
    T_D=fit_res[0][0]
    print('Debye temperature {0} K.'.format(T_D))
    #
    # ok, and now write the calculated alpha_V values to a savefile
    mh=math_helpers()
    mh.data_x=T
    mh.data_y=alpha
    debye.set_property('alpha_TD',mh.get_y(T_D))
    debye.set_property('debye_temperature',T_D)
    #
    with open('savefile-fit-te-'+arg_list[1].split('/')[len(arg_list[1].split('/'))-1],'w') as fp:
        fp.write('# datafile: {0}\n# obtained Debye temperature: {1} K\n#\n# T[K]\talpha_V[1e6/K]\n'.format(arg_list[1],T_D))
        for i in T:
            fp.write('{0}\t{1}\n'.format(i,debye.get_alpha_v_Garai2006(i)))
    #
    print('All written to {0}\nHave a nice day old fellow!'.format('savefile-fit-te-'+arg_list[1].split('/')[len(arg_list[1].split('/'))-1]))
else:
    print('No datafile given!')
