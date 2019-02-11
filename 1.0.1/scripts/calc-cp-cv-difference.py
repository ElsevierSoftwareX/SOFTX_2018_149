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

    DESCRIPTION: calculates the cp-cv difference from cp data and a given Debye temperature

    USAGE: python calc-cp-cv-difference.py DATAFILE

    DATAFILE: tab-separated T-cp data

"""

# imports
import sys
sys.path.append('classes')
#
# own
from class_debye import debye
from class_datafile import datafile
#
TD=float(input('Give Debye temperature: '))
savefile=input('Give savefile: ')

debye=debye(TD)
df=datafile()
#
df.set_property('datafile',sys.argv[1]) # set datafile
data=df.get_data_values('x:y','ask',0) # and read the datapoints
# ok start calculation
fp=open(savefile,'w')
fp.write('#\n# File: {0}\n# TD = {1} K\n#\n#T[K]\tcp-cv[J/K*mol]\n'.format(sys.argv[1],TD))
#
for i in range(0,len(data[0])):
    fp.write('{0}\t{1}\n'.format(data[0][i],data[1][i]-debye.get_debye_integral(data[0][i])))
#
fp.close()
print('All written to {0}, have a nice day old fellow!'.format(savefile))
