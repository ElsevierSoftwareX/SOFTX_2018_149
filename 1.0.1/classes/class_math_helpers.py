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

# all the misc math stuff for different functions
#
# imports
import numpy as np
from math import exp
from scipy.integrate import quad
from scipy.optimize import leastsq
#
class math_helpers:
    def __doc__(self):
        """
            Definition of some mathematical helper functions
        """
    def __init__(self):
        self.data_x=[]
        self.data_y=[]

    ###########################
    # interpolating functions #
    ###########################

    def get_slope(self,x1,x2,y1,y2):
        """
            calculates the slope between the two points P1(x1,y1) and P2(x2,y2)
        """
        if x1==x2:
            print('x1=x2 = {0}'.format(x1))
            return 0
        if (y2-y1)/(x2-x1):
            return (y2-y1)/(x2-x1)
        else :
            return 0

    def get_y(self,x):
        """
            returns the linearly interpolated y-value for a given x using the data stored in self.data_x and self.data_y
        """
        y=0
        l=len(self.data_x)
        for i in np.arange(0,l):
            if x==self.data_x[i]:
                y=self.data_y[i]
                break
            elif x>self.data_x[i] and x<self.data_x[i+1]:
                y=self.data_y[i]+(x-self.data_x[i])*self.get_slope(self.data_x[i],self.data_x[i+1],self.data_y[i],self.data_y[i+1])
                break
            elif i==(l-2) and x>self.data_x[i+1]:
                print('x outside of range!')
        return y

    ########################
    # integration routines #
    ########################

    def entropy_integrand(self,T):
        return self.get_y(T)/T

    def integrand(self,x):
        """
            integrand function for integration of stored values (data_x, data_y)
        """
        return self.get_y(x)

    def volume_integral(self,T,V0):
        return V0*exp(quad(self.get_y,0,T,limit=100)[0])

    # calculates the absolute integral
    def calc_abs_integral(self,zeropoints):
        erg=0
        for i in range(0,len(zeropoints)-1):
            if i==0:
                if zeropoints[i]==self.data_x[0]:
                    erg+=abs(quad(self.get_y,zeropoints[i],zeropoints[i+1],limit=100)[0])
                else:
                    erg+=abs(quad(self.get_y,self.data_x[0],zeropoints[i],limit=100)[0])
            else:
                erg+=abs(quad(self.get_y,zeropoints[i],zeropoints[i+1],limit=100)[0])
        #
        if zeropoints[len(zeropoints)-1] < self.data_x[len(self.data_x)-1]:
            erg+=abs(quad(self.get_y,zeropoints[len(zeropoints)-1],self.data_x[len(self.data_x)-1],limit=100)[0])
        #
        return erg

    ####################
    # fitting routines #
    ####################

    def data_fit(self,variant,data_x,data_y):
        """
            general interface to fit a given dataset
            using linear or quadratic function
        """
        data_x=np.array(data_x)
        data_y=np.array(data_y)
        #
        # first, guess a linear slope
        first_m=self.get_slope(data_x[0],data_x[len(data_x)-1],data_y[0],data_y[len(data_y)-1])
        if variant=='linear':
            p=[data_y[0],first_m]
            fit_res=leastsq(self.fit_linear_residuals,p,args=(data_x,data_y))
        elif variant=='linear+quad':
            p=[data_y[0],first_m,0]
            fit_res=leastsq(self.fit_linear_quad_residuals,p,args=(data_x,data_y))

        return fit_res[0]

    def fit_linear_residuals(self,p,x,y):
        """
            residuals function for linearly fitting
        """
        tmp=[]
        for i in x:
            tmp.append(p[0]+p[1]*i)

        return y-np.array(tmp)

    def fit_linear_quad_residuals(self,p,x,y):
        """
            residuals function for quadraticaly fitting
        """
        tmp=[]
        for i in x:
            tmp.append(p[0]+p[1]*i+p[2]*i**2)
        return y-np.array(tmp)
