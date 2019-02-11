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

# Definitions of all modules for calculating Debye-stuff
#
# imports
import gmpy2
import math as m
import numpy as np
from scipy.integrate import quad
from scipy.optimize import leastsq
#
# own classes
from class_consts import consts
from class_math_helpers import math_helpers
from class_datafile import datafile
#
NONE=0
#
class debye:
    def __doc__(self):
        """ Implementation of the cp-predicting algorithm described in

                Zienert,Tilo and Fabrichnaya,Olga
                Prediction of heat capacity for crystalline substances
                Calphad XX (2019) XX-XX

                will be referred to as [Zienert2019]

                based on ideas and concepts from

                Debye, P.
                Zur Theorie der spezifischen Wärmen
                Annalen der Physik 344 (1912) 789-839

                here, referred to as [Debye1912]

                Anderson, O. L.
                A simplified method for calculation the Debye temperature from elastic constants
                Journal of Physics and Chemistry of Solids 24 (1963) 909-917

                here, reffered to as [Anderson1963]

                Raju, S.; Sivasubramanian, K. & Mohandas, E.
                An integrated thermodynamic approach towards correlating thermal and elastic properties: development of some simple scaling relations
                Solid State Communications 124 (2002) 151-156
                                        
                here, referred to as [Raju2002b]

                Garai,J.
                Correlation between thermal expansion and heat capacity
                Calphad 30 (2006) 354-356

                here, referred to as [Garai2006]

                OUTDATED:

                Raju, S.; Sivasubramanian, K. & Mohandas, E.
                On the thermodynamic interrelationship between enthalpy, volume thermal expansion and bulk modulus
                Scripta Materialia 44 (2001) 269-274

                referred here as [Raju2001] , lambda-concept

            covers:
                - cp prediction itself
                - fit routines to find
                    - the Debye temperature based  cp, volume or alpha_V data
                    - alpha_V at the Debye temperature and V_0 from experimental volume data
                    - B(T) from experimental data
                - estimate the Debye temperature from elastic constants, see [Debye1912] and [Anderson1963]
        """
    def __init__(self,debye_temperature):
        self.__debye_temperature=0
        # init vars with standard values
        self.__B0=0 # bulk modulus at T = 0 K in Pa
        self.__T_B0=0 # temperature in K at which B(T) = 0.9998*B_0 is fulfilled, in [Zienert2018]: T_trans^B
        self.__B_T0=0 # bulk modulus in Pa at T = T_B0 , in [Zienert2018]: B(T_trans^B)
        self.__B_m=0 # slope of B(T)
        self.__V0=0 # molar volume at T = 0 K in m^3
        self.__cv_TD=0 # value of c_V at the Debye temperature in J/K*mol
        self.__alpha_v_TD=0 # value of alpha_V at the Debye temperature
        self.__lambda=0
        self.__Raju2002b_epsilon=0 # const epsilon introduced by [Raju2002b]
        #
        self.__cp_T_start=1 # start temperature in K from which calculation will be started
        self.__cp_T_end=1500 # end temperature for calculation in K 
        self.__cp_delta_T=1 # stepsize for calculation in K
        self.__cp_variant='' # one of ['savefile','B(T)','V','fit-alpha_TD,V0','save-V(H)-TD']
        #
        """ vars Debye-fit
            equation: T_D(E_s)=a+b*E_s
            doc: to calculate the Debye temperature from values of the specific energy (see [Zienert2018]) """
        # values from [Zienert2019, equation (42)]
        self.__debye_fit_a = 38.912
        self.__debye_fit_b = 14.5631
        #
        """ vars alpha_v-fit
            equation: alpha_v(T_m)=a/(exp(b/T_m))+c+d*T_m
            doc: to calculate alpha_V^{T_D} from the correlation between alpha_V^{T_D} and the melting temperature (see [Zienert2018]) """
        # values from [Zienert2019, equation (44)]
        self.__alpha_fit_a = 233.72
        self.__alpha_fit_b = -191.496
        self.__alpha_fit_c = -234.183
        self.__alpha_fit_d = 0.00143
        #
        # filename to save results in
        self.__savefile='savefile'
        self.__fp_savefile=''
        # load necessary constants from 'class_consts'
        self.__c=consts()
        # final init: set the Debye temperature and c_V at the Debye temperature
        self.set_property('debye_temperature',debye_temperature)

    def set_property(self,variant,value):
        """ all global variables are set using the set_property function
            
            variant: the variable or name of set of variables
            value: value, type depends on variant
        """
        if variant=='debye_temperature': # set the Debye temperature and c_V at the Debye temperature
            self.__debye_temperature=value
            self.__cv_TD=self.get_debye_integral(value)
            return self.__cv_TD
        elif variant=='debye_temperature-only': # set only the Debye temperature
            self.__debye_temperature=value
        elif variant=='B0': # set bulk modulus at T = 0 K
            self.__B0=value
        elif variant=='V0': # set molar volume at T = 0 K
            self.__V0=value
        elif variant=='T_B0':
            self.__T_B0=value
        elif variant=='B_T0':
            self.__B_T0=value
        elif variant=='B_m': # set the linear slope of B(T)
            self.__B_m=value
        elif variant=='debye-fit': # set the variables to calculate the Debye temperature from the specific energy
            self.__debye_fit_a=value[0]
            self.__debye_fit_b=value[1]
        elif variant=='alpha_v-fit': # set the variables to calculate alpha_V at the Debye temperature from the melting temperature
            self.__alpha_fit_a=value[0]
            self.__alpha_fit_b=value[1]
            self.__alpha_fit_c=value[2]
            self.__alpha_fit_d=value[3]
        elif variant=='alpha_TD_fit': # set alpha_V at the Debye temperature using the T_ correlation, value=T_m
            self.__alpha_v_TD=self.get_alpha_v_fit(value)
            return self.__alpha_v_TD
        elif variant=='alpha_TD': # set alpha_V at the Debye temperature
            self.__alpha_v_TD=value
        elif variant=='lambda':
            """ OUTDATED OPTION
            
                set lambda, using the lambda-concept introduced in [Raju2001], value=B
            """
            V=self.get_volume(self.__V0,0,self.__debye_temperature,NONE,NONE,NONE) # get_volume at T_D
            print(self.__alpha_v_TD,self.__cv_TD,self.__debye_temperature,V,value)
            self.__lambda=gmpy2.mpfr(self.__alpha_v_TD)/(gmpy2.mpfr(self.__cv_TD)+(gmpy2.mpfr(self.__debye_temperature)*(gmpy2.mpfr(self.__alpha_v_TD)**2)*gmpy2.mpfr(V)*gmpy2.mpfr(value))) # see [Raju2001] for details
            return self.__lambda
        elif variant=='lambda_fix': # value=lambda
            self.__lambda=value
        elif variant=='epsilon': # set the epsilon value from [Raju2002b]
            self.__Raju2002b_epsilon=value
            print('Raju2002b_epsilon: {0}'.format(value))
        elif variant=='cp_calc': # set the start and end temperature in K for cp prediction
            self.__cp_T_start=value[0]
            self.__cp_T_end=value[1]
        elif variant=='savefile': # set the filename for saving the results
            self.__savefile=value
        elif variant=='fp_savefile': # set a filepointer to the savefile, return of open(savefile,'w')
            self.__fp_savefile=value
        elif variant=='cp_variant': # set the calculation variant, see __init__ description
            self.__cp_variant=value
        # generel cfg
        elif variant=='cp_T_end': # set only the end temperature for calculation
            self.__cp_T_end=value

    ########################
    # set cp-parameter cfg #
    ########################

    def set_cp_parameters_manual(self,variant):
        """ set all necessary variables for the chosen cp_variant calculation
            by answering questions

            variant = cp_variant , see __init__ description
        """
        self.set_property('cp_variant',variant)
        # ask questions
        B0=float(input('Give B0 in GPa: '))*1e9 # [B0]=Pa
        m=float(input('Give m in GPa/K: '))*1e9 # [m]=Pa/K
        V0=float(input('Give V0 in cm^3/mol: '))*1e-6 # [V0]=m^3
        self.set_property('V0',V0)
        self.set_property('B0',B0)
        self.set_property('B_m',m)
        #
        """ how to get alpha_V^{T_D} and the Debye temperature
            free: set manually
            fit: use the prediction correlation
        """
        fit=input('Using free cp? (y=free/n=fit): ')
        if 'n' in fit:
            # FIT
            if variant not in ['V','fit-alpha_TD,V0'] :
                T_m=float(input('Give T_m in K: ')) # melting temperature
                alpha_TD=self.set_property('alpha_TD_fit',T_m)
            E_s=float(input('Give E_s in kJ/g: '))
            T_D=self.get_debye_temperature(E_s)
            self.set_property('debye_temperature',T_D)
            print('T_D: {0} K, alpha_TD: {1} 1e6/K'.format(T_D,self.__alpha_v_TD*1e6))
        else:
            # FREE
            T_D=float(input('Give Debye temperature in K: '))
            self.set_property('debye_temperature',T_D)
            if variant not in ['V','fit-alpha_TD,V0','fit-alpha_TD,TD']:
                alpha_TD=float(input('Give alpha_v at the Debye temperature in 1e6/K: '))*1e-6
                self.set_property('alpha_TD',alpha_TD)
        #
        B_T0,T_B0=self.get_T_B0() # get B(T_trans^B) and T_trans^B according to [Zienert2019]
        self.set_property('T_B0',T_B0)
        self.set_property('B_T0',B_T0)
        #
        if variant=='savefile' or variant=='B(T)' or variant=='save-V(H)-TD':
            savefile=input('Give Savefile: ')
            self.set_property('savefile',savefile)
        #
        if variant not in ['V','fit-alpha_TD,V0','fit-alpha_TD,TD']:
            self.set_property('epsilon',self.get_Raju2002b_epsilon())
        #

    def set_cp_parameters_auto(self,V0,B0,m,T_D,alpha_TD,variant,savefile):
        """ set all variables for cp calculation automatically
            
            variant = cp_variant , __init__ description
        """
        self.set_property('V0',V0)
        self.set_property('B0',B0)
        self.set_property('B_m',m)
        self.set_property('debye_temperature',T_D)
        self.set_property('alpha_TD',alpha_TD)
        B_T0,T_B0=self.get_T_B0() # get B(T_trans^B) and T_trans^B according to [Zienert2019]
        self.set_property('T_B0',T_B0)
        self.set_property('B_T0',B_T0)
        self.set_property('epsilon',self.get_Raju2002b_epsilon()) # set epsilon [Raju2002b] using the V(T)-H(T) relation
        self.set_property('cp_variant',variant)
        self.set_property('savefile',savefile)

    #########################
    # routines bulk modulus #
    #########################

    def get_T_B0(self):
        """ calculates the temperature were B < 0.998 * B0 is fulfilled and the corresponding value of B
            which is T_trans^B and B(T_trans^B) in [Zienert2019]
            see equation (26) in [Zienert2019]
        """
        V_old=self.__V0
        T_old=0
        B=0
        i=0
        for i in np.arange(1,self.__debye_temperature,1): # T_max is the Debye temperature
            cv_tmp=self.get_debye_integral(i) # get c_V at T = i
            V_i=self.get_volume(V_old,T_old,i,NONE,NONE,cv_tmp) # calculate V at T=i with the alpha_V-c_V-correlation
            B=self.__B0*self.__V0/V_i
            #
            if self.__cp_variant=='B(T)': # if only the trend of B(T) is needed
                self.__fp_savefile.write('{0}\t{1}\n'.format(i,B*1e-9))
            #
            if ((100*B)/self.__B0) <= 99.98: # break when T = i = T_trans^B
                break
            V_old=V_i
            T_old=i
        #
        return B,i

    def get_B(self,T,V):
        """ calculates B(T) using equation (26) in [Zienert2019]
        """
        B=0
        if T<self.__T_B0:
            B=(self.__B0*self.__V0)/V # first assumption
        else:
            #m=-0.04014 # Al-10043
            #m=-0.02012 # Al-10071, 10087
            #m=-0.02813 # Au
            #m=-0.02227 # Ag
            #m=-0.0225 # Cu-10391, 10425
            #m=-0.01464 # Pb-10683
            #m=-0.015333 # Mo-10996
            #m=-1.818e-3 # K-11227
            #m=-1.37931e-3 # K-Schouten1974
            #self.set_property('B_m',-0.0204988*1e9) # a mean value
            #self.set_property('B_m',m*1e9)
            B=self.__B_T0+self.__B_m*(T-self.__T_B0)
        #
        return B

    def get_B_T(self):
        """ TEST ROUTINE
            
            saves B(T) after calculating m from a given B at one T
        """
        T=float(input('Give exp. B(T) \n\t for the Temperature in K: '))
        B=float(input('\tB in GPa: '))*1e9
        #
        mh=math_helpers()
        m=mh.get_slope(self.__T_B0,T,self.__B_T0,B)
        print('T_B0 = {0} [K], B_T0 = {1} [GPa], m = {2} [GPa/K]'.format(self.__T_B0,self.__B_T0*1e-9,m*1e-9))
        self.set_property('B_m',m)
        # save
        fp=open(self.__savefile,'w')
        fp.write('# T[K]\tB[GPa]\n')
        self.set_property('fp_savefile',fp)
        self.get_T_B0() # saves B(T) up to T = T_B0
        for i in np.linspace(self.__T_B0,self.__cp_T_end,num=(self.__cp_T_end-self.__T_B0)):
            fp.write('{0}\t{1}\n'.format(i,(self.__B_T0+m*(i-self.__T_B0))*1e-9))
        fp.close()
        return 'All done!' # IPC :-)
    
    ############################
    # routines volume, alpha_V #
    ############################

    def get_alpha_v_fit(self,T_m):
        """ calculates alpha_V at the Debye temperature from a given melting temperature
            using equation (42) in [Zienert2019]
        """
        return 1e-6*((self.__alpha_fit_a/(gmpy2.exp(self.__alpha_fit_b/T_m)))+self.__alpha_fit_c+self.__alpha_fit_d*T_m)

    def get_alpha_v_Garai2006(self,T):
        """ calculates alpha_V based on the relation described by [Garai2006]
            see equation (24) in [Zienert2019]
        """
        return (self.__alpha_v_TD/self.__cv_TD)*self.get_debye_integral(T)

    def get_alpha_v(self,T,B,V,cv_tmp):
        """ OUTDATED ROUTINE (partially)

            was intended to be a general function to calculate alpha_v
            using the concepts of [Garai2006] up to the Debye temperature
            and the lambda-model from [Raju2001] for higher temperatures

            after evaluation of the lambda-model: it is not exact enough to be used here

            The routine is used by the cp-calculations in the range 0 <= T <= Debye temperature ( based on [Garai2006] )
        """
        alpha_v=0 # for T=0 K
        if T>0 and T<=self.__debye_temperature:
            """ USED OPTION
                only [Garai2006] thermal expansion
            """
            alpha_v=self.get_alpha_v_Garai2006(T)
        elif T>self.__debye_temperature:
            """ OUTDATED OPTION
                lambda-concept
            """
            if cv_tmp==0:
                cv_tmp=self.get_debye_integral(T)
            # ok calculate
            p=-1/(T*B*V*self.__lambda)
            q=cv_tmp/(T*B*V)
            alpha_v=-0.5*p-m.sqrt(((p**2)/4)-q)
        # and back
        return alpha_v

    def get_volume(self,V0,T0,T,B,V,cv_tmp):
        """ OUTDATED ROUTINE (partially)

            general interface to calculate the volume between T_0 and T
            should (and is for cp calculation) only be used for temperatures <= Debye temperature

            warning: B is implemented as constant for lambda calculation ( for T > Debye temperature )
        """
        return V0*gmpy2.exp(quad(self.get_alpha_v,T0,T,args=(B,V,cv_tmp))[0])

    def get_volume_Raju2002b(self,H):
        """
            calculates the volume based on the volume-enthalpy relation by [Raju2002b]
            equation (28) in [Zienert2019]
        """
        return self.__V0+self.__Raju2002b_epsilon*H

    ####################
    # FITTING ROUTINES #
    ####################

    def fit_alpha_TD(self):
        """
            routine to fit alpha_v at the Debye temperature based on experimental volume data
            
            Option 1: parameter alpha_V^TD
            
            Option 2: parameter alpha_V^TD and V_0 --> set cp_variant = 'fit-alpha_TD,V0'
                        optimise also for V_0 (e.g. if it is based on DFT)

            Option 3: parameter alpha_V^TD and debye_temperature --> set cp_variant = 'fit-alpha_TD,TD'
                        aka 'find Debye temperature from volume instead of cp'
        """
        T=[]
        V=[]
        # ask for datafile
        str_datafile=input('Give datafile with experimental data,\nleave empty to manually input only one datapoint\n: ')
        if str_datafile: # read experimental data
            # columns
            df=datafile()
            columns=df.ask_property('columns') # ask which columns should be used
            del(df)
            #
            T_max=float(input('Give maximum temperature for fitting in K: '))
            with open(str_datafile,'r') as fp:
                for line in fp:
                    if line[0]!='#':
                        data=line.split()
                        if float(data[0])>T_max:
                            break
                        T.append(float(data[columns[0]]))
                        V.append(float(data[columns[1]])*1e-6)
        else: # manually input a datapoint (T,V) for fitting
            T.append(float(input('Give temperature in K: ')))
            V.append(float(input('Give V in cm^3/mol: '))*1e-6)
        #
        data_x=np.array(T)
        data_y=np.array(V)
        self.set_property('cp_calc',[1,data_x[len(data_x)-1]]) # set the start and end temperature for cp calculation
        #
        p=[50e-6] # first guess for alpha_V^TD
        #
        # adding more fitting parameter if necessary
        if self.__cp_variant=='fit-alpha_TD,V0':
            p.append(self.__V0)
        if self.__cp_variant=='fit-alpha_TD,TD':
            p.append(self.__debye_temperature)
        #            
        fit_res=leastsq(self.fit_alpha_TD_residuals,p,args=(data_x,data_y)) # do fitting

        return fit_res

    def fit_alpha_TD_residuals(self,p,x,y):
        """
            residuals routine for fitting alpha_V^TD from experimental volume data
            for calling, use the fit_alpha_TD() function
        """
        self.set_property('alpha_TD',p[0]) # set value of alpha_V^TD
        # depending on cp_variants, set other fitting parameters
        if self.__cp_variant=='fit-alpha_TD,V0':
            self.set_property('V0',p[1])
        if self.__cp_variant=='fit-alpha_TD,TD':
            self.set_property('debye_temperature',p[1])
        #
        self.set_property('epsilon',self.get_Raju2002b_epsilon()) # calculate epsilon [Raju2002b] from the given set of parameters
        #
        mhf=math_helpers()
        mhf.data_x,mhf.data_y=self.calc_cp() # do a cp calculation
        V=[]
        for i in x:
            V.append(mhf.get_y(i)) # x could be a float, cp_calc returns a set of integers in data_x --> mhf.get_y(i) returns a linearly interpolated value
        V=np.array(V)
        #
        # print the fitting results
        if self.__cp_variant=='V':
            print('alpha_TD: {0}'.format(p[0]*1e6))
        elif self.__cp_variant=='fit-alpha_TD,V0':
            print('alpha_TD: {0}, V0: {1}'.format(p[0]*1e6, p[1]*1e6))
        elif self.__cp_variant=='fit-alpha_TD,TD':
            print('alpha_TD: {0}, TD: {1}'.format(p[0]*1e6, p[1]))
        #
        return y-V # return residuals

    def Garai2006_residuals(self,p,x,y,T_min,T_max):
        """
            residuals function for fitting alpha_V data using the alpha_V-cv relation by [Garai2006]
        """
        mh = math_helpers()
        mh.data_x=x
        mh.data_y=y
        #
        alpha_TD=0
        self.set_property('debye_temperature',p[0]) # set the Debye temperature and calculate cv at T_D
        #
        # get alpha_V at T_D from experimental data
        if self.__debye_temperature < T_min:
            alpha_TD=1e6
        elif self.__debye_temperature >= T_min:
            alpha_TD=mh.get_y(self.__debye_temperature)
        #
        tmp=[]
        #  calculate for each T an alpha_V value
        for i in x:
            cv=self.get_debye_integral(i)
            tmp.append(alpha_TD*(cv/self.__cv_TD))
        tmp=np.array(tmp)
        #
        if self.__debye_temperature > T_max or self.__debye_temperature < 0:
            tmp=tmp+(self.__debye_temperature - T_max)*1e6 # the fitting routine should find its minimun within the given temperature boundary
        #
        return y-tmp

    def get_Raju2002b_epsilon(self):
        """ 
            fit epsilon [Raju2002b] from calculated enthalpy and volume
            using equation (28) in [Zienert2019]
        """
        data_H=[np.float64(0)]
        data_V=[np.float64(self.__V0)]
        mh=math_helpers()
        mh.data_x.append(0)
        mh.data_y.append(0)
        V_old=self.__V0
        T_old=0
        H_old=0
        # preliminary calculation of cp up to the Debye temperature, section 3.1.4 in [Zienert2019]
        for i in np.arange(1,self.__debye_temperature,1):
            cv_tmp=self.get_debye_integral(i)
            V_i=self.get_volume(V_old,T_old,i,NONE,NONE,cv_tmp)
            alpha_i=self.get_alpha_v_Garai2006(i)
            B=self.get_B(i,V_i)
            cp_i=cv_tmp+(i*(alpha_i**2)*V_i*B)
            mh.data_x.append(i)
            mh.data_y.append(cp_i)
            d_H=quad(mh.integrand,T_old,i,limit=100)[0]
            H_old+=d_H
            data_V.append(np.float64(V_i))
            data_H.append(np.float64(H_old))
            # save old
            T_old=i
            V_old=V_i
        #
        # now fit with linear equation
        res=mh.data_fit('linear',data_H,data_V)
        #
        if self.__cp_variant=='save-V(H)-TD': # save calculated H and V values, and quit
            with open(self.__savefile,'w') as fp:
                fp.write('# V(H) = {0} + {1}*H\n#\n#H[J]\tV[m^3]\n'.format(res[0],res[1]))
                for i in range(0,len(data_H)):
                    fp.write('{0}\t{1}\n'.format(data_H[i],data_V[i]))
            #
            print('V(H) data written to {0}\nHave a nice day old fellow!'.format(self.__savefile))
            exit(0)
        #
        return res[1]

    def debye_residuals(self,p,phases,x,y):
        """
            the residuals functions to obtain the Debye temperature from fitting cp data

            Option 1: phases=NONE --> only one phase

            Option 2: phases = [[phase_amount,Theta_D],[],..]
                      phase_amount of fit is 1-sum(phase_amount_i)

                      --> cp is sum of cp from different phases, with one unknown
        """
        #
        if not phases: # standard: single phase
            tmp=[]
            self.set_property('debye_temperature-only',p[0])
            for i in x:
                tmp.append(self.get_debye_integral(i))
            return y-np.array(tmp)
        else: # more than one phase
            tmp=[]
            for i in x:
                n=0
                cv=0
                for l in phases: # sum up values (n, c_V) of known phases
                    self.set_property('debye_temperature-only',l[1])
                    n=n+l[0]
                    cv=cv+l[0]*self.get_debye_integral(i)
                # and fit the residual amount
                self.set_property('debye_temperature-only',p[0])
                cv=cv+((1-n)*self.get_debye_integral(i))
                tmp.append(cv)
            #
            return y-np.array(tmp) # return residuals

    #####################################
    # Debye temperature, Debye function #
    #####################################

    def TD_Debye1912(self,molar_mass,density,K,nu):
        """
            calculates the Debye temperature from elastic constants
            using equation (23) from [Debye1912]
        """
        A=(self.__c._const_h/self.__c._const_k)*m.pow((9*self.__c._const_N)/(4*m.pi),1./3.)
        fnu=2*m.pow((2*(1+nu))/(3*(1-2*nu)),1.5)+m.pow((1+nu)/(3*(1-nu)),1.5)
        #
        return A*(1/m.pow(molar_mass,1./3.))*(1/m.pow(density,1./6.))*(1/m.pow(fnu,1./3.))*(1/m.sqrt((1/(K*1e-2))*1e-12))

    def TD_Anderson1963(self,molar_mass,density,K,nu,q):
        """
            calculates the Debye temperature from elastic constants
            using equations (1) and (3) from [Anderson1963]
        """
        G=(3*K*(1-(2*nu)))/(2*(1+nu))
        vs=m.sqrt(G/density) # shear sound velocity
        vl=m.sqrt((K+(4/3)*G)/density) # longitudinal sound velocity
        vm=m.pow((1./3.)*((2./m.pow(vs,3))+(1./m.pow(vl,3))),(-1./3.)) # averaged sound velocity
        #
        print('\nsound velocities:\n\tv_s= {0} m/s\n\tv_l= {1} m/s\n\n\tv_m={2} m/s\n'.format(vs*1e3,vl*1e3,vm*1e3))
        #
        return (self.__c._const_h/self.__c._const_k)*m.pow((3*q*self.__c._const_N*density)/(4*m.pi*molar_mass),1./3.)*vm,vm

    def get_debye_temperature(self,E_s):
        """
            calculates the Debye temperature from the specific energy relation
            equation (42) in [Zienert2019]
        """
        return self.__debye_fit_a+self.__debye_fit_b*E_s

    def debye_integrand(self,x):
        """
            integrand of the Debye function
        """
        return (gmpy2.exp(x)*(x**4))/((gmpy2.exp(x)-1)**2)

    def get_debye_integral(self,T):
        """
            calculates c_V using the Debye function [Debye1912] at a given temperature
        """
        return 9*self.__c._const_N*self.__c._const_k*((T/self.__debye_temperature)**3)*quad(self.debye_integrand,0,self.__debye_temperature/T)[0]

    ###################
    #     Finally     #
    # c_p calculation #
    ###################

    def calc_cp(self):
        """
            calculates the thermophysical properties as described in [Zienert2019]
            results are saved in a textfile
        """
        # init variables
        B_old=self.__B0
        T_old=0
        V_old=self.__V0
        H_old=0 # enthalpy
        H_est=0
        d_cp=[0,0]
        d_H_cp=0
        #
        # for saving last three values, see equation (37) in [Zienert2019]
        delta_cp_x=[0,0,0]
        delta_cp_y=[0,0,0]
        #
        # for saving the results
        cp_results={}
        #
        # necessary math helpers
        cv_mh=math_helpers()
        cv_mh.data_x=[0]
        cv_mh.data_y=[float(0)]
        mh=math_helpers()
        mh.data_x.append(0)
        mh.data_y.append(0)
        #
        # ok start
        array_T,self.__cp_delta_T=np.linspace(self.__cp_T_start,self.__cp_T_end,num=self.__cp_T_end/self.__cp_delta_T,retstep=True) # creates the temperature list for calculation
        for i in array_T: # here, i is the temperature not the index of the calculation step (as described in [Zienert2019], section 3.1.6)
            cv_tmp=self.get_debye_integral(i)
            cv_mh.data_x.append(i)
            cv_mh.data_y.append(cv_tmp)
            # calc V with Raju2002b_epsilon
            # --> H estimation with integral(cv)+d_cp
            cv_d_H_est=quad(cv_mh.integrand,T_old,i,limit=100)[0]
            #
            # estimate dcp from fit of last 3 d_cp values (if already calculated)
            if i >= self.__cp_T_start+3*self.__cp_delta_T and i <= self.__debye_temperature:
                res=mh.data_fit('linear+quad',delta_cp_x,delta_cp_y) # equation (37) in [Zienert2019]
                dcp=res[0]+res[1]*i+res[2]*i**2
                d_H_cp=0.5*self.__cp_delta_T*dcp # equation (35) in [Zienert2019]
            if i>self.__debye_temperature and len(delta_cp_x)==3:
                # delete the first entry
                delta_cp_x.pop()
                delta_cp_y.pop()
            if i>self.__debye_temperature:
                d_H_cp=0.5*self.__cp_delta_T*(delta_cp_y[0]+(delta_cp_y[0]-delta_cp_y[1])) # equation (35) in [Zienert2019]
            #
            H_est=H_old+cv_d_H_est+d_H_cp # equation (38) in [Zienert2019]
            #
            V_i=self.get_volume_Raju2002b(H_est) # use enthalpy estimation to calculate the volume
            #
            # and alpha_V
            slope=mh.get_slope(T_old,i,V_old,V_i) # equation (31) in [Zienert2019]
            alpha_i=slope/V_i # equation (30) in [Zienert2019]
            #
            B=self.get_B(i,V_i) # calculate bulk modulus at T = i, equation (26) in [Zienert2019]
            #
            cp_i=cv_tmp+(i*(alpha_i**2)*V_i*B)
            #
            # start with lambda calculation after dcp estimation
            # note: lambda calculation is kept for evaluation purpose
            if i>=self.__cp_T_start+4*self.__cp_delta_T:
                l_i=alpha_i/cp_i
            else:
                l_i=0
            #
            # enthalpy
            mh.data_x.append(i)
            mh.data_y.append(cp_i)
            H_i=quad(mh.integrand,T_old,i,limit=100)[0]
            H_old+=H_i
            #
            # save the last three values
            d_cp.insert(0,cp_i-cv_tmp) # d_cp=[new,old,old_old]
            d_cp.pop() # d_cp=[new,old,old]
            delta_cp_x.insert(0,i)
            delta_cp_x.pop()
            delta_cp_y.insert(0,float(d_cp[0]-d_cp[1]))
            delta_cp_y.pop()
            #
            # save calculated values in cp_results
            cp_results[i]=[V_i,alpha_i,l_i,B,cv_tmp,cp_i,H_old,delta_cp_y[0],H_est]
            # save old values
            T_old=i
            V_old=V_i
        #
        if self.__cp_variant=='savefile': # all cp results are saved to savefile
            with open(self.__savefile,'w') as fp:
                fp.write('# cp parameter:\n#\tV0 = {0} cm^3/mol\n#\tB0 = {1} GPa\n#\tB_m = {2} GPa/K\n#\tT_B0 = {3} K\n#\tTheta_D = {4} K\n#\talpha_TD = {5} 1e6/K\n#\tRaju2002b_epsilon = {6} cm^3/J\n#\n'.format(self.__V0*1e6,self.__B0*1e-9,self.__B_m*1e-9,self.__T_B0,self.__debye_temperature,self.__alpha_v_TD*1e6,self.__Raju2002b_epsilon*1e6))
                fp.write('#T[K]\tV[1e6*m^3]\talpha_v[1e6/K]\tlambda\tB[GPa]\tcv[J/K*mol]\tcp[J/K*mol]\tH[J/mol]\td_cp[J/K*mol]\tH_est[J/K*mol]\n')
                for i in sorted(cp_results.keys()):
                    #fp.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(i,V_i*1e6,alpha_i*1e6,l_i*1e6,B*1e-9,cv_tmp,cp_i,H_old,delta_cp_y[0],H_est))
                    fp.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(i,cp_results[i][0]*1e6,cp_results[i][1]*1e6,cp_results[i][2]*1e6,cp_results[i][3]*1e-9,cp_results[i][4],cp_results[i][5],cp_results[i][6],cp_results[i][7],cp_results[i][8]))
            #
        elif self.__cp_variant in ['V','fit-alpha_TD,V0','fit-alpha_TD,TD']: # only volume data is needed for calculation of residuals
            x=[]
            y=[]
            for i in sorted(cp_results.keys()):
                x.append(i) # T
                y.append(cp_results[i][0]) # V
            #
            return x,y
