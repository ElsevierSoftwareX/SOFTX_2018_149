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

    Description: starts the calculation of thermophysical properties using the algorithm described in [Zienert2018], see class_debye

    USAGE: set the variables in the config section below and call the script using 'python cp-calc-auto.py'
"""
#
# imports
import sys
sys.path.append('classes')
import os
import pathlib
from concurrent.futures import ProcessPoolExecutor, wait # for parallel computing
# own
from class_debye import debye
from class_datafile import datafile
#
#########
#

#
# config # change to your needs
#
cp_variant='free' # ['free','prediction']
TD_variant='fit' # ['fit','calc','prediction']
#
m = -0.0204988*1e9 # mean Bulk modulus in Pa


#PCT_vars=['alpha','B0','m','V0','TD']
PCT_vars=[]

PCT_percents=[-10,-5,-2,-1,0,1,2,5,10]
#PCT_percents=[0]

CFG_compute_only_new=True # [True,False] , checks if the saving file already exists
CFG_n_cpu = 4 # number of your CPUs that should be used for parallel computing
#

cfg_dir='cfg-cp-calc/' # folder with config files
cfg_file='data-cp-calc.cfg'

datafile=cfg_dir+cfg_file # config file for auto-cp calculation
#                           structure datafile (tab separated): name,V0[cm^3/mol],B0[GPa],m[GPa/K],Theta_fit[K],Theta_calc[K],E_s[KJ/g],alpha_D[1e6/K],T_m[K]
#

file_PREFIX='cp-calc/' # folder where the calculated data will be saved in

#
##########
#

#
# functions that are called by ProcessPoolExecutor #
#

def cp_calc(V0,B0,m,T_D,alpha_TD,savefile,name,pct,var,value):
    """
        routine that actually starts the cp-calculation for each phase with the respective PHASE config
    """
    #
    # print information
    if var:
        print('starting: '+name+'-{0}_{1}_'.format(var,pct)+' , {0}'.format(value))
    else:
        print('starting: '+name)
    #
    from class_debye import debye
    debye=debye(T_D) # set Debye temperature
    debye.set_cp_parameters_auto(V0,B0,m,T_D,alpha_TD,'savefile',savefile) # set all other parameters
    #
    # look for extra cfg
    if name+'.cfg' in files_cfg:
        from class_datafile import datafile
        df=datafile()
        df.set_property('datafile',cfg_dir+name+'.cfg') # set the filename
        cfg_data=df.get_cfg('cfg_cp_name') # and read the cfg data
        #
        cfg_keys=cfg_data.keys()
        if 'cp_T_end' in cfg_keys: # set a new value for the end temperature of cp-calculation, std: T = 2000 K
            debye.set_property('cp_T_end',cfg_data['cp_T_end'])
            print('{0}: Setting cp_T_end = {1}\n'.format(name,cfg_data['cp_T_end']))
        if 'epsilon' in cfg_keys: # NOTE: testing purpose only, set an arbitary epsilon value
            debye.set_property('epsilon',cfg_data['epsilon']*1e-6)
    #
    # ok start cp calculation
    debye.calc_cp() # results will be written into 'savefile'

def pct_calc(base_file,pct_file,pct_variant,column_x):
    """
        calculates the percentage (here the difference) of one dataset (pct_file) with another (base_file)
        it is used here to calculate the influence of the single cp-parameters on the calculated thermophysical properties
    """
    from class_datafile import datafile
    df=datafile()
    df.set_property('pct_variant',pct_variant) # sets the variant, probably 'delta' is here meaningful
    df.set_property('datafile',base_file) # data to compare with
    df.set_property('multidata_file',pct_file) # data for percentage calculation
    df.set_property('column_x',column_x) # sets the index of the x value (here of the temperature)
    df.calc_percentage() # start the calculation

#
####################################################
#

# init
debye=debye(1)
debye.set_property('alpha_v-fit',[233.847471388462,-191.509884803456,-235.603810791667,0.00195317441266783]) # set the paramters for alpha_V-T_m relation
debye.set_property('debye-fit',[38.9117,14.5631]) # set the paramters for T_D-E_s relation

PCT_files=[]
#
T_D=0
#
# check if the directories exist, and if not create them (needs python 3.5+)
pathlib.Path(file_PREFIX).mkdir(parents=True, exist_ok=True)
pathlib.Path(cfg_dir).mkdir(parents=True, exist_ok=True)
#
# list files of the results and cfg directories
files_result=os.listdir(file_PREFIX)
files_cfg=os.listdir(cfg_dir)

# multiprocessor stuff
pool=ProcessPoolExecutor(CFG_n_cpu)
futures=[]



fp=open(datafile,'r') # open the config file for automatic cp-calculation
n=0
for l in fp:
    if l[0]!='#': # ingore comments
        data=l.split('\t') # we have tab separated values
        #
        # get the respective values
        name=data[0] # name of the calculations, is part of the savefile
        V0=float(data[1])*1e-6 # volume at T = 0 K
        B0=float(data[2])*1e9 # bulk modulus at T = 0 K
        if data[3]!='0':
            m=float(data[3])*1e9 # mean linear slope of B(T)
        T_D_fit=float(data[4]) # calorimetric Debye temperature (fit from cp)
        T_D_calc=float(data[5]) # calculated Debye temperature based on elastic properties
        E_s=float(data[6]) # specific energy needed for TD_variant='prediction'
        alpha_TD=float(data[7])*1e-6 # alpha_V at T = T_D
        T_m=float(data[8]) # melting temperature
        #
        ## set the value of the Debye temperature depending on TD_variant
        if TD_variant=='calc': # using calculated T_D from elastic properties
            T_D=T_D_calc
        elif TD_variant=='fit': # using calorimetric Debye temperature
            T_D=T_D_fit
        elif TD_variant=='prediction': # calculate the Debye temperature from specific energy
            T_D=debye.get_debye_temperature(E_s)
        #
        # set the value of alpha_V at T = T_D depending on cp_variant
        if cp_variant=='prediction': # calculate alpha_V at T = T_D from melting temperature
            alpha_TD=debye.get_alpha_v_fit(T_m)
        #
        if len(PCT_vars)>0: # for evaluation of parameters
            old_value=0
            for var in PCT_vars:
                # set start values
                parameters={
                    'alpha':{'tmp':alpha_TD,'old':alpha_TD},
                    'B0':{'tmp':B0,'old':B0},
                    'm':{'tmp':m,'old':m},
                    'V0':{'tmp':V0,'old':V0},
                    'TD':{'tmp':T_D,'old':T_D}
                }
                #
                # calculate one percent of the percentage parameter
                old_value=parameters[var]['old']
                one_pct=old_value/100
                #
                tmp_files=[] # filenames where the data is save that will be compared with the basefile data
                tmp_base='' # filename of the basefile datafile
                #
                for i in PCT_percents: # iterate over all percentages
                    new_value=old_value+i*one_pct # calculate the new value of the percentage parameter
                    parameters[var]['tmp']=new_value # and save it
                    #
                    if i==0:
                        savefile=file_PREFIX+'cp-calc-'+cp_variant+'-'+TD_variant+'-'+name+'.data' # basefile
                    else:
                        savefile=file_PREFIX+'cp-calc-'+cp_variant+'-'+TD_variant+'-'+name+'-{0}_{1}_'.format(var,i)+'.data' # percentage file
                    #
                    if CFG_compute_only_new: # do not recalculate if the results file exists
                        if savefile.split('/')[len(savefile.split('/'))-1] not in files_result:
                            if i==0:
                                tmp_base=savefile
                            else:
                                tmp_files.append(savefile)
                            # ok, start cp-calc
                            futures.append(pool.submit(cp_calc,parameters['V0']['tmp'],parameters['B0']['tmp'],parameters['m']['tmp'],parameters['TD']['tmp'],parameters['alpha']['tmp'],savefile,name,i,var,new_value))
                    else: # calculate all and overwrite existing results
                        if i==0:
                            tmp_base=savefile
                        else:
                            tmp_files.append(savefile)
                        # ok, start cp-calc
                        futures.append(pool.submit(cp_calc,parameters['V0']['tmp'],parameters['B0']['tmp'],parameters['m']['tmp'],parameters['TD']['tmp'],parameters['alpha']['tmp'],savefile,name,i,var,new_value))
                #
                if tmp_base:
                    PCT_files.append([tmp_base,tmp_files]) # save the set of percentage files for later percentage calculation
        #
        else: # no variation of cp parameters
            savefile=file_PREFIX+'cp-calc-'+cp_variant+'-'+TD_variant+'-'+name+'.data'
            ## now, we can start :-)
            #
            if CFG_compute_only_new:
                if savefile.split('/')[len(savefile.split('/'))-1] not in files_result: # check if the savefile already exists
                    futures.append(pool.submit(cp_calc,V0,B0,m,T_D,alpha_TD,savefile,name,0,'',0))
            else: # we do all calculations and overwrite existing results
                futures.append(pool.submit(cp_calc,V0,B0,m,T_D,alpha_TD,savefile,name,0,'',0))
#
fp.close() # close cfg file
wait(futures) # wait until all calls of cp-calc are finished
#
# finally, compute percentage (if any)
for i in PCT_files:
    if len(i[1])>0:
        for files in i[1]:
            futures.append(pool.submit(pct_calc,i[0],files,'delta',0))
#
wait(futures)
#
print('All done, have a nice day old fellow!')
