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
#
# import
from class_math_helpers import math_helpers
#
class datafile:
    def __doc__(self):
        """
            This class collects inpout/output routines for working with datafiles.

            NOTE: Not all routines are used by the scripts of cp-tools. However, I hope they are helpful for someone.
        """
    def __init__(self):
        # classes
        self.__mh=math_helpers()
        # init of vars
        self.__datafile='' # filename of the datafile
        self.__multidata_datafile=[] # list of multiple datafile filenames
        self.__multidata_equal=True # are columns for reading data are the same for all multiple datafiles?
        self.__multidata_columns=[] # if the columns are not equal, they are saved here for each datafile
        self.__column_x=0 # standard column for x-data
        # configs
        self.__pct_all=True
        self.__pct_variant='delta' # ['std','delta']

    def set_property(self,variant,value):
        """
            set the variables of the class using this routine
        """
        if variant=='datafile': # sets the filename of the datafile
            self.__datafile=value
        elif variant=='equal_data': # equal columns to be read from the multiple datafiles? True/False
            self.__multidata_equal=value
        elif variant=='column_x': # standard column for x-data
            self.__column_x=value
        elif variant=='multidata_file': # append filenames to the list of multiple datafile filenames
            self.__multidata_datafile.append(value)
        elif variant=='reset_multidata_file': # resets the list of multiple datafile filenames
            self.__multidata_datafile=[]
        elif variant=='multidata_columns':
            self.__multidata_columns=value # set the columns to be read from the single datafiles of the multiple datafile filenames
        elif variant=='pct_all':
            self.__pct_all=value
        elif variant=='pct_variant': # set the variant for a percentage calculation
            self.__pct_variant=value

    def get_property(self,variant):
        """
            read some of the variables of the class
        """
        result=0
        if variant=='multidata_columns': # return the columns to be read from the single files in the list of multiple datafile filenames
            result=self.__multidata_columns
        #
        return result

    def ask_property(self,variant):
        """
            use this routine to set and/or ask necessary informations during running of the program by users input
        """
        result=0
        if variant=='equal_data': # read equal columns from multiple datafiles?
            equal=input('Equal datafiles? ')
            if equal in ('y','yes'):
                self.set_property('equal_data',True)
            else:
                self.set_property('equal_data',False)
            #
            result=self.__multidata_equal
        #
        elif variant=='columns': # which columns should be read ?
            tmp=input('Give column numbers [x:y]: ').split(':')
            result=[int(tmp[0])-1,int(tmp[1])-1]
        #
        elif variant=='savefile': # get the filename of the savefile
            result=input('Give savefile: ')
        #
        return result

    def get_column_number(self):
        """
            returns the amount of data columns of a datafile
        """
        fp=open(self.__datafile,'r')
        result=0
        for line in fp:
            if line[0]!='#':
                tmp_data=line.split()
                result=len(tmp_data)
                break
        #
        fp.close()
        return result

    def ask_column_number(self):
        """
            asks and sets the [x,y] column numbers to be read in case of multiple datfiles
        """
        columns=[]
        #
        for i in self.__multidata_datafile:
            tmp=input('Give [x:y] columns for {0}: '.format(i))
            columns.append([int(tmp.split(':')[0])-1,int(tmp.split(':')[1])-1])
            if self.__multidata_equal: # we only need to ask once if the data is to be read the same from all files
                break
        #
        self.set_property('multidata_columns',columns)
        #
        return columns

    #######################
    # Data handling stuff #
    #######################

    # extract cfg stuff from __datafile
    def get_cfg(self,variant):
        """
            extracts the config informations for cp-tools

            variant:    1) 'cfg_cp_name' --> reads extra options from PHASE.cfg files
                        2) 'cfg-cp-calc-name' --> reads the config for a single phase from the auto-calc config file 

            NOTE: This routine is maybe ported in future to use the cfg routines of the SciPy package.
        """
        fp=open(self.__datafile,'r')
        #
        if variant=='cfg_cp_name': # reads PHASE.cfg
            data_cfg={}
            # read
            for line in fp:
                if line[0]!='#':
                    tmp=line.split('=')
                    data_cfg[tmp[0]]=float(tmp[1]) # standard value=float
            # over and out
            fp.close()
            return data_cfg
        #
        elif variant=='cfg-cp-calc-name': # read cfg from cp-calc-auto config file
            name=input('Give phase name in cfg file: ')
            res={}
            for line in fp:
                if line[0]!='#':
                    tmp=line.split('\t')
                    if tmp[0] == name:
                        res['V0']=float(tmp[1])*1e-6
                        res['B0']=float(tmp[2])*1e9
                        res['B_m']=float(tmp[3])*1e9
                        res['TD_fit']=float(tmp[4])
                        res['TD_calc']=float(tmp[5])
                        res['E_s']=float(tmp[6])
                        res['alpha_TD']=float(tmp[7])*1e-6
                        res['T_m']=float(tmp[8])
            #
            fp.close()
            if not res:
                print('{0} not found in cfg file {1}'.format(name,self.__datafile))
                return 0
            else:
                return res

    # extract data from __datafile
    def get_data_values(self,variant,column,args):
        """
            general routine to read tab-separated data from self.__datafile

            variant:    'first,last' --> reads the first and last value saved in column 'column', args=NONE
                        
                        'last' --> similar as 'first,last', but only the last value
                        
                        'column' --> returns a list with the values of column 'column', args=NONE
                        
                        'x:y' --> returns a list [[x],[y]] with the columns given in 'columns', args=NONE
                                    if 'columns'='ask' --> the columsn numbers will be asked

                        'fitdata' --> returns [[x],[y]] for the columns of 'columns' for all x >= x_last-arg[0], if only the last part of a dataset is needed

            column: [x,y] with x and y are integers of the column numbers to be read

            args:   list of additional values depending on chosen 'variant'
        """
        result=0
        #
        if variant=='first,last':
            data_start=0
            data_end=0
            firstline=1
            tmp=[]
            fp=open(self.__datafile,'r')
            for line in fp:
                if line[0]!='#':
                    tmp=line.split('\t')
                    if firstline:
                        data_start=float(tmp[column])
                        firstline=0
            fp.close()
            data_end=float(tmp[column])
            result=[data_start,data_end]
        #
        elif variant=='last':
            fp=open(self.__datafile,'r')
            for line in fp:
                if line[0]!='#':
                    result=float(line.split('\t')[column])
            fp.close()
        #
        elif variant=='column':
            result=[]
            fp=open(self.__datafile,'r')
            for line in fp:
                if line[0]!='#':
                    result.append(float(line.split('\t')[column]))
            fp.close()
        #
        elif variant=='x:y':
            data_x=[]
            data_y=[]
            if column=='ask':
                print('FILE: {0}'.format(self.__datafile))
                column=self.ask_property('columns')
            #
            fp=open(self.__datafile,'r')
            for line in fp:
                if line[0]!='#':
                    tmp=line.split('\t')
                    data_x.append(float(tmp[column[0]]))
                    data_y.append(float(tmp[column[1]]))
            fp.close()
            result=[data_x,data_y]
        #
        elif variant=='fitdata':
            # column=[x,y]
            # args=[dx]
            data_end=self.get_data_values('last',column[0],0)
            data_x=[]
            data_y=[]
            fp=open(self.__datafile,'r')
            for line in fp:
                if line[0]!='#':
                    tmp=line.split('\t')
                    tmp_x=float(tmp[column[0]])
                    if tmp_x >= (data_end-args[0]):
                        data_x.append(tmp_x)
                        data_y.append(float(tmp[column[1]]))
            fp.close()
            result=[data_x,data_y]
        #
        return result

    def convert_to_tab(self):
        """
            as the name says, the routine converts self.__datafile to a one with tab-separated values
            asks for the name of the savefile
        """
        fp=open(self.__datafile,'r')
        savefile=self.ask_property('savefile')
        fp_save=open(savefile,'w')
        for line in fp:
            if line[0]=='#':
                fp_save.write(line)
            else:
                tmp=line.split()
                newline=''
                for i in tmp:
                    if not newline:
                        newline=newline+i
                    else:
                        newline=newline+'\t'+i
                #
                newline=newline+'\n'
                fp_save.write(newline)
        #
        fp.close()
        fp_save.close()

    ############################
    # Some datafile math stuff #
    ############################

    def calc_percentage(self):
        """
            calculates percentages between a basefile (self.__datafile) and percentage_file (self.__multidata_datafile)

            NOTE: It works only if all data have the same structure (self.__multidata_equal=True), which is the case for the class_debye.calc_cp() output.
        """
        # first save base-data
        data_base=[]
        fp=open(self.__datafile,'r')
        for line in fp:
            if line[0]!='#':
                line_data=line.split()
                if len(data_base)==0:
                    for i in line_data:
                        data_base.append([float(i)]) # creates a list array
                else:
                    i=0
                    for d in line_data:
                        data_base[i].append(float(d)) # saves the value from each column
                        i+=1
                #
            #
        #
        fp.close()
        #
        for multi_datafile in self.__multidata_datafile: # ok, now go through all multidata files
            fp=open(multi_datafile,'r')
            #
            line='#'
            while(line[0]=='#'): # ignore all comments
                line=fp.readline()
            #
            columns=len(line.split()) # first data line
            line=line.split()
            #
            # open savefile
            savefile=multi_datafile.split('.data')[0]+'-pct-'+self.__pct_variant+'.data' # the savefile name is based on the datafile name
            fp_save=open(savefile,'w')
            fp_save.write('# percentage calculation between:\n#\tBase-file: {0}\n#\tPct-file: {1}\n#\n'.format(self.__datafile,multi_datafile)) # write some infos
            #
            pos=0 # start with first dataline in data_base
            #
            for x in data_base[self.__column_x]: # iterate the x-values
                fp_save.write('{0}'.format(x)) # first column will be the x value
                #
                for i in range(0,columns): # index of y values
                    if i != self.__column_x: # ignore the x column
                        tmp_pct=self.__mh.get_percentage(self.__pct_variant,data_base[i][pos],float(line[i])) # calculate the percentage value
                        fp_save.write('\t{0}'.format(tmp_pct)) # and write it to savefile
                #
                fp_save.write('\n')
                pos+=1 # next index of the dataline of the data_base
                line=fp.readline().split() # read new line from percentage file
            #
            # close files
            fp_save.close()
            fp.close()
            # and print information
            print('PCT: [{0}] finished.'.format(savefile))

    def reduce_datapoints(self,data_x,data_y,err):
        """
            reduces the number of datapoints of a given dataset without smoothing the data
            
            NOTE: it is helpful if you need to plot the original experimental values without blowing up the size of the corresponding plotfile

            METHOD: all data pairs will be dropped that can be linearly interpolated between two points with an error <= 'err'
        """
        # vars
        res_x=[]
        res_y=[]
        res_x.append(data_x[0])
        res_y.append(data_y[0])
        #
        last=len(data_x)-1 # index of the last datapoint
        if not err:
            err=float(input('Give error [%]: '))
        # ok drop all datapoints that can be linearly interpolated by the start and the nth next point
        start=0
        while(True):
            if (last-start) > 1: # at least one datapoint must be between start and end
                for i in range(start+2,last+1):
                    tmp_y=[]
                    m=self.__mh.get_slope(data_x[start],data_x[i],data_y[start],data_y[i])
                    for n in range(start+1,i):
                        tmp_y.append(data_y[start]+m*(data_x[n]-data_x[start]))
                    go=True
                    # get last k value
                    k=0
                    for k in range(0,len(tmp_y)):
                        if abs(tmp_y[k]-data_y[start+1+k]) > abs(data_y[start+1+k]/100)*err: # can it be interpolated within the error boundary?
                            break
                    #
                    if k != len(tmp_y)-1: # if the last checked datapoint can not be interpolated within error boundary
                        # save the point of i-1 in the results list
                        res_x.append(data_x[i-1])
                        res_y.append(data_y[i-1])
                        start=i-1 # set the new start index
                        break
                    else:
                        if i==last:
                            start=last
                            break
                #
            # end conditions
            if (last-start) == 1:
                break
            if start == last:
                break
        # add last data point
        res_x.append(data_x[last])
        res_y.append(data_y[last])
        #
        return [res_x,res_y]

    def merge_datafiles(self,columns):
        """
            merge x-y files into one file with sorted x-values

            USAGE: give datafiles in self.__multidata_datafile and the respective column numbers to be read in columns together with self.__multidata_equal
        """
        data={}
        for i in range(0,len(self.__multidata_datafile)):
            fp=open(self.__multidata_datafile[i],'r')
            for line in fp:
                if line[0]!='#':
                    tmp_data=line.split()
                    tmp_x=0
                    tmp_y=0
                    if self.__multidata_equal:
                        tmp_x=float(tmp_data[columns[0][0]])
                        tmp_y=float(tmp_data[columns[0][1]])
                    else:
                        tmp_x=float(tmp_data[columns[i][0]])
                        tmp_y=float(tmp_data[columns[i][1]])
                    if tmp_x in data:
                        if 'list' in str(type(data[tmp_x])):
                            data[tmp_x].append(tmp_y)
                        else:
                            tmp=data[tmp_x]
                            data[tmp_x]=[tmp,tmp_y]
                    else:
                        data[tmp_x]=tmp_y
            #
            fp.close()
        # ok sort and save
        self.set_property('datafile',input('Savefile?: '))
        fp=open(self.__datafile,'w')
        header_x=input('Header x?: ')
        header_y=input('Header y?: ')
        fp.write('# merged data from :\n')
        for i in self.__multidata_datafile:
            fp.write('#\t{0}\n'.format(i))
        fp.write('#\n#{0}\t{1}\n'.format(header_x,header_y))
        # write data
        for i in sorted(data.keys()):
            if 'list' in str(type(data[i])):
                for d in data[i]:
                    fp.write('{0}\t{1}\n'.format(i,d))
            else:
                fp.write('{0}\t{1}\n'.format(i,data[i]))
        #
        fp.close()
        print('All data written to [{0}].'.format(self.__datafile))

    def merge_datafiles_all(self,header):
        """
            merge all datafiles into one savefile by simply copying each dataline of each file (in self.__multidata_datafile) into the savefile (self.__datafile)
        """
        fp_save=open(self.__datafile,'w')
        if not header:
            header=input('Give header #[header]: ')
        fp_save.write('#{0}\n'.format(header))
        for i in self.__multidata_datafile:
            fp=open(i,'r')
            for line in fp:
                if line[0]!='#':
                    fp_save.write(line)
            fp.close()
        fp_save.close()

    ####################
    ## Some misc stuff #
    ####################
    
    def create_latex_table(self):
        """
            creates the LaTex code for a table filled with your x-y data saved in self.__datafile for a given amount of columns
            asks for the filename of the savefile
        """
        savefile=input('Give savefile: ')
        header_x=input('Give header for x: ')
        header_y=input('Give header for y: ')
        # read data
        data=self.get_data_values('x:y','ask',0)
        print('{0} datapoints found\n'.format(len(data[0])))
        spalten=int(input('How many columns? '))
        zeilen=round(len(data[0])/spalten)
        if zeilen-(len(data[0])/spalten)<0:
            zeilen+=1
        print('Amount of lines: {0}\n'.format(zeilen))
        # create latex stuff
        fp=open(savefile,'w')
        fp.write(r'\b'+'egin{table}\n\t\centering\n')
        fp.write('\t\caption{Put your caption here.}\n')
        # tabular
        tmp_str='\t'+r'\b'+'egin{tabular}{'
        for i in range(0,spalten,1):
            tmp_str+='cc'
        tmp_str+='}\n'
        fp.write(tmp_str)
        fp.write('\t\t'+r'\t'+'oprule\n')
        tmp_str='\t'
        # header
        for i in range(0,spalten,1):
            if i==0:
                tmp_str+='\t'+r'\t'+'extbf{'+'{0}'.format(header_x)+'}\t&'+r'\t'+'extbf{'+'{0}'.format(header_y)+'}'
            else:
                tmp_str+='\t&'+r'\t'+'extbf{'+'{0}'.format(header_x)+'}\t&'+r'\t'+'extbf{'+'{0}'.format(header_y)+'}'
        tmp_str+=r'\\'+'\n'
        fp.write(tmp_str)
        #
        fp.write('\t\t\midrule\n')
        # nun tabelle ausgeben
        for i in range(0,zeilen,1):
            tmp_str=''
            for l in range(0,spalten,1):
                index=((l*zeilen)+i)
                #print(index)
                if index>=len(data[0]):
                    tmp_str=tmp_str+'\t&\t&'
                    print('empty field [{0}]'.format(index))
                else:
                    if l==0:
                        tmp_str=tmp_str+'\t\t{0}\t&{1}'.format(int(float(data[0][index])),round(float(data[1][index]),2))
                    else:
                        tmp_str=tmp_str+'\t&{0}\t&{1}'.format(int(float(data[0][index])),round(float(data[1][index]),2))
            tmp_str=tmp_str+r'\\'+'\n'
            fp.write(tmp_str)
        #
        fp.write('\t\t'+r'\b'+'ottomrule\n\t\end{tabular}\n\t\label{tbl:label}\n\end{table}')
        fp.close()
        print('\nLaTex code written to {0}!'.format(savefile))
