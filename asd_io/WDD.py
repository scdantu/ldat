#!/usr/bin/python
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sagar Khavnekar, Mudgha Dhurandhar, Sarath Chandra Dantu
    

     Copyright (C) 2018-2019 Sarath Chandra Dantu 

     This file is part of: Ligand Densities Analysis Tool for MD trajectories

     LDAT is a free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     LDAT is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with LDAT.  If not, see <http://www.gnu.org/licenses/>.

'''
import os,logging,timeit
import numpy as np
import asd.Utilities as asda_utils
import asd.FileIO as fileIO



class WDD(object):
    '''
    classdocs
    '''
    def __init__(self,out_Label,wdd_File):
        '''
        Constructor
        '''
        a=timeit.default_timer();
        self._logger=logging.getLogger('WDD')
        self._file_wdd=wdd_File;    self._out_label=out_Label;
        self._file_check();
        self._initialize();
        self.read_wdd()
        b=timeit.default_timer()
        self._logger.info('%-15s : %s (s)','WDD PT',b-a)
    
    def _initialize(self):
        self.DICT_NWATERS={};   
        self._time_stamp=0.0;   self._nwaters=0;
    def read_wdd(self):
        self._logger.info('%-15s : %s','Loading WDD...',self._file_wdd)
        data_nwaters="%12s%12s\n"%('dTime','nW');
        fOpen=open(self._file_wdd,'r')
        fLine=fOpen.readlines();
        for line in fLine:
            split_line=line.strip().split()
            if(line[0]=='#'):
                self._time_stamp=float(split_line[1])/1000.0; self._nwaters=int(split_line[2]); 
                data_nwaters+='%12.5f%12d\n'%(self._time_stamp,self._nwaters)
        fName='%s.nw'%(self._out_label);
        fileIO.saveFile(fName, data_nwaters)
            
    def _file_check(self):
        asda_utils.checkfile(self._file_wdd);
