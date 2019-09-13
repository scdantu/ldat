#!/usr/bin/python3
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sagar Khavnekar, Mugdha Dhurandhar, Sarath Chandra Dantu
    

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
import os,logging,datetime
from asd.WASD import WASD
from asd.PASD import PASD
logger="";
def main():    
    global logger
    '''
        logging stuff
    '''
    now         =   datetime.datetime.now()
    label       =   now.strftime("%Y-%m-%d-%H-%M-%S")
    log_file    =   'ASD-'+label+'.log'
    #logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG);
    logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO);
    logger=logging.getLogger('Main ASDA');
    WDIR="/media/scdantu/dsarath/PhD_Data/gottingen/TIM/amber99sb/InputGen/crys_sim/analysis/asda";
    DDIR="/media/scdantu/dsarath/PhD_Data/gottingen/TIM/amber99sb/InputGen/crys_sim/analysis/in/xtc/sys_xtc";
    geo_extr_pdb="ref.pdb"; asdType='WASD'
    label='1-WCS';     tpr_file="%s/1-WCS.tpr"%(DDIR);   xtc_file="%s/1-WCS-Sys.xtc"%(DDIR);     nframes=15001;    
    
    #ref_group="resname CA"; 
    #selection_group="(resid 21 and (name OD1 or name OD2)) or (resid 23 and name OG) or (resid 30 and (name OE1 or name OE2))";
    #selection_group="(resid 21 and (name OD1 or name OD2))";
    #selection_group="(resid 31 and (name OE1 or name OE2))";
    selection_group="(resid 23 and name OG)";
     
    os.chdir(WDIR);
    
    ref_group=[1690];     selection_group=['resname SOL and name OW'];  selection_group_dist=['5']
    object_asda="";
    
    if(asdType=='PASD'):
        object_asda=PASD(label,tpr_file,xtc_file,nframes,ref_group,selection_group);
    elif(asdType=='WASD'):
        object_asda=WASD(label,tpr_file,xtc_file,nframes,ref_group,selection_group);
    

    #ref_struct_group='resid 12';    ref_Struct_selection_group=['resid 1', 'resid 2', 'resid 3'];
    #object_ref_asda=PASD(folder,label,tpr_file,geo_extr_pdb,1,ref_struct_group,ref_Struct_selection_group)
    
if __name__ == '__main__':
    main()
    pass
