'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sarath Chandra Dantu
    

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
import os,logging
logger=logging.getLogger('FileIO')
def saveFile(fName,data):
    global logger
    f=open(fName,"w")
    f.write(data)
    logger.info('Saved file %s',fName)
def checkfile(fileName):
    global logger
    if(os.path.isfile(fileName))==False:
        logger.error('%s does not exist',fileName)
        exit();
    else:
        logger.debug('%s exists',fileName)
        
def checkfile_with_return(fileName):
    return os.path.isfile(fileName)
def checkfile_with_message(fileName,message):
    global logger
    if(os.path.isfile(fileName))==False:
        logger.info('%s does not exist. %s',fileName,message)
        exit();