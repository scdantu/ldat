#!/usr/bin/python
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sagar Khavnekar, Sarath Chandra Dantu
    

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
import os,logging,math,timeit
import numpy as np

'''
    matplotlib 
'''
import matplotlib as mplot
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as P 
'''
'''
import asd.Utilities as asda_utils
import asd.FileIO as fileIO

class WDD_Plot(object):
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
        self.X_COORDS=[];   self.Y_COORDS=[];   self.Z_COORDS=[];
    def read_wdd(self):
        self._logger.info('%-15s : %s','Loading WDD...',self._file_wdd)
        fOpen=open(self._file_wdd,'r')
        fLine=fOpen.readlines();
        for line in fLine:
            split_line=line.strip().split()
            if(line[0]!='#'):
                self.X_COORDS.append(float(split_line[2]));
                self.Y_COORDS.append(float(split_line[3]));
                self.Z_COORDS.append(float(split_line[4]));
        self._plot_positions()
    def _plot_positions(self):
        
        #NOTE: plot function should read a file and plot iteratively like other plotting scripts
        #NOTE: advised that it be called in main/do_analysis
        fig = P.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("equal")
        myc=['#3BB9FF','#4E9258','#FBB917','#1657FA','#FF00FF','#0000FF','#FF0000','#FF88FF'];
        N=len(self.X_COORDS);
        self._logger.info('Found total %d waters',N)
        for i in range(N):
            ax.scatter(self.X_COORDS[i],self.Y_COORDS[i],self.Z_COORDS[i],c='#000000')
            
        #ax.scatter(x1, y1, z1,c="#FF0000")
        
        xmin=-4;xmax=4
        ymin=-4;ymax=4
        zmin=-4;zmax=4
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_zlim(zmin,zmax)
            
        #draw center sphere : How to specify the sphere radius?
        u,v = np.mgrid[0:2*np.pi:12j,0:np.pi:5j]
        mg_x=np.cos(u)*np.sin(v)
        mg_y=np.sin(u)*np.sin(v)
        mg_z=np.cos(v)
        ax.plot_wireframe(mg_x,mg_y,mg_z,color='#FF0000')
        #P.legend(bbox_to_anchor=(0., 1.02, 1., .102),prop=mplot.font_manager.FontProperties(size=12),frameon=False,loc=3,ncol=3,mode="expand")
        #P.legend(prop=mplot.font_manager.FontProperties(size=16),frameon=False,ncol=3,mode="expand")
        
        #print notification total number of points ploted
        
        
        P.show()
        #save="Pos-%s.png"%(self._out_label)
        #P.savefig(save,dpi=600)
        exit()
    
    def _get_densities(self,in_x,in_y,in_z,binsize):
        #Calculates densities from 3D data
        #Range should also be defined by input parameter file
        i=0
        input_coords=np.random.randn(len(in_x),3)
        
        for i in range (len(in_x)):
            input_coords[i,0]=in_x[i]
            input_coords[i,1]=in_y[i]
            input_coords[i,2]=in_z[i]
        
        H,edges = np.histogramdd(input_coords, bins=(binsize,binsize,binsize),range=[[-4,4],[-4,4],[-4,4]])
        edges1 = np.delete(np.array(edges[0]),H.shape[0])
        edges2 = np.delete(np.array(edges[1]),H.shape[1])
        edges3 = np.delete(np.array(edges[2]),H.shape[2])
        x = edges1
        y = edges2
        z = edges3
        return H, x, y, z
    
    def _check_densities(self,H_in,x_in,y_in,z_in):
        CX=[]; CY=[]; CZ=[]; CH=[]
       
        i= 0
        j=0
        k=0
        for i in range(len(x_in)):
            for j in range(len(y_in)):
                for k in range(len(z_in)):
                    if H_in[i,j,k] > 0:
                        CH.append(H_in[i,j,k])
                        CX.append(x_in[i])
                        CY.append(y_in[j])
                        CZ.append(z_in[k])
                        #print x_in[i], "   ",y_in[j],"    ", z_in[k],"    ", H_in[i,j,k]
    
        #print CX,CY,CZ, len(CH), len(x_in)
        fig = P.figure()
        ay = fig.gca(projection='3d')
        ay.set_aspect("equal")
        
        i=0
        '''
        for i in range(len(CX)):
            print CX[i], "   ",CY[i],"    ", CZ[i],"    ", CH[i]
            ay.scatter(CX[i], CY[i], CZ[i],c="#0000FF")
        '''
        '''
        xmin=-4;xmax=4
        ymin=-4;ymax=4
        zmin=-4;zmax=5
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_zlim(zmin,zmax)
        '''
        #draw MG sphere : How to specify the sphere radius?
        u,v = np.mgrid[0:2*np.pi:24j,0:np.pi:10j]
        mg_x=np.cos(u)*np.sin(v)
        mg_y=np.sin(u)*np.sin(v)
        mg_z=np.cos(v)
        ay.plot_wireframe(mg_x,mg_y,mg_z,color='#FFA500')
        P.show()
    
    def _convert_densities_to_xyzv(self,H_in,x_in,y_in,z_in):
        #converts densities to x,y,z,value format
        i=0; j=0; k=0;
        value =[]; x_out = []; y_out = []; z_out = []
        for i in range(len(x_in)):
            for j in range(len(y_in)):
                for k in range(len(z_in)):
                    value.append(H_in[i,j,k])
                    x_out.append(x_in[i])
                    y_out.append(y_in[j])
                    z_out.append(z_in[k])
                    #if H_in[i,j,k] > 5:
                        #print x_in[i], "  ", y_in[j], "  ", z_in[k], "  ", H_in[i,j,k], "\n"
        return x_out, y_out, z_out, value


    def _write_gaussCUBE(self,H_in,x_in,y_in,z_in):
        #write gaussian cube formate
        #Write the comment/free text
        cubefile = open("density.cube", "w")
        cubefile.write("Active site densities calculated from ---- program\n")
        cubefile.write("Units are in Angstorms\n")
        #write num of atoms and the origin
        num_atoms = 1
        ooc = [x_in[0]/0.529177249, y_in[0]/0.529177249,z_in[0]/0.529177249]
        cubefile.write("%5i%18.6e%18.6e%18.6e\n" % (num_atoms,ooc[0],ooc[1],ooc[2]))
        #write number of axis info 
        x_info = [len(x_in),abs((x_in[0] - x_in[1])/0.529177249),0.0,0.0]
        y_info = [len(y_in),0.0,abs((y_in[0] - y_in[1])/0.529177249),0.0]
        z_info = [len(z_in),0.0,0.0,abs((z_in[0] - z_in[1])/0.529177249)]
        cubefile.write("%5i%18.6e%18.6e%18.6e\n" % (x_info[0],x_info[1],x_info[2],x_info[3]))
        cubefile.write("%5i%18.6e%18.6e%18.6e\n" % (y_info[0],y_info[1],y_info[2],y_info[3]))
        cubefile.write("%5i%18.6e%18.6e%18.6e\n" % (z_info[0],z_info[1],z_info[2],z_info[3]))
        #write atomic coordinates 
        atomic_number = 8; zero = 0.0; c = [0,0,0];
        cubefile.write("%5i%4.1f%18.6e%18.6e%18.6e\n" % (atomic_number,zero,c[0],c[1],c[2]))
        i=0; j=0; k=0;
        values =[]
        for i in range(len(x_in)):
            for j in range(len(y_in)):
                for k in range(len(z_in)):
                    value = H_in[i,j,k]
                    
                    cubefile.write("%18.6e" %value)
                    values.append(value)
                    
                    if (k %100 == 99):
                        #print k % 10
                        cubefile.write("\n")
                cubefile.write("\n")
        cubefile.close()
        self._logger.info('Wrote Densities as gaussian cube format')
                    #if H_in[i,j,k] > 5:
                        #print x_in[i], "  ", y_in[j], "  ", z_in[k], "  ", H_in[i,j,k], "\n"
        #print len(values)                
        #print ("\t".join(map(str,values)))               


        
        
    

    def _file_check(self):
        asda_utils.checkfile(self._file_wdd);
