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
from MDAnalysis.analysis.align import *
#from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mdana  
from MDAnalysis.lib.mdamath import triclinic_vectors as trivec
from numpy.distutils.from_template import list_re
#from tkinter.constants import ACTIVE
'''
    python defautlt
'''
import os,sys,logging,math,timeit
import numpy as np
from numpy import matrix
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

#from fileIO import Parse_Options as parse_options
#import MDAnalysis.coordinates as coord
#from MDAnalysis.core.Timeseries import Timeseries
#import fileIO.Write_1dData as w1d


class WASD(object):
    '''
    classdocs
    '''
    def __init__(self,out_label,tpr_file,trajectory_file,frames,ref_group,sel_group="resname SOL and name OW",cut_off=4.0,rotateM=False):
        '''
        Constructor
        '''
        self._logger=logging.getLogger('WASD')
        self._file_tpr=tpr_file;         self._file_trajectory=trajectory_file; self._out_label=out_label;
        '''
            tune these two variables to select which portions of trajectory are analyzed
        '''        
        self._first_frame=0;                self._frame_limit=frames;       
        self._reference_group=ref_group;    self._selection_group=sel_group;
        self._universe="";                  self._uni_ref="";               self._box="";
        
        self._cut_off=cut_off;
        
        
        self._file_checks();
        self._initialize();
        self._process_selections()
        self._analysis_manager();
    
    def _initialize(self):
        a=timeit.default_timer();
        self._list_frame_data=[];   self._str_frame_data="";
        self._logger.info('%-15s : %s','Loading XTC...',self._file_trajectory)
        self._universe  =   mdana.Universe(self._file_tpr,self._file_trajectory);
        #list_resids=list(self._universe.atoms.resids)
        #print(list_resids[0],list_resids[-1],len(list_resids))
        b=timeit.default_timer()
        self._logger.info('%-15s : %s (s)','XTC LT',b-a)
        self._box       =   self._universe.dimensions;
        self._nframes   =   len(self._universe.trajectory);
        self._dt        =   self._universe.trajectory.dt;
        self._sim_len   =   self._dt*self._nframes
        self._logger.info('%-15s : %s','box dim.',self._box);
        self._logger.info('%-15s : %s','nframes',self._nframes);
        self._logger.info('%-15s : %s','dt',self._dt);
        self._logger.info('%-15s : %s (ps)','Length',self._sim_len);
        if(self._frame_limit>self._nframes):
            self._logger.info('No. of frames requested for Analysis : %8d > nframes. Setting frame limit to nframes',self._frame_limit);
            self._frame_limit=self._nframes;
        else:
            self._logger.info('%-15s : %d ','Num frames',self._frame_limit);
        
        '''
            Since a pbc trajectory is not being used, we need a protein trajectory with pbc corrected to calculate the rotation matrix
            #self._uni_ref   =   Universe(self._file_tpr,trjfile_ref)
            for CaMD since the xtc is pbc converted, same trajectory file can be used
        '''
        
    def _analysis_manager(self):
        for current_selection in self._list_of_selections:
            self._logger.info('%-15s : %s','Processing',current_selection[0])
            self._get_active_site(current_selection[0],current_selection[1],current_selection[2]);
            
    def _process_selections(self):
        self._list_of_selections=[];    ref_group=[];
        self._logger.info('%-15s : %s','Num Selections',len(self._reference_group))
        for resid in self._reference_group:
            reference_group_selection='resid %d and name N'%(resid-1)
            #ref_group.append(ref);
            current_selection='(%s) and around %f (%s)'%(self._selection_group[0],self._cut_off,reference_group_selection)
            self._list_of_selections.append([resid,reference_group_selection,current_selection]);
        
    def _get_active_site(self,residue,reference_selection,current_selection):
        #gets density of selected atoms in the active site
        '''
            RESIDUE NUMBERING IS RE-INTIALIZED TO ZERO by MDANALYSIS
            WE NEED TO WORK AROUND THIS
        '''
        frameCount=0;
        for ts in self._universe.trajectory[self._first_frame:self._frame_limit]:
            self._current_time=frameCount*self._dt;
            print_t="WASD - INFO - %-15s : %.2f/%.2f"%('Frame No.',self._current_time,self._dt*self._nframes)
            sys.stdout.write("\r{0}".format(print_t))
            sys.stdout.flush()
            #ts.frame counter works only if the trajectory is analysed from the first frame
            #rm = R[ts.frame]
            self._ref_atom      =   self._universe.select_atoms(reference_selection);
            self._active_atoms  =   self._universe.select_atoms(current_selection);
            #list_of_atoms=active_atoms.atoms.unique.names;
            #list_of_atoms=self._active_atoms.atoms.names;
            #list_of_resid=self._active_atoms.residues
            '''
                this should be checked before this function
            '''
            if(len(self._active_atoms)>0):
                self._calculate_distances();
            self._store_data();
            frameCount+=1;
        self._logger.info('')    
        fName='%s-%d.wdd'%(self._out_label,residue)
        fileIO.saveFile(fName, self._str_frame_data)
        self._str_frame_data="";
    def _calculate_distances(self):
        self._list_frame_data=[];
        for active_atom in self._active_atoms:
            ref_pos = self._ref_atom.positions[0];
            active_pos = active_atom.position;
            (i_dist,x,y,z) = self._get_pbc_dist(ref_pos,active_pos)
            d=[(active_atom.resid)+1,i_dist,x,y,z]
            self._list_frame_data.append(d);
        
    def _store_data(self):
        if(len(self._list_frame_data)>0):
            self._str_frame_data+="#TSTAMP %20.5f%12d\n"%(self._current_time,len(self._list_frame_data))
            for iData in self._list_frame_data:
                self._str_frame_data+="%12d%12.5f%12.5f%12.5f%12.5f\n"%(iData[0],iData[1],iData[2],iData[3],iData[4])
        else:
            self._str_frame_data+="#TSTAMP %20.5f%12d\n"%(self._current_time,0)
        self._list_frame_data=[];
        
            
    def _box_check(self):
        """Take a box input and deduce what type of system it represents based
        on the shape of the array and whether all angles are 90.
        :Arguments:
          *box*
              box information of unknown format
        :Returns:
          * ``ortho`` orthogonal box
          * ``tri_vecs`` triclinic box vectors
          * ``tri_box`` triclinic box lengths and angles
        :Raises:
          * ``TypeError`` if box is not float32
          * ``ValueError`` if box type not detected
        """
        if self._box.dtype != np.float32:
            raise TypeError("Box must be of type float32")
    
        boxtype = 'unknown'
    
        if self._box.shape == (3,):
            boxtype = 'ortho'
        elif self._box.shape == (3, 3):
            if np.all([self._box[0][1] == 0.0,  # Checks that tri box is properly formatted
                          self._box[0][2] == 0.0,
                          self._box[1][2] == 0.0]):
                boxtype = 'tri_vecs'
            else:
                boxtype = 'tri_vecs_bad'
        elif self._box.shape == (6,):
            if np.all(self._box[3:] == 90.):
                boxtype = 'ortho'
            else:
                boxtype = 'tri_box'
    
        if boxtype == 'unknown':
            raise ValueError("box input not recognised"
                             ", must be an array of box dimensions")
    
        return boxtype
    
    def _triclinic_pbc(self,incoords,inbox,box_inverse):
        localcoords = incoords.copy()
        # z
        s = math.floor(localcoords[2] * box_inverse[2]);
        localcoords[2] -= s * inbox[2][2];
        localcoords[1] -= s * inbox[2][1];
        localcoords[0] -= s * inbox[2][0];
        # y
        s = math.floor(localcoords[1] * box_inverse[1]);
        localcoords[1] -= s * inbox[1][1];
        localcoords[0] -= s * inbox[1][0];
        # x
        s = math.floor(localcoords[0] * box_inverse[0]);
        localcoords[0] -= s * inbox[0][0];
        #print incoords, localcoords
        return localcoords
    
    def _get_pbc_dist(self,a,b):
        
        a_coord = a.copy()
        b_coord = b.copy()
        box = trivec(self._box)
        #print a_coord, a_coord[0] 
        
        dr = a_coord - b_coord
    
        r = np.linalg.norm(dr)
       
        #calculate box half and box inverse
        box_inv = np.zeros((3), dtype=np.float32)
        box_half= np.zeros((3), dtype=np.float32)
        
        box_inv[0] = 1.0 / box[0][0]
        box_inv[1] = 1.0 / box[1][1]
        box_inv[2] = 1.0 / box[2][2]
      
        box_half[0] = 0.5 * box[0][0];
        box_half[1] = 0.5 * box[1][1];
        box_half[2] = 0.5 * box[2][2];
        
        #get all coordinates inside the box
        new_a_coords=self._triclinic_pbc(a_coord, box, box_inv)
        new_b_coords=self._triclinic_pbc(b_coord, box, box_inv)
        
        #calculate distances considering PBC
        dx1=new_a_coords-new_b_coords
        dx=dx1.copy()
        #
        if (math.fabs(dx[2]) > box_half[2]):
            if (dx[2] < 0.0 ):
                dx[2] += box[2][2];
                dx[1] += box[2][1];
                dx[0] += box[2][0];
            else:
                dx[2] -= box[2][2];
                dx[1] -= box[2][1];
                dx[0] -= box[2][0];
        #y
        if (math.fabs(dx[1]) > box_half[1]):
            if (dx[1] < 0.0):
                dx[1] += box[1][1];
                dx[0] += box[1][0];
            else:
                dx[1] -= box[1][1];
                dx[0] -= box[1][0];
        #x
        if (math.fabs(dx[0]) > box_half[0]):
            if (dx[0] < 0.0):
                dx[0] += box[0][0];
            else:
                dx[0] -= box[0][0];
        
        d = np.linalg.norm(dx)
        x=dx[0];        y=dx[1];        z=dx[2];
        
        #return the normal value and differences in X Y and Z dimensions
        return d, x, y, z
    
    def _get_pos(self,input1):
        type1 = str(type(input1))
        if "AtomGroup.AtomGroup" in type1:
            pos = input1.positions[0]
        else:
            pos = input1.position
        return pos
    
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

    def _dump_shift_PDB(self):
        self._uni_ref.trajectory[0]   # make sure to be on initial frame
        ref_com     = self._uni_ref.select_atoms(self._reference_group).center_of_mass()
        entire      = self._uni_ref.select_atoms("protein or nucleic or resname CA")
        entire.translate(-ref_com)  
        #new_ref_com = u_ref.select_atoms("resid 495 and name MG").center_of_mass()
        #print ref_com, new_ref_com
        pdbfile = "all-reference.pdb"
        entire.write(pdbfile)
        self._logger.info('Writting %s',pdbfile)
        
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

    def get_centroid(self,in_x,in_y,in_z):
        #Calculates the centroid of given cluster
        sum_x = np.sum(in_x)
        sum_y = np.sum(in_y)
        sum_z = np.sum(in_z)
        length = float(len(in_x))
        c_x = float(sum_x/length)
        c_y = float(sum_y/length)
        c_z = float(sum_z/length)    
        return c_x,c_y,c_z
    
    def generate_clusters(self):
        # based on variace
        #based on iterative variance
        #may be 3D to 2D conversion after that 2D classifications and back conversion
        pass
        
    def plot_contour(self,x,y,z,x1,y1,z1):
        fig = P.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("equal")
        
        ax.scatter(x, y, z,c="#0000FF")
        ax.scatter(x1, y1, z1,c="#FF0000")
        
        xmin=-4;xmax=4
        ymin=-4;ymax=4
        zmin=-4;zmax=5
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_zlim(zmin,zmax)
        
        #draw MG sphere : How to specify the sphere radius?
        u,v = np.mgrid[0:2*np.pi:24j,0:np.pi:10j]
        mg_x=np.cos(u)*np.sin(v)
        mg_y=np.sin(u)*np.sin(v)
        mg_z=np.cos(v)
        ax.plot_wireframe(mg_x,mg_y,mg_z,color='#FFA500')
        
        #print notification total number of points ploted
        N = len(x) #total number of waters
        self._logger.info('Found total %d waters',N)
        P.show()
        '''
        olab="/home/sagar/Desktop/Water-PLOT"
        save="%s.png"%(olab)
        P.savefig(save,dpi=600)
        print "Saved %s"%(save)
        save="%s.ps"%(olab)
        P.savefig(save,dpi=600)
        print "Saved %s"%(save)
        '''
        
        
    def _plot_positions(self):
        
        #NOTE: plot function should read a file and plot iteratively like other plotting scripts
        #NOTE: advised that it be called in main/do_analysis
        fig = P.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("equal")
        myc=['#3BB9FF','#4E9258','#FBB917','#1657FA','#FF00FF','#0000FF','#FF0000','#FF88FF'];
        count=0;
        for residue,coords in self.RESIDUE_DATA.items():

            
            ax.scatter(coords[0], coords[1], coords[2],c=myc[count],label=residue)
            count+=1;
        #ax.scatter(x1, y1, z1,c="#FF0000")
        
        xmin=-4;xmax=4
        ymin=-4;ymax=4
        zmin=-4;zmax=5
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_zlim(zmin,zmax)
            
        #draw center sphere : How to specify the sphere radius?
        u,v = np.mgrid[0:2*np.pi:24j,0:np.pi:10j]
        mg_x=np.cos(u)*np.sin(v)
        mg_y=np.sin(u)*np.sin(v)
        mg_z=np.cos(v)
        ax.plot_wireframe(mg_x,mg_y,mg_z,color='#FF0000')
        #P.legend(bbox_to_anchor=(0., 1.02, 1., .102),prop=mplot.font_manager.FontProperties(size=12),frameon=False,loc=3,ncol=3,mode="expand")
        P.legend(prop=mplot.font_manager.FontProperties(size=16),frameon=False,ncol=3,mode="expand")
        
        #print notification total number of points ploted
        #N = len(x) #total number of waters
        #self._logger.info('Found total %d waters',N)
        #P.show()
        #save="Pos-%s.png"%(self._out_label)
        #P.savefig(save,dpi=600)
    def joint_probabilities(self):
        dist_kk_list=[3,4,3,3,5]
        all_in_one_list=[];
        for key,values in self.RESIDUE_DATA.items():
            #print len(values[3])
            self.all_in_one_list.append(values[3])
            #print all_in_one_list
            #print len(all_in_one_list)
            #print key
            hit=0.0
            for i in range(len(values[3])):
                if (values[3][i]<=3):
                    hit=hit+1
            probability=hit/1000
            self._logger.info('%f',probability)
            #print key, probability
        JP=0.0
        count_list=[];
        #print len(all_in_one_list[1])
        #print len(all_in_one_list)
        for i in range(len(all_in_one_list[1])):
            count=0.0
            for kk in range(len(all_in_one_list)):
                if (all_in_one_list[kk][i]<=dist_kk_list[kk]):
                    count=count+1
                count_list.append(count)
                if count==len(all_in_one_list):
                    JP=JP+1
        #print count_list
        joint_prob=JP/1000
        self._logger.info('%f',joint_prob)
        #print "joint_prob=",joint_prob
    def Geometry_prob_check_external(self):
        #from the first frame
        first_frame=True
        if first_frame==True:
            print ("djb")
    pass

    def individual_probabilities(self):
        for keys,values in self.RESIDUE_DATA.items():
            hit=0
            for i in len(values[3]):
                print (values[3])

    def _get_rotation_matrix(self):
        #calculates rotation matrix
        #uses the PBC fixed system either entire or just the protein
        self._logger.info("Calculating Rotation matrix ....");
        
        self._rotation_matrix=[];
        self._uni_ref.trajectory[0];   # make sure to be on initial frame
        ref= self._uni_ref.select_atoms("name CA");
        ref_com = ref.center_of_mass();
        
        ref0 = ref.positions - ref_com;   
        '''
        
        '''
        for ts in self._uni_ref.trajectory[self._first_frame:self._frame_limit]:

            #calculate time of the frame
            its=ts.frame*ts.dt/1000
            #calculate rotation matrix with respect to CA aroms
            mobile      =   self._uni_ref.select_atoms("name CA")
            mobile_com  =   mobile.center_of_mass()
            mobile0     =   mobile.positions - mobile_com
            Rm, rmsd    =   rotation_matrix(mobile0, ref0); #Rm is the rotation matrix
            self._rotation_matrix.append(Rm)        
        self._logger.info('Length of Rotation matrix .... %8d',len(self._rotation_matrix));

    def _file_checks(self):
        asda_utils.checkfile(self._file_tpr);
        asda_utils.checkfile(self._file_trajectory);
