# LDAT (v0.1b)
Python: Ligand densities Analsyis tool for Molecular Dynamics trajectories.

#### **Requirements**:
- Python3.6, numpy 1.15
- mdanalysis

At present: works only for water molecules

Searches for water centered around a residue 'x' with radius 'y' nm and reports densities at each point in the grid.

#### **Usage**:

Download and place the code in: /home/`whoami`/ldat

LDATPATH="/home/`whoami`/ldat"
export PATH=$PATH:$LDATPATH/bin
export PYTHONPATH=$PYTHONPATH:$LDATPATH

in terminal: 
gmx_ldat.py -h
