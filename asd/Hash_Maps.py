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
def aa2res(aa):
    aadict = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "CA": "CA",
        "ZN": "MG",
        "ZN": "ZN",
        "FE": "FE"                 
        }
    return aadict[aa]
def res_type(residue):
    dict={"ALA":"PRO",
          "ARG":"PRO",
          "ASN":"PRO",
          "ASP":"PRO",
          "CYS":"PRO",
          "GLN":"PRO",
          "GLU":"PRO",
          "GLY":"PRO",
          "HIS":"PRO",
          "HIP":"PRO",
          "HIE":"PRO",
          "ILE":"PRO",
          "LEU" :   "PRO",
          "LYS" :   "PRO",
          "MET" :   "PRO",
          "PHE" :   "PRO",
          "PRO" :   "PRO",
          "SER" :   "PRO",
          "THR" :   "PRO",
          "TYR" :   "PRO",
          "TRP" :   "PRO",
          "VAL" :   "PRO",
          "DC"  :   "DNA",
          "DG"  :   "DNA",
          "DT"  :   "DNA",
          "DA"  :   "DNA",
          "A"   :   "DNA",
          "G"   :   "DNA",
          "T"   :   "DNA",
          "C"   :   "DNA",
          "U"   :   "RNA",
          "RC"  :   "RNA",
          "RG"  :   "RNA",
          "RU"  :   "RNA",
          "RA"  :   "RNA",
          "CA"  :   "ION",
          "ZN"  :   "ION",
          "FE"  :   "ION",
          "CU"  :   "ION",
          "SOL"  :  "WATER",
                    
          }
    if(residue in dict):
        return dict[residue]
    else:
        return "LIG"

def res_selection(residue):
    dict={"ALA" :   "(name CB)",
          "ARG" :   "(name NH1 or name NH2)",
          "ASN" :   "(name OD1 or name ND2)",
          "ASP" :   "(name OD1 or name OD2)",
          "CYS" :   "(name SG)",
          "GLN" :   "(name OE1 or name NE2)",
          "GLU" :   "(name OE1 or name OE2)",
          "GLY" :   "(name CA)",
          "HIS" :   "(name ND1 or name NE2)",
          "HIP" :   "(name ND1 or name NE2)",
          "HIE" :   "(name ND1 or name NE2)",
          "HID" :   "(name ND1 or name NE2)",
          "ILE" :   "(name CG2)",
          "LEU" :   "(name CD1 or name CD2)",
          "LYS" :   "(name NZ)",
          "MET" :   "(name CE)",
          "PHE" :   "(name CG or name CD1 or name CE1 or name CZ or name CE2 or name CD2)",
          "PRO" :   "(name CB or name CG or name CD)",
          "SER" :   "(name OG)",
          "THR" :   "(name OG1 or name CG2)",
          "TYR" :   "(name OH)",
          "TRP" :   "PRO",
          "VAL" :   "(name CG1 or name CG2)",
          "CA"  :   "(name CA)",
          "ZN"  :   "(name ZN)",
          "FE"  :   "(name FE)",
          "CU"  :   "(name CU)",
          "SOL" :   "(name OW)",
          }
    if(residue in dict):
        return dict[residue]
    else:
        return "UNK"
def res_selection_atoms(residue):
    dict={"ALA" :   ["CB"],
          "ARG" :   ["NH1", "NH2"],
          "ASN" :   ['OD1', 'ND2'],
          "ASP" :   ["OD1","OD2"],
          "CYS" :   ["SG"],
          "GLN" :   ["OE1","NE2"],
          "GLU" :   ["OE1","OE2"],
          "GLY" :   ["CA"],
          "HIS" :   ["ND1","NE2"],
          "HIP" :   ["ND1","NE2"],
          "HIE" :   ["ND1","NE2"],
          "HID" :   ["ND1","NE2"],
          "ILE" :   ["CG2"],
          "LEU" :   ["CD1","CD2"],
          "LYS" :   ["NZ"],
          "MET" :   ["CE"],
          "PHE" :   ["CG","CD1","CE1","CZ","CE2","CD2"],
          "PRO" :   ["CB","CG","CD"],
          "SER" :   ["OG"],
          "THR" :   ["OG1","CG2"],
          "TYR" :   ["OH"],
          "TRP" :   ["PRO"],
          "VAL" :   ["CG1","CG2"],
          "CA"  :   ["CA"],
          "ZN"  :   ["ZN"],
          "FE"  :   ["FE"],
          "SOL" :   ["OW"],
                    
          }
    if(residue in dict):
        return dict[residue]
    else:
        return "UNK"
def res2aa_upgrade(residue):
    dict={"ALA":"A",
          "ARG":"R",
          "ASN":"N",
          "ASP":"D",
          "CYS":"C",
          "GLN":"Q",
          "GLU":"E",
          "GLY":"G",
          "HIS":"H",
          "HIP":"H",
          "HIE":"H",
          "ILE":"I",
          "LEU":"L",
          "LYS":"K",
          "MET":"M",
          "PHE":"F",
          "PRO":"P",
          "SER":"S",
          "THR":"T",
          "TYR":"Y",
          "TRP":"W",
          "VAL":"V",
          "DC":"NA",
          "DG":"NA",
          "DT":"NA",
          "DA":"NA",
          "C":"NA",
          "G":"NA",
          "T":"NA",
          "A":"NA",
          "U":"NA",
          "RC":"NA",
          "RG":"NA",
          "RU":"NA",
          "RA":"NA"          
          }
    if(residue in dict):
        return dict[residue]
    else:
        return "UNK"
'''
          "DC":"DC",
          "DG":"DG",
          "DT":"DT",
          "DA":"DA",
          "C":"C",
          "G":"G",
          "T":"T",
          "A":"A",
          "U":"U",
          "RC":"RC",
          "RG":"RG",
          "RU":"RU",
          "RA":"RA",
          "KCX":"K",
          "MSE":"S"

'''

def atom_mass(atom):
    atommass = {
        "H":1.008,
        "C":12.01,
        "N":14.0067,
        "O":16.001,
        "S":32.065,
        "P":30.97,
        "Na":22.989,
        "Mg":24.3050,
        "Cl":35.453,
        "Ca":40.078
    }
    return atommass[atom]
def experiments(exp):
    expt={
        "X-RAY-DIFFRACTION":"X-Ray",
        "SOLUTION-NMR":"NMR",
        "ELECTRON-MICROSCOPY":"EM",
        "ELECTRON-CRYSTALLOGRAPHY":"EC",
        "SOLID-STATE-NMR":"SSNMR"
        }
    if(exp in expt):
        return expt[exp]
    else:
        print(exp)
        exit(1)