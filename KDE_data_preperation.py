'''
In this file the data required to run --- code for getting KDE rama tables for rama_prepro score is generated.
'''

########## modules ##########
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
from pyrosetta.rosetta.core.scoring import *
from pyrosetta import *
init()

########## Functions ##########
def get_phi(pose,resi,resj):
    C = pose.residue(resi).atom("C").xyz()
    N = pose.residue(resj).atom("N").xyz()
    CA = pose.residue(resj).atom("CA").xyz()
    C2 = pose.residue(resj).atom("C").xyz()
    phi = rosetta.numeric.dihedral_degrees(C,N,CA,C2)
    return phi

def get_psi(pose,resi,resj):
    N = pose.residue(resi).atom("N").xyz()
    CA = pose.residue(resi).atom("CA").xyz()
    C = pose.residue(resi).atom("C").xyz()
    N2 = pose.residue(resj).atom("N").xyz()
    psi = rosetta.numeric.dihedral_degrees(N,CA,C,N2)
    return psi

def get_phi_psi(pose: Pose) -> list:
    torsion_list = []
    size = pose.size()
    for i in range(1, size+1):
        if i == 1:
            phi = get_phi(pose,size,1)
            psi = get_psi(pose,i,i+1)
        elif i == size:
            phi = get_phi(pose,i-1,i)
            psi = get_psi(pose,size,1)
        else:
            phi = get_phi(pose,i-1,i)
            psi = get_psi(pose,i,i+1)
            
        torsion_list.append(phi)
        torsion_list.append(psi)
        
    return torsion_list
  
  ########## Loading loops ##########
  df = pd.read_csv('loops-dataset.csv')
