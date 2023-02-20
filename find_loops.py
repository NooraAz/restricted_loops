## ==== importing modules =====
import os
import itertools as it
import math
import pandas as pd
from pyrosetta import *
init()

## ==== searching for loops =====

# Getting a list of all PDB ids
df = pd.read_csv('resolu.idx',header=None)

# searching for loops
for name in df[0]:
    print(name, end =' ')
    # getting pose
    pose = toolbox.rcsb.pose_from_rcsb(name)
    # getting secondary structure of pose
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose)
    secondary_structure = pose.secstruct()
    # finding all continuous loops
    targets = []
    c = 1
    for k,g in it.groupby(secondary_structure):
        le = len(list(g))
        if k == 'L' and le <= 10 and le >= 4 and c != 1 :
            print(k,c,le,end=' ')
    # getting the distance between CA of the beginning and end of each loop
            a = pose.residue(c).xyz(pose.residue(c).atom_index('CA'))
            x = c+le-1
            b = pose.residue(x).xyz(pose.residue(x).atom_index('CA'))
            d = math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)
    # if distance < 3.9 A then add ['L',starting_residue_number,length_of_loop] to targets
            if d <= 3.9:
                targets.append([k,c,le])
                print('yes')
            else:
                print('no')
        c += le
    # if there were any restricted loops found, add them as a value to "lib" dictionary with pdb_name key
    if len(targets) > 0:
        lib[name] = targets

# === dumping loops ====
for pdb_name,loop in lib.items():
    pose = toolbox.rcsb.pose_from_rcsb(pdb_name)
    for a,b,c in loop:
        cloned_pose = pose.clone()
        # deleting pre-loop residues
        drm = rosetta.protocols.grafting.simple_movers.DeleteRegionMover()
        drm.region('1',str(b-1))
        drm.apply(cloned_pose)
        size = str(cloned_pose.total_residue())
        # deleting pro-loop reidues
        drm.region(str(c+1),size)
        drm.apply(cloned_pose)
        # dump the remaining
        cloned_pose.dump_pdb(f'{pdb_name}-{b}-{c}.pdb')
