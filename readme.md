# Averaging Scripts

## Requirements for installation:

1. Working and up to date install of python, ProDy, and networkx
    - pip install prody
    - pip install networkx
2. For do_minimization.py, a working installation of gromacs 4.5.

## The scripts

All trajectories should be in multi-frame PDB format.


cartesian_average.py trajectory.pdb average.pdb

Iteratively superimpose all frames of the trajectory onto the average and then
output the cartesian average of all heavy atoms.


find_centroid.py trajectory.pdb centroid.pdb

Compute the average as for cartesian_average.py and then output the frame from
the trajectory that is closest to the average.


network_average.py --cutoff 1.0 trajectory.pdb average.pdb

Compute the NetAvg for trajectory. A cutoff of 1.0 works best on our benchmark
set.


do_minimization.py start.pdb target.pdb final.pdb

Apply position restraints to move start closer to target, while maintaining a
physical structure. A working version of gromacs must be in your path for this
step. Uses implicit solvent and amber99sb-ildn forcefield, but the forcefield
shouldn't matter much as the atom positions are largely dicated by the average.


## Potential Problems

Other than issues with installtion, the two most likely causes of probelms is
atom naming. These tools work with my benchmark tests, but I suspect there will
be problems when using other data sets due to differences in heavy atom naming
conventions.

