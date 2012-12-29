#!/usr/bin/env python
# encoding: utf-8


import argparse
import prody
import os
import shutil
import subprocess
import numpy


mdp_string = '''
define = -DPOSRES
integrator = {integrator}
nsteps = 1000
emtol = 1
nstlist = 1
coulombtype = Cut-off
vdwtype = Cut-off
ns_type = simple
rlist = 1.8
rcoulomb = 1.8
rvdw = 1.8
pbc = xyz
implicit_solvent = GBSA
gb_algorithm = OBC
sa_algorithm = ACE-approximation
rgbradii = 1.8
;nstxout = 1
'''


def parse_args():
    parser = argparse.ArgumentParser(description='Generate trajectory with gaussian flucutations.')
    parser.add_argument('pdb_file', metavar='INPUT_PDB_FILE', help='path to input pdb file')
    parser.add_argument('trajectory', metavar='TRAJECTORY', help='path to input trajectory')
    parser.add_argument('out_file', metavar='OUTPUT_PDB_FILE', help='path to input pdb file')
    parser.add_argument('--skip-first', action='store_true', default=False)
    args = parser.parse_args()

    return args.pdb_file, args.trajectory, args.out_file, args.skip_first


def get_closest_frame(trajectory, average_structure, skip_first):
    output = prody.AtomGroup('Cartesian average coordinates')
    output.setCoords( trajectory.getCoords() )
    output.setNames( trajectory.getNames() )
    output.setResnums( trajectory.getResnums() )
    output.setResnames( trajectory.getResnames() )

    ensemble = prody.PDBEnsemble(trajectory)
    if skip_first:
        ensemble.delCoordset(0)
    ensemble.setCoords( average_structure.getCoords() )
    ensemble.superpose()
    rmsds = ensemble.getRMSDs()
    min_index = numpy.argmin(rmsds)

    output.setCoords( ensemble.getCoordsets(min_index) )
    return output



def load_pdb(in_file):
    protein = prody.parsePDB(in_file)
    return protein


def create_no_h_file():
    # make the index file
    cmd = 'make_ndx -f min_round_2.gro -o no_h.ndx'
    p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
    p1.communicate('q\n')

    # run editconf
    edit_cmd = 'editconf -f min_round_2.gro -o no_h.gro -n no_h.ndx'
    p2 = subprocess.Popen(edit_cmd, shell=True, stdin=subprocess.PIPE)
    p2.communicate('2\n')


def re_order():
    # create a new index file
    lines = open('index.ndx').read().splitlines()
    header = lines[0]
    indices = []
    for line in lines[1:]:
        cols = line.split()
        for col in cols:
            indices.append( int(col) )
    resorted = [ indices.index(val)+1 for val in range( 1, max(indices)+1 ) ]
    with open('resort.ndx', 'w') as out:
        print >>out, header
        for val in resorted:
            print >>out, val

    # resort
    edit_cmd = 'editconf -f no_h.gro -o min.pdb -n resort.ndx'
    subprocess.check_call(edit_cmd, shell=True)


def  main():
    in_file, trajectory, out_file, skip_first = parse_args()
    average = load_pdb(in_file)
    trajectory = load_pdb(trajectory)
    starting_model = get_closest_frame(trajectory, average, skip_first=skip_first)
    prody.writePDB(out_file, starting_model)


if __name__ == '__main__':
    main()

