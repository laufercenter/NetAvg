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
    args = parser.parse_args()

    return (args.pdb_file, args.trajectory, args.out_file)


def get_closest_frame(trajectory, average_structure):
    output = prody.AtomGroup('Cartesian average coordinates')
    output.setCoords( trajectory.getCoords() )
    output.setNames( trajectory.getNames() )
    output.setResnums( trajectory.getResnums() )
    output.setResnames( trajectory.getResnames() )

    ensemble = prody.PDBEnsemble(trajectory)
    ensemble.setCoords(average_structure)
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
    #edit_cmd = 'editconf -f min_round_2.gro -o no_h.gro -n no_h.ndx'
    #p2 = subprocess.Popen(edit_cmd, shell=True, stdin=subprocess.PIPE)
    #p2.communicate('2\n')


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
    #edit_cmd = 'editconf -f no_h.gro -o min.pdb -n resort.ndx'
    edit_cmd = 'editconf -f min_round_2.gro -o min.pdb -n resort.ndx'
    subprocess.check_call(edit_cmd, shell=True)


def run_minimization(average, start):
    # create temp dir
    os.mkdir('Temp')
    os.chdir('Temp')

    # write the average file
    prody.writePDB('average.pdb', average)
    pdb_cmd = 'pdb2gmx -f average.pdb -ff amber99sb-ildn -water none -n index.ndx -posrefc 1000 -o ref.gro -his'
    p = subprocess.Popen(pdb_cmd, shell=True, stdin=subprocess.PIPE)
    p.communicate('0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n')
    # put it in a bigger box
    box_cmd = 'editconf -f ref.gro -o ref_box.gro -c -box 999 999 999'
    subprocess.check_call(box_cmd, shell=True)

    # write pdb file
    prody.writePDB('start.pdb', start)

    # pdb2gmx
    pdb_cmd = 'pdb2gmx -f start.pdb -ff amber99sb-ildn -water none -n index.ndx -posrefc 1000 -his'
    p = subprocess.Popen(pdb_cmd, shell=True, stdin=subprocess.PIPE)
    p.communicate('0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n')

    # put it in a bigger box
    box_cmd = 'editconf -f conf.gro -o box.gro -c -box 999 999 999'
    subprocess.check_call(box_cmd, shell=True)

    #
    # Round 1
    #

    # write mdp file
    with open('min_round_1.mdp', 'w') as min_file:
        min_file.write( mdp_string.format(integrator='steep') )

    # run grompp
    grompp_cmd = 'grompp -f min_round_1.mdp -c box.gro -p topol.top -o min_round_1 -r ref_box.gro'
    subprocess.check_call(grompp_cmd, shell=True)

    # run mdrun
    md_cmd = 'mdrun -deffnm min_round_1 -v -nt 1'
    subprocess.check_call(md_cmd, shell=True)

    #
    # Round 2
    #

    # write mdp file
    with open('min_round_2.mdp', 'w') as min_file:
        min_file.write( mdp_string.format(integrator='l-bfgs') )

    # run grompp
    grompp_cmd = 'grompp -f min_round_2.mdp -c min_round_1.gro -p topol.top -o min_round_2 -maxwarn 1 -r ref_box.gro'
    subprocess.check_call(grompp_cmd, shell=True)

    # run mdrun
    md_cmd = 'mdrun -deffnm min_round_2 -v -nt 1'
    subprocess.check_call(md_cmd, shell=True)

    #
    # gather results
    #
    create_no_h_file()
    re_order()

    # load the pdb
    #protein = prody.parsePDB('min.pdb').select('not hydrogen')
    protein = prody.parsePDB('min.pdb').select('all')

    # clean up
    os.chdir('..')
    shutil.rmtree('Temp')

    return protein


def  main():
    in_file, trajectory, out_file = parse_args()
    average = load_pdb(in_file)
    trajectory = load_pdb(trajectory)
    starting_model = get_closest_frame(trajectory, average)
    minimized_protein = run_minimization(average, starting_model)
    prody.writePDB(out_file, minimized_protein)


if __name__ == '__main__':
    main()

