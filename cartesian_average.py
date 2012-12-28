#!/usr/bin/env python
# encoding: utf-8


import argparse
import prody
import numpy


def parse_args():
    parser = argparse.ArgumentParser(description='Perform cartesian average of trajectory.')
    parser.add_argument('trajectory', metavar='TRAJECTORY', help='input trajectory in PDB format')
    parser.add_argument('out_file', metavar='OUTPUT_PDB_FILE', help='output PDB file')

    args = parser.parse_args()
    return args


def load_pdb(input_file):
    return prody.parsePDB(input_file)


def calc_average(trajectory):
    output = prody.AtomGroup('Cartesian average coordinates')
    output.setCoords( trajectory.getCoords() )
    output.setNames( trajectory.getNames() )
    output.setResnums( trajectory.getResnums() )
    output.setResnames( trajectory.getResnames() )

    ensemble = prody.PDBEnsemble(trajectory)
    ensemble.iterpose()

    coords = ensemble.getCoordsets()
    average_coords = numpy.mean(coords, axis=0)
    output.setCoords(average_coords)
    return output


def main():
    args = parse_args()

    trajectory = load_pdb(args.trajectory)
    average = calc_average(trajectory)

    prody.writePDB(args.out_file, average)


if __name__ == '__main__':
    main()
