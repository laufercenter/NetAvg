#!/usr/bin/env python

import argparse
import prody
import numpy


prody.confProDy(verbosity='none')


def parse_args():
    parser = argparse.ArgumentParser(description='Find the RMSD, taking the closest position for each atom')
    parser.add_argument('native', help='native structure in pdb format')
    parser.add_argument('trajectory', help='input trajectory in pdb format')
    parser.add_argument('--skip-frames', type=int, default=0, help='number of frames to skip from start')
    args = parser.parse_args()
    return args


def find_close(native_name, traj_name, skip_frames):
    native = prody.parsePDB(native_name)
    traj = prody.parsePDB(traj_name)

    ensemble = prody.Ensemble('ensemble')
    ensemble.setCoords(native.getCoords())
    ensemble.addCoordset(traj.getCoordsets()[skip_frames:, ...])  # skip the first 10 frames
    ensemble.superpose()

    native_coords = native.getCoords()
    ensemble_coords = ensemble.getCoordsets()

    diff2 = (ensemble_coords - native_coords) ** 2
    diff2 = numpy.sum(diff2, axis=2)
    min_dev = numpy.min(diff2, axis=0)
    return numpy.sqrt(numpy.sum(min_dev) / float(min_dev.shape[0]))


def main():
    args = parse_args()
    result = find_close(args.native, args.trajectory, args.skip_frames)
    print result


if __name__ == '__main__':
    main()
