#!/usr/bin/env python
# encoding: utf-8


import argparse
import prody
import numpy
from matplotlib import pyplot


def get_traj_rmsds(reference, trajectory):
    ref_backbone = reference.select('backbone or name OC2')
    traj_backbone = trajectory.select('backbone or name OC2')

    ensemble = prody.Ensemble('trajectory ensemble')
    ensemble.setCoords( ref_backbone.getCoords() )
    ensemble.addCoordset( traj_backbone.getCoordsets() )

    ensemble.superpose()
    return ensemble.getRMSDs()


def get_single_rmsd(reference, model):
    ref_backbone = reference.select('backbone or name OC2')
    mod_backbone = model.select('backbone or name OC2')

    prody.superpose(mod_backbone, ref_backbone)
    return prody.calcRMSD(mod_backbone, ref_backbone)


def parse_args():
    parser = argparse.ArgumentParser(description='Perform cartesian average of trajectory.')
    parser.add_argument('reference')
    parser.add_argument('trajectory')
    parser.add_argument('cart')
    parser.add_argument('cart_min')
    parser.add_argument('net')
    parser.add_argument('net_min')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    reference = prody.parsePDB(args.reference)
    trajectory = prody.parsePDB(args.trajectory)
    cart = prody.parsePDB(args.cart)
    cart_min = prody.parsePDB(args.cart_min)
    net = prody.parsePDB(args.net)
    net_min = prody.parsePDB(args.net_min)

    cart_rmsd = get_single_rmsd(reference, cart)
    cart_min_rmsd = get_single_rmsd(reference, cart_min)
    net_rmsd = get_single_rmsd(reference, net)
    net_min_rmsd = get_single_rmsd(reference, net_min)
    traj_rmsds = get_traj_rmsds(reference, trajectory)

    pyplot.hist(traj_rmsds, range=(0, 5), bins=26, histtype='stepfilled',
            color='lightgrey', normed=True, edgecolor='none')
    pyplot.axvline(cart_rmsd, color='blue', linewidth=2, linestyle='--')
    pyplot.axvline(cart_min_rmsd, color='blue', linewidth=2, linestyle='-')
    pyplot.axvline(net_rmsd, color='orange', linewidth=2, linestyle='--')
    pyplot.axvline(net_min_rmsd, color='orange', linewidth=2, linestyle='-')

    pyplot.show()


if __name__ == '__main__':
    main()

