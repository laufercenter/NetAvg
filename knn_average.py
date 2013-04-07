#!/usr/bin/env python
# encoding: utf-8


import argparse
import prody
import numpy
from numpy import matlib
import networkx as nx
import operator


def parse_args():
    parser = argparse.ArgumentParser(description='Perform cartesian average of trajectory.')
    parser.add_argument('--knn', type=int, default=32, help='number of nearest neighbors to consider')
    parser.add_argument('trajectory', metavar='TRAJECTORY', help='input trajectory in PDB format')
    parser.add_argument('out_file', metavar='OUTPUT_PDB_FILE', help='output PDB file')

    args = parser.parse_args()
    return args


def load_pdb(input_file):
    return prody.parsePDB(input_file)


def get_dist_matrix(points):
    numPoints = len(points)
    distMat = numpy.sqrt(numpy.sum((matlib.repmat(points, numPoints, 1) - matlib.repeat(points, numPoints, axis=0))**2,
                         axis=1))
    return distMat.reshape((numPoints, numPoints))


def get_network_average(coords, knn):
    # create all nodes
    print '\tCreating graph'
    g = nx.DiGraph()
    n_frames = coords.shape[0]
    g.add_nodes_from(range(n_frames))

    # get distances
    print '\tCalculating distance matrix'
    distances = get_dist_matrix(coords)

    # add edges
    print '\tAdding edges'
    for i in range(n_frames):
        # find knn + 1 nearest neighbors
        neighbors = numpy.argsort(distances[i, :])[:knn+1]
        for j in neighbors:
            if not i == j:
                g.add_edge(i, j)

    # turn the graph into a directed graph
    g = g.to_undirected(reciprocal=True)

    # find the largest connected component
    connected_components = nx.connected_component_subgraphs(g)
    largest_connected_component = connected_components[0]
    print '\t{} connected components.'.format( len(connected_components) )
    print '\tLargest connected component contains {} nodes.'.format( len(largest_connected_component) )

    # find most central point
    print '\tFinding most central point'
    centrality = nx.eigenvector_centrality_numpy(largest_connected_component)
    index_of_max = max( centrality.iteritems(), key=operator.itemgetter(1) )[0]

    # calculate the average of central point and all it's neighbors
    neighbors = largest_connected_component.neighbors(index_of_max)
    neighbors.append(index_of_max)
    average = numpy.mean( coords[neighbors,:], axis=0 )

    # return it
    return average


def calc_average(trajectory, knn):
    output = prody.AtomGroup('Cartesian average coordinates')
    output_coords = trajectory.getCoords()
    output.setCoords( trajectory.getCoords() )
    output.setNames( trajectory.getNames() )
    output.setResnums( trajectory.getResnums() )
    output.setResnames( trajectory.getResnames() )

    ensemble = prody.PDBEnsemble(trajectory)
    ensemble.iterpose()

    print 'Using knn of {}'.format(knn)

    input_coords = ensemble.getCoordsets()

    n_atoms = output_coords.shape[0]
    for i in range(n_atoms):
        print 'Computing residue {} of {}:'.format(i+1, n_atoms)
        average = get_network_average( input_coords[:,i,:], knn )
        output_coords[i,:] = average

    output.setCoords(output_coords)
    return output


def main():
    args = parse_args()

    trajectory = load_pdb(args.trajectory)
    average = calc_average(trajectory, args.knn)

    prody.writePDB(args.out_file, average)


if __name__ == '__main__':
    main()
