#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 24.10.18"""
import argparse
import bz2
import pickle
import sys
from scipy.spatial import ConvexHull


class Parser(argparse.ArgumentParser):
    """

    """

    def __init__(self):
        super().__init__(prog='hullanalis.py')
        self.parser_setup()

    def parser_setup(self):
        """

        """
        self.add_argument('-i', '--input', type=str, default='', help='Input file name (.dat)')

def loadstate(file: str):
    """

    :param file:
    """
    with bz2.open(file) as f:
        global_state = pickle.load(f)
    X = global_state['X']
    # pdist = global_state['pdist']
    labels = global_state['labels']
    # sparse_n = global_state['sparse_n']
    # noise_filter = global_state['noise_filter']
    core_samples_mask = global_state['core_samples_mask']
    n_clusters = global_state['n_clusters']
    # s_array = global_state['s_array']
    # htable = global_state['htable']
    # parse_results = global_state['parse_results']
    # auto_params = global_state['auto_params']
    # score = global_state['score']
    # eps = global_state['eps']
    # min_samples = global_state['min_samples']
    # metric = global_state['auto_params'][-1]
    # weight_array = global_state['weight_array']
    # aa_list = global_state['aa_list']
    # states = global_state['states']
    # figs = global_state['figs']
    return X, labels, core_samples_mask, n_clusters


def convexhull(points, ncl=1):
    """

    :param points:
    :return:
    """

    if ncl == 1:
        ch = ConvexHull(points)
        v = ch.volume
        n = len(points)
        c = n/v
        return v, c, n
    elif ncl > 1:
        from sklearn.cluster import k_means
        vs = []
        while ncl > 1:
            try:
                centroid, labels, inetrtia = k_means(points, ncl)
                for i in set(labels):
                    ch = ConvexHull(points[labels == i])
                    v = ch.volume
                    vs.append(v)
            except Exception:
                ncl -= 1
                vs.clear()
            else:
                break
        v = sum(vs)
        n = len(points)
        c = n/v
        return v, c, n

def main():
    parser = Parser()
    namespace = parser.parse_args()
    input_fn = namespace.input
    if not input_fn:
        print("Input file name is not defined!")
        sys.exit(-1)
    X, labels, core_samples_mask, n_clusters = loadstate(input_fn)
    Vcores = []
    Ncores = []
    for i in range(n_clusters):
        filterXYZ = X[labels == i]
        filterXYZcore = X[(labels == i) & core_samples_mask]
        V, C, _ = convexhull(filterXYZ)
        Vcore, Ccore, Ncore = convexhull(filterXYZcore)
        Vcores.append(Vcore)
        Ncores.append(Ncore)
        print("For cluster No.{:d}: V= {:.2f} \u212B\u00B3, C= {:.4f} \u212B\u207B\u00B3, "
              "Vcore= {:.2f} \u212B\u00B3, Ccore = {:.3f} \u212B\u207B\u00B3".format(i+1, V, C, Vcore, Ccore))
    V, C, _ = convexhull(X)
    print("For all hydrophobic residues: V= {:.2f} \u212B\u00B3, C= {:.4f} \u212B\u207B\u00B3".format(V, C))
    print("Mean for core hydrophobic C= {:.4f} \u212B\u207B\u00B3".format(sum(Ncores)/sum(Vcores)))


if __name__ == '__main__':
    main()