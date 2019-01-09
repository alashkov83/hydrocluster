#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 16.10.18"""
import argparse
import bz2
import pickle
import sys

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

metrics_name = {'calinski': 'Calinski-Harabaz score', 'si_score': 'Silhouette score', 'dbcv': 'DBCV score'}

epsilon = 0.00000001


class Parser(argparse.ArgumentParser):
    """

    """

    def __init__(self):
        super().__init__(prog='hyperanalis.py')
        self.parser_setup()

    def parser_setup(self):
        """

        """
        self.add_argument('-i', '--input', type=str, default='', help='Input file name (.dat)')
        self.add_argument('-o', '--output', type=str, default='hiper.png', help='Output figure (.png)')
        self.add_argument('-eps', '--eps', type=float, default=None, help='EPS value (A)')
        self.add_argument('-min_samples', '--min_samples', type=int, default=None, help='MIN SAMPLES')


def loadstate(file: str):
    """

    :param file:
    """
    with bz2.open(file) as f:
        global_state = pickle.load(f)
    # X = global_state['X']
    # pdist = global_state['pdist']
    # labels = global_state['labels']
    # sparse_n = global_state['sparse_n']
    # noise_filter = global_state['noise_filter']
    # core_samples_mask = global_state['core_samples_mask']
    # n_clusters = global_state['n_clusters']
    # s_array = global_state['s_array']
    htable = global_state['htable']
    # parse_results = global_state['parse_results']
    # auto_params = global_state['auto_params']
    # score = global_state['score']
    # eps = global_state['eps']
    # min_samples = global_state['min_samples']
    metric = global_state['auto_params'][-1]
    # weight_array = global_state['weight_array']
    # aa_list = global_state['aa_list']
    states = global_state['states']
    # figs = global_state['figs']
    return htable, metric, states


def colormap(x, y, htable, metric, xparametr, sa: str, const_str):
    """

    :param log:
    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """

    try:
        fig = Figure(figsize=(12, 6))
        ax1 = fig.add_subplot(111)
        ax1.set_title(metrics_name[metric] + ' vs ' + xparametr + '\nhtable: ' + htable + ", " + const_str)
        ax1.set_xlabel(xparametr)
        ax1.set_ylabel(metrics_name[metric])
        ax1.grid(True)
        if xparametr == 'EPS (\u212B)':
            ax1.plot(x, y)
        elif xparametr == 'MIN_SAMPLES':
            ax1.bar(x, y)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError as e:
        print('Error! Failed to plot!!', e)
        sys.exit(-1)


def calculate_xy(states, param, xparm):
    if xparm == 'eps':
        x = [state[5] for state in states if abs(state[4] - param) <= epsilon]
        y = [state[3] for state in states if abs(state[4] - param) <= epsilon]
    elif xparm == 'min_samples':
        x = [state[4] for state in states if state[5] == param]
        y = [state[3] for state in states if state[5] == param]
    x, y = list(zip(*sorted((zip(x, y)), key=lambda tup: tup[0])))
    return x, y


def main():
    parser = Parser()
    namespace = parser.parse_args()
    input_fn = namespace.input
    if not input_fn:
        print("Input file name is not defined!")
        sys.exit(-1)
    output_fn = namespace.output
    htable, metric, states = loadstate(input_fn)
    eps = namespace.eps
    min_samples = namespace.min_samples
    if (eps or min_samples) and not (eps and min_samples):
        if eps:
            x, y = calculate_xy(states, eps, 'eps')
            xparametr = 'MIN_SAMPLES'
            const_str = 'EPS: {:.2f} \u212B'.format(eps)
        else:
            x, y = calculate_xy(states, min_samples, 'min_samples')
            xparametr = 'EPS (\u212B)'
            const_str = 'MIN_SAMPLES: {:d}'.format(min_samples)

    else:
        print('Only one parameter (eps or min_samples) will bee given!')
        sys.exit(-1)
    colormap(x, y, htable, metric, xparametr, output_fn, const_str)


if __name__ == '__main__':
    main()
