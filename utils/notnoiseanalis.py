#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 30.10.18"""

import argparse
import bz2
import numpy as np
import sys
import pickle

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


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
        self.add_argument('-o', '--output', type=str, default='notnoise.png', help='Output figure (.png)')
        self.add_argument('-min_samples', '--min_samples', type=str, default=None, help='MIN SAMPLES list')



def notnoise_percent(labels):
    labflat = labels.flatten()
    n = len(labflat)
    noise_n = len([x for x in labflat if x == -1])
    return 100-(noise_n * 100 / n)


def loadstate(file: str):
    """

    :param file:
    """
    with bz2.open(file) as f:
        global_state = pickle.load(f)
    htable = global_state['htable']
    states = global_state['states']
    return htable, states


def plot_iso(states, sa, min_samples, htable):
    """

    :param log:
    :param cls:
    :param newdir:
    :param basefile:
    """

    try:
        fig = Figure(figsize=(12, 12))
        ax1 = fig.add_subplot(211)
        ax1.set_title('Percent of not noise points\nPtable: {:s}'.format(htable))
        ax1.set_xlabel('EPS (\u212B)')
        ax1.set_ylabel('%')
        ax1.grid(True)
        ax2 = fig.add_subplot(212)
        ax2.set_title('Derivative percent of not noise points\nPtable: {:s}'.format(htable))
        ax2.set_xlabel('EPS (\u212B)')
        ax2.set_ylabel('d%/dEPS (%/\u212B)')
        ax2.grid(True)
        x_min = []
        x_max = []
        for ms in min_samples:
            x = [state[4] for state in states if state[5] == ms]
            y = [notnoise_percent(state[0]) for state in states if state[5] == ms]
            x, y = list(zip(*sorted((zip(x, y)), key=lambda tup: tup[0])))
            x_max.append(x[y.index(max(y))])
            x_min.append(x[len(x) - list(reversed(y)).index(min(y))])
            dy = np.diff(y)/np.diff(x)
            x0 = (np.diff(x))/2+np.array(x[:-1])
            xmax = x0[dy == np.max(dy)]
            ax1.plot(x, y, label='min_samples: {:d}'.format(ms))
            ax2.plot(x0, dy, label='min_samples: {:d}, EPS: {:.2f}'.format(ms, float(xmax)))
        ax1.legend(loc='best', frameon=False)
        ax1.set_xlim(min(x_min), max(x_max))
        ax2.legend(loc='best', frameon=False)
        ax2.set_xlim(min(x_min), max(x_max))
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError as e:
        print('Error! Failed to plot!!', e)
        sys.exit(-1)


def main():
    parser = Parser()
    namespace = parser.parse_args()
    input_fn = namespace.input
    if not input_fn:
        print("Input file name is not defined!")
        sys.exit(-1)
    output_fn = namespace.output
    htable, states = loadstate(input_fn)
    if namespace.min_samples:
        min_samples = list(map(int, namespace.min_samples.split(',')))
        plot_iso(states, output_fn, min_samples, htable)


if __name__ == '__main__':
    main()

