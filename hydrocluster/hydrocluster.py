#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Thu Apr  5 13:19:39 2018.

@author: lashkov
"""

import argparse
import sys

from .dbcreator.testlist import db_main
from .ui.cli import cli
from .ui.tkgui import TkGui


class Parser(argparse.ArgumentParser):
    """

    """

    def __init__(self):
        super().__init__(prog='hydrocluster.py')
        self.parser_setup()

    def parser_setup(self):
        """

        """
        self.add_argument('-i', '--input', type=str, default='', help='Input file name (pdb)')
        self.add_argument('-emin', '--emin', type=float, default=3.0, help='Minimum EPS value (A)')
        self.add_argument('-emax', '--emax', type=float, default=15.0, help='Maximum EPS value (A)')
        self.add_argument('-es', '--estep', type=float, default=0.1, help='Step of EPS (A)')
        self.add_argument('-smin', '--smin', type=int, default=3, help='Minimum MIN SAMPLES')
        self.add_argument('-smax', '--smax', type=int, default=50, help='Minimum MIN SAMPLES')
        self.add_argument('-g', '--gui', choices=['tkgui', 'cli', 'testlist'], type=str, default='tkgui',
                          help='UI modes')
        self.add_argument('-o', '--output', type=str, default='', help='Output directory name')
        self.add_argument('-c', '--chains', type=str, default=None, help='Selected chains')
        self.add_argument('-rl', '--reslist', type=str, default=None, help='Selected amino acid residues')
        self.add_argument('-pt', '--ptable', choices=['hydropathy', 'menv', 'fuzzyoildrop', 'nanodroplet',
                                                      'aliphatic_core', 'hydrophilic', 'positive', 'negative'],
                          type=str, default='hydropathy', help='Property table for weighting')
        self.add_argument('-pH', '--pH', type=float, default=7.0,
                          help='pH value for calculation of net charges (positive or negative)')
        self.add_argument('-sc', '--score', choices=['si_score', 'calinski', 'dbcv'], type=str, default='calinski',
                          help='Score coefficient')
        self.add_argument('-nf', '--noise_filter', action='store_const', const=True, default=False,
                          help='Activate filter of noise for scoring function (Not recommended!!!')
        self.add_argument('-na', '--noauto', action='store_const', const=True, default=False,
                          help='No automatic mode. --eps and --min_samples options required')
        self.add_argument('-eps', '--eps', type=float, default=3.0, help='EPS value (A)')
        self.add_argument('-min_samples', '--min_samples', type=int, default=3, help='MIN SAMPLES')
        self.add_argument('-ss', '--save_state', action='store_const', const=True, default=False,
                          help='Save states on testlist mode (Required more disk space!)')
        self.add_argument('-pts', '--ptables', type=str, default='all',
                          help='Property tables list for testlist, separator - ","')
        self.add_argument('-scs', '--scores', type=str, default='all',
                          help='Score coefficients list for testlist, separator - ","')


def main():
    """

    """
    if sys.hexversion < 0x30400f0:
        raise RuntimeError("Must use at least python version 3.4")
    if sys.implementation.name != 'cpython':
        raise RuntimeError("Must use CPython implementation")
    parser = Parser()
    namespace = parser.parse_args()
    if namespace.gui == 'tkgui':
        gui = TkGui()
        gui.mainloop()
    elif namespace.gui == 'cli':
        cli(namespace)
    elif namespace.gui == 'testlist':
        db_main(namespace)


if __name__ == '__main__':
    main()
