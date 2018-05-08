#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Thu Apr  5 13:19:39 2018.

@author: lashkov
"""

import argparse


class Parser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__()
        self.parser_setup()

    def parser_setup(self):
        self.add_argument('-i', '--input', type=str, default='', help='Input file name (pdb)')
        self.add_argument('-emin', '--emin', type=float, default=3.0, help='Minimum EPS value (A)')
        self.add_argument('-emax', '--emax', type=float, default=15.0, help='Maximum EPS value (A)')
        self.add_argument('-es', '--estep', type=float, default=0.1, help='Step of EPS (A)')
        self.add_argument('-smin', '--smin', type=int, default=2, help='Minimum MIN SAMPLES')
        self.add_argument('-smax', '--smax', type=int, default=50, help='Minimum MIN SAMPLES')
        self.add_argument('-g', '--gui', choices=['tkgui', 'cli'], type=str, default='tkgui', help='UI modes')
        self.add_argument('-o', '--output', type=str, default='', help='Output directory name')
        self.add_argument('-sc', '--score', choices=['si_score', 'calinski'], type=str, default='calinski',
                          help='Score coefficient')
        self.add_argument('-ht', '--htable', choices=['hydropathy', 'nanodroplet'], type=str, default='hydropathy',
                          help='Hydrophobity table for weighting')
        self.add_argument('-na', '--noauto', action='store_const', const=True, default=False,
                          help='No automatic mode. --eps and --min_samples options are required')
        self.add_argument('-eps', '--eps', type=float, default=0, help='EPS value (A)')
        self.add_argument('-min_samples', '--min_samples', type=int, default=0, help='MIN SAMPLES')


if __name__ == '__main__':
    parser = Parser()
    namespace = parser.parse_args()
    if namespace.gui == 'tkgui':
        from hydrocluster.tkgui import TkGui

        gui = TkGui(namespace)
        gui.mainloop()
    elif namespace.gui == 'cli':
        from hydrocluster.cli import Cli

        cli = Cli(namespace)