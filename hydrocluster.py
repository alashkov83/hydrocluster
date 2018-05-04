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
        self.add_argument('-i', '--input', type=str,
                          default='', help='Input file name (pdb)')
        self.add_argument('-s1', '--segment1', type=str,
                          help='Ranges of residues in domain No 1')
        self.add_argument('-s2', '--segment2', type=str,
                          help='Ranges of residues in domain No 2')
        self.add_argument('-g', '--gui', choices=['tkgui', 'cli'], type=str, default='tkgui',
                          help='UI modes')
        self.add_argument('-hd', '--hydrofob', action='store_const', const=True, default=False,
                          help='Only hydrophobic residues.')
        self.add_argument('-n', '--n_cluster', type=int, default=0,
                          help='Number of clusters for clustering. 0 (default) for MeanShift algorithm')
        self.add_argument('-o', '--output', type=str,
                          default='', help='Output file name (dat, xsl, xslx)')
        self.add_argument('-oc', '--ocluster', type=str,
                          default='', help='Output file name for clustering histogram'
                                           ' (eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff)')
        self.add_argument('-of', '--ofigure', type=str,
                          default='', help='Output file name for graphic'
                                           ' (eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff)')


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
