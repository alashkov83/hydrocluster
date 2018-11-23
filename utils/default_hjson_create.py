#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 22.11.18"""

import argparse
from collections import OrderedDict

import hjson


class Parser(argparse.ArgumentParser):
    """

    """

    def __init__(self):
        super().__init__(prog='default_hjson_create.py')
        self.parser_setup()

    def parser_setup(self):
        self.add_argument('-o', '--output', type=str, default='default.hjson', help='Output hjson filename')


options_dict = OrderedDict()
options_dict['Input file name'] = 'input.txt'
options_dict['Project name'] = 'project'
options_dict['Property tables list'] = ('hydropathy',
                                        'nanodroplet',
                                        'menv',
                                        'rekkergroup',
                                        'fuzzyoildrop',
                                        'aliphatic_core',
                                        'hydrophilic',
                                        'positive',
                                        'negative',
                                        'ngroup',
                                        'pgroup')
options_dict['pH'] = 7.0
options_dict['Project name'] = 'project'
options_dict['Minimum eps value (A)'] = 3.0
options_dict['Maximum eps value (A)'] = 15.0
options_dict['Step of eps value (A)'] = 0.1
options_dict['Minimum min_samples'] = 3
options_dict['Maximum min_samples'] = 50
options_dict['Scoring coefficients list'] = ('si_score',
                                             'calinski',
                                             's_dbw')
options_dict['Save states'] = False
parser = Parser()
namespace = parser.parse_args()
output_fn = namespace.output
with open(output_fn, 'w') as f:
    hjson.dump(options_dict, f)
