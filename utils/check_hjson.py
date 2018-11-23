#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 22.11.18"""

import argparse
import os.path
import sys

import hjson


class Parser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__(prog='check_hjson.py')
        self.parser_setup()

    def parser_setup(self):
        self.add_argument('-i', '--input', type=str, default='default.hjson', help='Input hjson filename')


parser = Parser()
namespace = parser.parse_args()
input_fn = namespace.input
if not input_fn:
    print("Input file name is not defined!")
    sys.exit(-1)
if not os.path.exists(input_fn):
    print("Input file was not found!")
    sys.exit()
with open('default.hjson') as f:
    options_dict = hjson.load(f)
for key in options_dict:
    print("key: ", key, ", value: ", str(options_dict[key]), ", type:", str(type(options_dict[key])))
