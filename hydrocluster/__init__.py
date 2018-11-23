#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Tue Nov  1 23:02:39 2016.

@author: lashkov

"""

from urllib import error, request
import json

__version__ = '0.2.0'
__license__ = 'GPLv3'
try:
    from sklearn import __version__ as __skversion__
except ImportError:
    raise ImportError('Error! Scikit-learn not installed!')


def newversioncheck():
    try:
        download_version = json.load(request.urlopen("https://pypi.org/pypi/hydrocluster/json"))['info']['version']
    except (request.HTTPError, error.URLError):
        download_version = "Not available"
    return download_version


def tohexversion(ver: str):
    major, minor, revision = ver.split(',')
    hexvariant = hex(major)*0x1000 + hex(minor)*0x100 + hex(revision)*0x10
    return hexvariant

