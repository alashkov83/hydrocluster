#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Tue Nov  1 23:02:39 2016.

@author: lashkov

"""

import json
import re
from urllib import error, request

__version__ = '0.2.0d56'
__license__ = 'GPLv3'
try:
    from sklearn import __version__ as __skversion__
except ImportError:
    raise ImportError('Error! Scikit-learn not installed!')


def newversioncheck() -> str:
    try:
        download_version = json.load(request.urlopen("https://pypi.org/pypi/hydrocluster/json"))['info']['version']
    except (request.HTTPError, error.URLError):
        download_version = "Not available"
    return download_version


def tohexversion(ver: str) -> hex:
    """
    :param ver:: Version in string format: a.b.cde, where a - integer number of the major version (0-255),
    b - integer number of the minor version (0-255),
    c - integer number of the revision (0-255),
    d - char of the developer status (d - developer, a -alpha, b - beta) for release version this field is missing,
    e - integer number of variant of dev status (optional, for release field is missing) (0-63)
    :return: Version in hex number: 0xAABBCCDD, where AA - hex number of the major version (0-FF),
    BB - hex number of the minor version (00-FF),
    CC - hex number of the revision (00-FF),
    DD - base + hex number of variant of dev status (base: 0x00 - dev, 0x80 - alpha, 0xC0 - beta, 0xFF - release)
    """

    smajor, sminor, srevision, sdev, sdevn = re.match(r'^(\d+)\.(\d+)\.(\d+)([a-z]?)(\d*)$', ver).groups()
    if sdev == 'd':
        base = 0x0
    elif sdev == 'a':
        base = 0x80
    elif sdev == 'b':
        base = 0xc0
    else:
        base = 0xff
    devn = int(sdevn) if sdevn else 0x0
    hexvariant = int(smajor) * 0x1000000 + int(sminor) * 0x10000 + int(srevision) * 0x100 + base + devn
    return hexvariant
