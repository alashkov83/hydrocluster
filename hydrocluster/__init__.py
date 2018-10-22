#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Tue Nov  1 23:02:39 2016.

@author: lashkov

"""

__version__ = '0.2.0'
__license__ = 'GPLv3'
try:
    from sklearn import __version__ as __skversion__
except ImportError:
    raise ImportError('Error! Scikit-learn not installed!')
