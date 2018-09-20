#!/usr/bin/env python
# coding=utf-8
from setuptools import setup, find_packages

setup(
    name='hydrocluster',
    version='1.0',
    packages=find_packages(exclude=['Publication', "PDB_LISTS", 'PDB_TESTS']),
    url='https://github.com/alashkov83/hydrocluster',
    license='GPL v.3',
    author='Aleksandr Lashkov',
    maintainer='Aleksandr Lashkov',
    author_email='alashkov83@gmail.com',
    classifiers=['Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4'],
    description='',
    install_requires=['scipy', 'progressbar', 'python-utils', 'psutil', 'progressbar2', 'matplotlib', 'numpy',
                      'scikit_learn', 'biopython', 'mmtf-python', 'msgpack'])
