#!/usr/bin/env python3
# coding=utf-8
from setuptools import setup

setup(
    name='hydrocluster',
    version='0.1',
    python_requires='>=3.4',
    packages=['hydrocluster'],
    url='https://github.com/alashkov83/hydrocluster',
    license='GPL v.3',
    platforms=['any'],
    author='Aleksandr Lashkov',
    author_email='alashkov83@gmail.com',
    maintainer='Aleksandr Lashkov',
    classifiers=['Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4'],
    description='',

    install_requires=['scipy', 'progressbar', 'python-utils', 'psutil', 'progressbar2', 'matplotlib', 'numpy',
                      'scikit_learn', 'biopython', 'mmtf-python', 'msgpack'],
    entry_points={'console_scripts': ['hydrocluster = hydrocluster.hydrocluster:main']})
