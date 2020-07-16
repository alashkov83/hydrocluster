#!/usr/bin/env python3
# coding=utf-8
import re

from setuptools import setup

with open("README.md") as fh:
    long_description = fh.read()
    long_description = re.sub(r'!\[image\]\([a-z0-9/_]+\.png\)\n', '', long_description)

with open('requirements.txt') as f:
    requirements = f.readlines()

setup(
    name='hydrocluster',
    version='0.3.0d2',
    python_requires='>=3.4',
    packages=['hydrocluster', 'hydrocluster.ui', 'hydrocluster.core', 'hydrocluster.dbcreator'],
    url='https://github.com/alashkov83/hydrocluster',
    license='GPL v.3',
    platforms=['any'],
    author='Alexander Lashkov, Sergey Rubinsky, Polina Eistrikh-Heller',
    author_email='alashkov83@gmail.com',
    maintainer='Alexander Lashkov',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Environment :: Win32 (MS Windows)',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    description='Cluster analysis of hydrophobic or charged regions of macromolecules. '
                'The program is based on the DBSCAN algorithm.',
    long_description=long_description,
    keywords='molecular modeling, bioinformatic, protein structure, hydrophobic core, hydrophobic cluster, DBSCAN',
    long_description_content_type="text/markdown",
    install_requires=requirements,
    entry_points={'console_scripts': ['hydrocluster = hydrocluster.hydrocluster:main']})
