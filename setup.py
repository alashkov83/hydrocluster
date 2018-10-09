#!/usr/bin/env python3
# coding=utf-8
from setuptools import setup

with open("README.md") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.readlines()

setup(
    name='hydrocluster',
    version='0.1.0',
    python_requires='>=3.4',
    packages=['hydrocluster', 'hydrocluster.ui', 'hydrocluster.core', 'hydrocluster.dbcreator'],
    url='https://github.com/alashkov83/hydrocluster',
    license='GPL v.3',
    platforms=['any'],
    author='Aleksandr Lashkov',
    author_email='alashkov83@gmail.com',
    maintainer='Aleksandr Lashkov',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: MacOS',
        'Environment :: Win32 (MS Windows)',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    description='Cluster analysis of hydrophobic or charged regions of macromolecules',
    long_description=long_description,
    keywords='',
    long_description_content_type="text/markdown",
    install_requires=requirements,
    entry_points={'console_scripts': ['hydrocluster = hydrocluster.hydrocluster:main']})
