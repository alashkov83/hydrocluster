# HYDROCLUSTER

## Requirements:
Python 3.4 or higher (CPython only support)
psutil
progressbar2
matplotlib>=1.5.1  
numpy>=1.14.2  
scikit_learn>=0.19.1  
biopython>=1.71  
mmtf-python>=1.1.0  
msgpack>=0.5.6

To easily browse through db files you will need a DB Browser for SQLite (<https://sqlitebrowser.org>).
It is recommended to install Pymol molecular viewer (version: 1.7+).


##### For MS Windows:
Use Anaconda (<https://anaconda.org>) for Windows - it includes majority of the dependencies required.
But mmtf-python and msgpack not available on Anaconda - need to use pip.
Define environment variable PYTHONIOENCODING to UTF-8.
For correct display of the Angstrom symbol use console fonts including this symbol
(for example, SimSun font family).
