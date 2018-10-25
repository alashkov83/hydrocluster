Hydrocluster - Bioimformatics Tool
==================================

Short description
-----------------
The program Hydrocluster is designed to determine the position, size and content of hydrophobic,
hydrophilic and charged clusters in protein molecules. The program is based on the DBSCAN algorithm.

**Keywords:** molecular modeling, bioinformatic, protein structure,
hydrophobic core, hydrophobic cluster, DBSCAN

Installation
------------

```shell
pip install --upgrade hydrocluster
```
(or pip3 in distributive with default python 2 version)

User Interface
--------------

### Command line

The program is called with the command ’hydrocluster’ and following
parameters:
```shell
hydrocluster [-h][-i INPUT][-emin EMIN][-emax EMAX][-es ESTEP]
[-smin SMIN][-smax SMAX][-g {tkgui,cli,testlist}][-o OUTPUT][-c CHAINS]
[-rl RESLIST][-pt{hydropathy,menv,fuzzyoildrop,nanodroplet,aliphatic_core,hydrophilic,positive,negative}]
[-pH PH][-sc {si_score,calinski,dbcv}][-nf][-na][-eps EPS][-min_samples MIN_SAMPLES]
```
#### Arguments:

**-h, --help**  
show help message and exit

**-i INPUT, --input INPUT**  
Input file name (pdb, txt, cif, ent) - pdb file name, cif file name,
idpdb or list of ispdbs

**-emin EMIN, --emin EMIN**  
Minimum EPS value (A). Default=3.0

**-emax EMAX, --emax EMAX**  
Maximum EPS value (A). Default=15.0

**-es ESTEP, --estep ESTEP**  
Step of EPS (A). Default=0.1

**-pts,--ptables**  
Property table’s list for testlist.py. Separator: ’,’

**-scs, --scores**  
Score coefficients list for testlist. Separator: ’,’

**-smin SMIN, --smin SMIN**  
Minimum MIN SAMPLES. Default=3

**-ss --save\_state**  
Save states on testlist.py

**-smax SMAX, --smax SMAX**  
Minimum MIN SAMPLES. Default=50

**-g {tkgui,cli,testlist}, --gui**  
UI modes. Default=’tkgui’ (tkgui - graphic interface, cli - command
line, testlist - using testlist module for data processing (see -i
filename.txt and -o filename of data base).

**-o OUTPUT, --output OUTPUT**  
Output directory name/file name or db name

**-c CHAINS, --chains CHAINS**  
Selected chains. Default=None

**-rl RESLIST, --reslist RESLIST**  
Selected amino acid residues. Default=None

**-pt{hydropathy,menv,fuzzyoildrop,nanodroplet,aliphatic\_core,hydrophilic,positive,negative}, --ptable**  
Property table for weighting. Default=’hydropathy’

**-pH PH**  
pH value for calculatation of net charges (positive or negative). Default=7.0

**-sc {si\_score,calinski,dbcv}, --score {si\_score,calinski,dbcv}**  
Score coefficient. Default=’calinski’

**-nf, --noise\_filter**  
Activate filter of noise for scoring function (**Not recommended!!!**).

**-na, --noauto**  
No automatic mode.

**-eps EPS**  
EPS value (A). Default=3.0

**-min\_samples MIN\_SAMPLES**  
MIN SAMPLES value. Default=3

**At startup of hydrocluster without any parameters the program opens
with graphics interface.**

### Examples:
```shell
hydrocluster -i subpdb.txt -g testlist -o subpdb
```

Processing of file\_name.txt by testlist.py, file\_name.db file and
file\_name\_data folder consisting of tree structure with data files on
return

```shell
hydrocluster -i 1atg.pdb -g cli -o 1atg
```
Processing of file\_name.pdb by command line interface and file\_name
folder on return

File\_name folder consists of file\_name.py file for processing by
pymol, binary file (.dat) with saved session state, file\_name.log file
with saved log-data and two png files with pictures.

Graphical User Interface
------------------------

GUI was realized using Tkinter. It consists of a panel for selecting the
operation mode, a window for graphical representation of clustering
results Cluster analysis, and a window for displaying log file.

![image](screenshots/main_window.png)

At the beginning of working with the graphical interface, it is
necessary to select the desired hydrophobicity/hydrophilicity table in
the sub-window of the mode selection window, select the method for
scoring of clustering in the metrics window and run on Manual (Manual mode -&gt;
Start) or automatic mode of operation (Auto mode -&gt; Start) in one of
the underlying windows. In the automatic mode, the optimal parameters
eps and min\_samples are selected by enumeration within the given
boundaries and with the given step. Upon completion of the work in the
automatic mode, when you click Options -&gt; Dispplay colormap, you can
get a graphical interpretation of the process of selecting the optimal
values namely dependencies min\_samples (eps) and min\_samples (eps<sup>3</sup>).
The point corresponding to the optimal parameters is marked in
color.

![image](screenshots/colormap.png)

The Cluster analysis window presents a three-dimensional image of
clusters selected by the program in a protein molecule. Appropriate mtnu
sections allow you to make a coordinate grid in the image and get a brief comment
on the picture.

The Log window shows the numerical results of clustering, namely the
number of chains and clusters, the percentage of noise and the optimal
values of the hyperparameters (eps,min\_samples) and the metric used.
Further study of the macromolecule can be carried out using the PyMol
program (Options-&gt; OpenPyMol).

![image](screenshots/pymol_view.png)

### Menu options:

**File-&gt;**

Open PDB - opens pdb file on a disk  
Open CIF - opens CIF file on a disk  
Open IDPDB - opens file from RSCB PDB data base with IDPDB  
Load state - loads program state, saved in file  
Save PyMOL script - saves script (.py) for further processing with PyMOL  
Save state - saves the current state of program in file  
Save picture - saves the clustering result in png format file  
Save LOG - saves log file of the current session  
Quit - quit from the program

**Options-&gt;**

Plot grid - makes coordinate grid in the Cluster analysis window  
Plot legend - displays the brief description of the picture  
Select clustering solution - display and choice for other solutions cluster analysis  
Display colormap - shows graphs obtained as a result of clustering
parameters selection. Marked point corresponds optimal values of eps and min\_samples  
Clear log - clears log file in the appropriate window  
Open PyMol - opens PyMol for further data display

Theory
------

Hydrophobic cores and hydrophobic clusters play an important role in the
folding of the protein, being the skeleton for functionally important
amino acid residues of enzyme proteins. In the cases of ligands of
amphiphilic nature, the hydrophobic clusters themselves are included in
the functionally important regions of the molecules. The interaction
with them should be taken into account, for example, when evaluating
molecular docking solutions. Hydrocluster programm is based on
ensity-Based Spatial Clustering of Applications with Noise (DBSCAN).
Atomic coordinates, their type and description of amino acid residues
(a. r.) are loaded from a file of the PDB, CIF formats, or directly from
the Protein Data Bank. For each a.r. from the table of relative
hydrophobicity center of mass of non-H atoms is calculated. As weights
in the cluster analysis, various tables of a.r. hydrophobicity known in
the literature are used. \[1-4]. Separately, for clustering
electrically charged amino acid residues, the function of calculating
weighting coefficients as modules of partial charges of side groups
according to the formulas, which are derived from the
Henderson-Hasselbach equation, is implemented \[5]. As
hyperparameters DBSCAN uses the epsilon neighborhood radius (eps) and
the minimum number of neighbors (min\_samples). Eps is defined as the
maximum distance (in Angstrom) between the centers of mass of
hydrophobic a.r. which are adjacent in one cluster. The
min\_samples/eps\^{3} ratio is proportional to the maximum distribution
density of the centers of mass of the hydrophobic a.r.. Silhouette score
\[6] and Calinski and Harabaz score \[7] and DBCV
\[8] are used as the quality criteria for cluster analysis. For
clusters of complex shape, it is better to use the silhouette
coefficient. At the same time, Calinski and Harabaz score, which uses
the distance between the element and the center of the cluster,
correctly estimates the areas of clusters with the highest density. This
areas are of interest from the point of view of the structural
organization of proteins. A feature of the DBSCAN algorithm is the
strong dependence of clustering results on the parameters - eps and
min\_samples. Hydrocluster implemented the selection of these parameters
by simply iterating over their values at user-defined boundaries,
followed by sorting the results according to the criterion of maximizing
the value of the corresponding estimated coefficient.

Requirements
------------

* Python 3.4 or higher (CPython only support)
* psutil
* progressbar2
* matplotlib>=1.5.1
* numpy>=1.14.2
* scikit_learn>=0.19.1
* biopython>=1.71
* mmtf-python>=1.1.0
* msgpack>=0.5.6

To easily browse through db files you will need a DB Browser for SQLite
(<https://sqlitebrowser.org>). It is recommended to install Pymol
molecular viewer (version: 1.7+).

**For MS Windows:** Use Anaconda (<https://anaconda.org>) for Windows -
it includes majority of the dependencies required. But mmtf-python and
msgpack not available on Anaconda - need to use pip. Define environment
variable PYTHONIOENCODING to UTF-8. For correct display of the Angstrom
symbol use console fonts including this symbol (for example, SimSun font
family).


References
----------
1. J. Kyte, R. F. Doolittle. J Mol Biol. 1982. 157, 105-132.
2. B. Kalinowska, M. Banach, Z. Wisniowski, L. Konieczny, I. Roterman. J
Mol Model. 2017. 23 , 205.
3. D. Bandyopadhyay .E. L. Mehler.Proteins 2008.72.646-659
4. Zhu C. Q., Gao Y. R. , Li H. et.al.// Proc. NAS. 2016.113.12946.
5. Ikai, A.J. 1980. J. Biochem. 88, 1895-1898.
6. Rousseeuw P. Comput. Appl. Math. 1987. 20. 53.
7. Calinski T., Harabasz J. // Communications in Statistics. 1974. 3 . 1.
8. Moulavi, Davoud, et al. Proc.2014 SIAM International Conf. on Data
Mining. 839-847. 2014.
