#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import os
import re
import shutil
import sys
import time
from urllib.error import HTTPError

from matplotlib.backends.backend_agg import FigureCanvasAgg

try:
    from ..core.pdbcluster import ClusterPdb
except ImportError:
    print('Error! Scikit-learn not installed!')
    sys.exit()
import progressbar as progressbar


class Log:
    """

    """

    def __init__(self, log_name):
        self.log_name = log_name

    def append(self, line: str, ascitime=False):
        """
    
        :param ascitime:
        :param line:
        """
        if ascitime:
            line = time.asctime() + ': ' + line
        print(line, end='')
        with open(self.log_name, 'at', encoding='utf-8') as f:
            f.write(line)


def run(cls: ClusterPdb, log: Log, namespace) -> None:
    """The main algorithm of the program."""

    metric = namespace.score
    min_eps = namespace.emin
    max_eps = namespace.emax
    step_eps = namespace.estep
    min_min_samples = namespace.smin
    max_min_samples = namespace.smax
    log.append(('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
        min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
    bar1 = progressbar.ProgressBar(maxval=cls.init_cycles(
        min_eps, max_eps, step_eps, min_min_samples, max_min_samples, metric=metric), redirect_stdout=True).start()
    try:
        for n, j, i, n_clusters, score in cls.auto_yield():
            log.append(('Step No. {0:d}: EPS = {1:.2f} \u212B, min_samples = {2:d}, No. of clusters = {3:d}, '
                        '{4:s} = {5:.3f}\n').format(
                n, j, i, n_clusters, cls.metrics_name[cls.metric], score))
            bar1.update(n)
        eps, min_samples = cls.auto()
    except ValueError:
        log.append('Error! Could not parse file or clustering failed\n')
        sys.exit(-1)
    else:
        log.append('Autoscan done... \n')
        bar1.finish()
        log.append(('\nNumber of clusters = {0:d}\n{4:s} = {1:.3f}\n'
                    'EPS = {2:.3f} \u212B\nMIN_SAMPLES = {3:d}\n'
                    'Percent of noise = {5:.2f} %{6:s}\n').format(
            cls.n_clusters, cls.score, eps, min_samples,
            cls.metrics_name[cls.metric], cls.noise_percent(),
            (' !!WARNING!!!' if cls.noise_percent() > 30 else '')))


def noauto(cls: ClusterPdb, log: Log, eps: float, min_samples: int, metric: str):
    """

    :param log:
    :param cls:
    :param metric:
    :param eps:
    :param min_samples:
    """
    if min_samples <= 0 or eps <= 0:
        print("--eps and --min samples options are required with values > 0")
        sys.exit(-1)
    try:
        cls.cluster(eps, min_samples, metric=metric)
    except ValueError:
        log.append('Error! Could not parse file or clustering failed\n')
        sys.exit(-1)
    else:
        log.append(('\nNumber of clusters = {0:d}\n{4:s} = {1:.3f}\n'
                    'EPS = {2:.3f} \u212B\nMIN_SAMPLES = {3:d}\n'
                    'Percent of noise = {5:.2f} %{6:s}\n').format(cls.n_clusters, cls.score, eps, min_samples,
                                                                  cls.metrics_name[cls.metric], cls.noise_percent(),
                                                                  ' !!WARNING!!!' if cls.noise_percent() > 30 else ''))


def graph(cls: ClusterPdb, log: Log, newdir: str, basefile: str):
    """

    :param log:
    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    grid, legend = True, True
    sa = os.path.join(newdir, '{:s}'.format(basefile + '.png'))
    try:
        fig, ax = cls.graph(grid, legend)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError:
        log.append('Error! Failed to plot!\n')


def open_file(cls: ClusterPdb, log: Log, filename: str):
    """

    :param log:
    :param cls:
    :param filename:
    """
    if not filename:
        log.append('Filename not defined\n')
        sys.exit(-1)
    if filename.split('.')[-1].strip().lower() == 'pdb':
        try:
            cls.open_pdb(filename)
        except FileNotFoundError:
            log.append('File {:s} not found!\n'.format(filename))
            sys.exit(-1)
        else:
            log.append('File {:s} successfully read!\n'.format(filename))
    elif filename.split('.')[-1].strip().lower() == 'cif':
        try:
            cls.open_cif(filename)
        except ImportError:
            log.append('Import error! BioPython not found! \nPlease install biopython and mmtf!\n')
            sys.exit(-1)
        except FileNotFoundError:
            sys.exit(-1)
        except ValueError:
            log.append('Error! Incorrect CIF file: {0:s}!\n'.format(filename))
            sys.exit(-1)
        else:
            log.append('File {:s} successfully read!\n')
    else:
        try:
            cls.open_url(filename)
        except ImportError:
            log.append('Import error! BioPython unavailable! \nPlease install biopython and mmtf\n')
            sys.exit(-1)
        except HTTPError as e:
            log.append('Error! ID PDB: {0:s} not found or refers to an incorrect file!'.format(
                filename, str(e)))
            sys.exit(-1)
        else:
            log.append('File ID PDB: {0:s} successfully downloaded!\n'.format(filename))


def chainsSelect(cls: ClusterPdb, log: Log, namespace):
    """

    :param log:
    :param cls:
    :param namespace:
    :return:
    """
    if namespace.chains is None:
        log.append("All chains selected!\n")
        return None
    selectChains = re.findall(r'[A-Z]', namespace.chains.strip().upper())
    allChains = cls.preparser()
    if set(selectChains).issubset(allChains):
        log.append("Selected chains: {:s}\n".format(', '.join(selectChains)))
        return selectChains
    else:
        log.append("Error! Chain(s): {:s} is not include in structure! Structure included: {:s} chain(s)!\n".format(
            ', '.join(set(selectChains).difference(set(allChains))), ','.join(allChains)))
        sys.exit(-1)


def parse_pdb(cls: ClusterPdb, log: Log, htable: str, pH: float, chains: list = None, res=''):
    """

    :param log:
    :param cls:
    :param res:
    :param htable:
    :param pH:
    :param chains:
    :return:
    """
    if htable == 'positive' or htable == 'negative':
        if pH < 0 or pH > 14:
            print("pH value range is 0-14")
            sys.exit(-1)
    try:
        parse_results = cls.parser(htable=htable, pH=pH, selectChains=chains, res=res)
    except ValueError:
        log.append('\nError! Invalid file format\nor file does not contain {:s} residues\n'.format(
            'hydrophobic' if htable in ('hydropathy', 'nanodroplet', 'menv', 'fuzzyoildrop')
            else 'negative' if htable == 'negative' else 'positive'))
    else:
        log.append("No. of residues: {:d}\nMinimum distance = {:.3f} \u212B\n"
                   "Maximum distance = {:.3f} \u212B\nMean distance = {:.3f} \u212B\n".format(*parse_results))


def save_state(cls: ClusterPdb, newdir: str, basefile: str):
    """

    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    st = os.path.join(newdir, '{:s}'.format(basefile + '.dat'))
    try:
        cls.savestate(st)
    except FileNotFoundError:
        return


def save_pymol(cls: ClusterPdb, newdir: str, basefile: str):
    """

    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    pymol = os.path.join(newdir, '{:s}'.format(basefile + '.py'))
    try:
        cls.save_pymol_script(pymol)
    except FileNotFoundError:
        return


def resi(cls: ClusterPdb, log: Log):
    """

    :return:
    """
    dict_aa = cls.get_dict_aa()
    if not dict_aa:
        return
    for k, aa_list in dict_aa.items():
        log.append('\n{:s} cluster No. {:d} contains: {:s}'.format(
            ("Core" if k[0] else "Non-core"), k[1],
            ", ".join(['{2:s}:{1:s}{0:d}'.format(*aac) for aac in aa_list])))
    log.append('\n\n')


def colormap(cls: ClusterPdb, log: Log, newdir: str, basefile: str):
    """

    :param log:
    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    sa = os.path.join(newdir, '{:s}'.format(basefile + '.cm.png'))
    try:
        fig = cls.colormap(grid_state=True)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError:
        log.append('Error! Failed to plot!!\n')


def cli(namespace) -> None:
    """

    :param namespace:
    """
    if namespace.output:
        inpt = namespace.output
    else:
        inpt = os.path.splitext(namespace.input)[0]
    newdir = os.path.join(os.path.dirname(inpt), os.path.basename(inpt))
    basefile = os.path.basename(inpt)
    if os.path.exists(newdir):
        shutil.rmtree(newdir)
    try:
        os.makedirs(newdir, exist_ok=True)
    except OSError:
        print('Unable to create folder ' + newdir)
        sys.exit(-1)
    log = Log(os.path.join(newdir, '{:s}'.format(basefile + '.log')))
    cls = ClusterPdb()
    open_file(cls, log, namespace.input)
    parse_pdb(cls, log, namespace.ptable, namespace.pH, chainsSelect(cls, log, namespace), namespace.reslist)
    cls.noise_filter = namespace.noise_filter
    if namespace.noauto:
        noauto(cls, log, namespace.eps, namespace.min_samples, namespace.score)
    else:
        run(cls, log, namespace)
        save_state(cls, newdir, basefile)
        colormap(cls, log, newdir, basefile)
    resi(cls, log)
    save_pymol(cls, newdir, basefile)
    graph(cls, log, newdir, basefile)
