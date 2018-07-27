#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import os
import shutil
import sys
from urllib.error import HTTPError

from matplotlib.backends.backend_agg import FigureCanvasAgg

try:
    from .pdbcluster import ClusterPdb
except ImportError:
    print('Error! Scikit-learn not installed!')
    sys.exit()
try:
    import progressbar2 as progressbar
except ImportError:
    import progressbar


class Cli:
    """

    """
    def __init__(self, namespace) -> None:
        self.namespace = namespace
        if namespace.output:
            newdir = namespace.output
            basefile = namespace.output
        else:
            if namespace.input.count('.') > 0:
                newdir = namespace.input.split('.')[-2]
                basefile = os.path.basename(namespace.input).split('.')[-2]
            else:
                newdir = namespace.input
                basefile = namespace.input
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        try:
            os.makedirs(newdir, exist_ok=True)
        except OSError:
            print('Unable to create folder ' + newdir)
            sys.exit(-1)
        self.log_name = os.path.join(newdir, '{:s}'.format(basefile + '.log'))
        self.cls = ClusterPdb()
        self.open_file(namespace.input)
        self.parse_pdb(namespace.ptable, namespace.pH)
        if namespace.noauto:
            self.noauto(namespace.eps, namespace.min_samples)
        else:
            self.run()
            self.save_state(newdir, basefile)
            self.colormap(newdir, basefile)
        self.resi()
        self.save_pymol(newdir, basefile)
        self.graph(newdir, basefile)

    def log_append(self, line: str):
        """

        :param line:
        """
        print(line, end='')
        with open(self.log_name, 'at',  encoding='utf-8') as f:
            f.write(line)

    def run(self) -> None:
        """The main algorithm of the program."""

        metric = self.namespace.score
        min_eps = self.namespace.emin
        max_eps = self.namespace.emax
        step_eps = self.namespace.estep
        min_min_samples = self.namespace.smin
        max_min_samples = self.namespace.smax
        self.log_append(('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                         'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
            min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
        bar1 = progressbar.ProgressBar(maxval=self.cls.init_cycles(
            min_eps, max_eps, step_eps, min_min_samples, max_min_samples), redirect_stdout=True).start()
        #       import time
        #       start_time = time.time()
        try:
            for n, j, i in self.cls.auto_yield():
                self.log_append(('Step No. {0:d}: EPS = {1:.2f} \u212B, min_samples = {2:d}, No. of clusters = {3:d}, '
                                 'Silhouette score = {4:.3f}; Calinski score = {5:.3f}\n').format(
                    n, j, i, self.cls.n_clusters, self.cls.si_score, self.cls.calinski))

                bar1.update(n)
            #           print("--- %s seconds ---" % (time.time() - start_time))
            eps, min_samples = self.cls.auto(metric=metric)
        except ValueError:
            self.log_append('Error! Could not parse file or clustering failed\n')
            sys.exit(-1)
        else:
            self.log_append('Autoscan done... \n')
            bar1.finish()
            self.log_append(('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                             'Calinski-Harabaz score = {4:.3f}\n'
                             'EPS = {2:.1f} \u212B\nMIN_SAMPLES = {3:d}\n'
                             'Percent of noise = {5:.2f} %{6:s}\n').format(
                self.cls.n_clusters, self.cls.si_score, eps, min_samples, self.cls.calinski, self.cls.noise_percent(),
                (' !!WARNING!!!' if self.cls.noise_percent() > 30 else '')))

    def noauto(self, eps: float, min_samples: int):
        """

        :param eps:
        :param min_samples:
        """
        if min_samples <= 0 or eps <= 0:
            print("--eps and --min samples options are required with values > 0")
            sys.exit(-1)
        try:
            self.cls.cluster(eps, min_samples)
        except ValueError:
            self.log_append('Error! Could not parse file or clustering failed\n')
            sys.exit(-1)
        else:
            self.log_append(('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                             'Calinski-Harabaz score = {4:.3f}\nEPS = {2:.1f} \u212B\n'
                             'MIN_SAMPLES = {3:d}\nPercent of noise = {5:.2f} %\n').format(
                self.cls.n_clusters, self.cls.si_score, eps, min_samples, self.cls.calinski, self.cls.noise_percent()))

    def graph(self, newdir: str, basefile: str):
        """

        :param newdir:
        :param basefile:
        :return:
        """
        grid, legend = True, True
        sa = os.path.join(newdir, '{:s}'.format(basefile + '.png'))
        try:
            fig, ax = self.cls.graph(grid, legend)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Error! Failed to plot!\n')

    def open_file(self, filename: str):
        """

        :param filename:
        """
        if not filename:
            self.log_append('Filename not defined\n')
            sys.exit(-1)
        if filename.split('.')[-1].strip().lower() == 'pdb':
            try:
                self.cls.open_pdb(filename)
            except FileNotFoundError:
                self.log_append('File {:s} not found!\n'.format(filename))
                sys.exit(-1)
            else:
                self.log_append('File {:s} successfully read!\n'.format(filename))
        elif filename.split('.')[-1].strip().lower() == 'cif':
            try:
                self.cls.open_cif(filename)
            except ImportError:
                self.log_append('Import error! BioPython not found! \nPlease install biopython and mmtf!\n')
                sys.exit(-1)
            except FileNotFoundError:
                sys.exit(-1)
            except ValueError:
                self.log_append('Error! Incorrect CIF file: {0:s}!\n'.format(filename))
                sys.exit(-1)
            else:
                self.log_append('File {:s} successfully read!\n')
        else:
            try:
                self.cls.open_url(filename)
            except ImportError:
                self.log_append('Import error! BioPython unavailable! \nPlease install biopython and mmtf\n')
                sys.exit(-1)
            except HTTPError as e:
                self.log_append('Error! ID PDB: {0:s} not found or refers to an incorrect file!'.format(
                    filename, str(e)))
                sys.exit(-1)
            else:
                self.log_append('File ID PDB: {0:s} successfully downloaded!\n'.format(filename))

    def parse_pdb(self, htable: str, pH: float):
        """

        :param htable:
        :param pH:
        :return:
        """
        if htable == 'positive' or htable == 'negative':
            if pH < 0 or pH > 14:
                print("pH value range is 0-14")
                sys.exit(-1)
        try:
            parse_results = self.cls.parser(htable=htable, pH=pH)
        except ValueError:
            self.log_append('Error! Invalid file format\nor file does not contain {:s} residues\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'nanodroplet')
                else 'negative' if htable == 'negative' else 'positive'))
        else:
            self.log_append("No. of residues: {:d}\nMinimum distance = {:.3f} \u212B\n"
                            "Maximum distance = {:.3f} \u212B\nMean distance = {:.3f} \u212B\n".format(*parse_results))

    def save_state(self, newdir: str, basefile: str):
        """

        :param newdir:
        :param basefile:
        :return:
        """
        st = os.path.join(newdir, '{:s}'.format(basefile + '.dat'))
        try:
            self.cls.savestate(st)
        except FileNotFoundError:
            return

    def save_pymol(self, newdir: str, basefile: str):
        """

        :param newdir:
        :param basefile:
        :return:
        """
        pymol = os.path.join(newdir, '{:s}'.format(basefile + '.py'))
        try:
            self.cls.save_pymol_script(pymol)
        except FileNotFoundError:
            return

    def resi(self):
        """

        :return:
        """
        dict_aa = self.cls.get_dict_aa()
        if not dict_aa:
            return
        for k, aa_list in dict_aa.items():
            self.log_append('\n{:s} cluster No. {:d} contains: {:s}'.format(
                ("Core" if k[0] else "Non-core"), k[1], ", ".join(['{2:s}:{1:s}{0:d}'.format(*aac) for aac in aa_list])))
        self.log_append('\n\n')

    def colormap(self, newdir: str, basefile: str):
        """

        :param newdir:
        :param basefile:
        :return:
        """
        sa = os.path.join(newdir, '{:s}'.format(basefile + '.cm.png'))
        try:
            fig = self.cls.colormap(grid_state=True)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Error! Failed to plot!!\n')
