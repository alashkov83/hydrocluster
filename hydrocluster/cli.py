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
    print('Error! Scikit-learn is not installed!')
    sys.exit()
try:
    import progressbar2 as progressbar
except ImportError:
    import progressbar


class Cli:
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
        self.graph(newdir, basefile)

    def log_append(self, line):
        print(line, end='')
        with open(self.log_name, 'at') as f:
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
        try:
            for n, j, i in self.cls.auto_yield():
                self.log_append(('Step No {0:d}: EPS = {1:.2f} \u212B, min_samples = {2:d}, No clusters = {3:d}, '
                                 'Silhouette score = {4:.3f} Calinski score = {5:.3f}\n').format(
                    n, j, i, self.cls.n_clusters, self.cls.si_score, self.cls.calinski))

                bar1.update(n)
            eps, min_samples = self.cls.auto(metric=metric)
            self.log_append('Autoscan done... \n')
        except ValueError:
            self.log_append('Error! File was not parse or clustering was fail\n')
            sys.exit(-1)
        else:
            self.log_append(('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                             'Calinski-Harabaz score = {4:.3f}\nEPS = {2:.1f} \u212B\nMIN_SAMPLES = {3:d}\n').format(
                self.cls.n_clusters, self.cls.si_score, eps, min_samples, self.cls.calinski))

    def noauto(self, eps, min_samples):
        if min_samples <= 0 or eps <= 0:
            print("--eps and --min samples options are required and it's values > 0")
            sys.exit(-1)
        try:
            self.cls.cluster(eps, min_samples)
        except ValueError:
            self.log_append('Error! File was not parse or clustering was fail\n')
            sys.exit(-1)
        else:
            self.log_append(('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                             'Calinski-Harabaz score = {4:.3f}\nEPS = {2:.1f} \u212B\nMIN_SAMPLES = {3:d}\n').format(
                self.cls.n_clusters, self.cls.si_score, eps, min_samples, self.cls.calinski))

    def graph(self, newdir, basefile):
        grid, legend = True, True
        sa = os.path.join(newdir, '{:s}'.format(basefile + '.png'))
        try:
            fig, ax = self.cls.graph(grid, legend)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Error! Plot was not created!\n')
            return

    def open_file(self, filename):
        if not filename:
            self.log_append('Filename was not defined\n')
            sys.exit(-1)
        if filename.split('.')[-1].strip().lower() == 'pdb':
            try:
                self.cls.open_pdb(filename)
            except FileNotFoundError:
                self.log_append('File {:s} not found!\n'.format(filename))
                sys.exit(-1)
            else:
                self.log_append('File {:s} was read!\n'.format(filename))
        elif filename.split('.')[-1].strip().lower() == 'cif':
            try:
                self.cls.open_cif(filename)
            except ImportError:
                self.log_append('Import error! Bio Python unavailable! \nPlease install biopython and mmtf!\n')
                sys.exit(-1)
            except FileNotFoundError:
                sys.exit(-1)
            except ValueError:
                self.log_append('Error! Incorrect CIF file: {0:s}!\n'.format(filename))
                sys.exit(-1)
            else:
                self.log_append('File {:s} read!\n')
        else:
            try:
                self.cls.open_url(filename)
            except ImportError:
                self.log_append('Import error! Bio Python unavailable! \nPlease install biopython and mmtf\n')
                sys.exit(-1)
            except HTTPError as e:
                self.log_append('Error! ID PDB: {0:s} not found or refers to an incorrect file!'.format(
                    filename, str(e)))
                sys.exit(-1)
            else:
                self.log_append('File ID PDB: {0:s} was downloaded!\n'.format(filename))

    def parse_pdb(self, htable, pH):
        if htable == 'positive' or htable == 'negative':
            if pH < 0 or pH > 14:
                print("pH value range is 0-14")
                sys.exit(-1)
        try:
            parse_results = self.cls.parser(htable=htable, pH=pH)
        except ValueError:
            self.log_append('Error! Invalid file format\nor file does not {:s} contain resides\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'nanodroplet')
                else 'negative' if htable == 'negative' else 'positive'))
            return
        else:
            self.log_append("No of residues: {:d}\nMinimum distance = {:.3f} \u212B\n"
                            "Maximum distance = {:.3f} \u212B\nMean distance = {:.3f} \u212B\n".format(*parse_results))

    def save_state(self, newdir, basefile):
        st = os.path.join(newdir, '{:s}'.format(basefile + '.dat'))
        try:
            self.cls.savestate(st)
        except FileNotFoundError:
            return

    def resi(self):
        aa_list = self.cls.aa_list
        if (not aa_list) or self.cls.labels is None:
            return
        aa_list = list(zip(aa_list, self.cls.labels, self.cls.core_samples_mask))
        for k in sorted(set(self.cls.labels)):
            if k != -1:
                core_aa_list = []
                uncore_aa_list = []
                for aa in aa_list:
                    if aa[1] == k and aa[2]:
                        core_aa_list.append(aa[0])
                    elif aa[1] == k and not aa[2]:
                        uncore_aa_list.append(aa[0])
                self.log_append('\nIn Core cluster No {:d} included: '.format(k + 1))
                for aac in core_aa_list:
                    self.log_append('{2:s}:{1:s}{0:d} '.format(*aac))
                if uncore_aa_list:
                    self.log_append('\nIn UNcore cluster No {:d} included: '.format(k + 1))
                    for aac in uncore_aa_list:
                        self.log_append('{2:s}:{1:s}{0:d} '.format(*aac))
        self.log_append('\n\n')

    def colormap(self, newdir, basefile):
        sa = os.path.join(newdir, '{:s}'.format(basefile + '.cm.png'))
        try:
            fig = self.cls.colormap(grid_state=True)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Error! Plot was not created !!\n')
            return
