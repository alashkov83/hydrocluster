#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import os
import shutil
import sys
from urllib.error import HTTPError

import matplotlib

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg
try:
    from .pdbcluster import ClusterPdb
except ImportError:
    print('Ошибка! Библиотека scikit-learn не установлена!')
    sys.exit()
try:
    import progressbar2 as progressbar
except ImportError:
    import progressbar


class Cli:
    def __init__(self, namespace) -> None:
        self.namespace = namespace
        if self.namespace.output:
            self.newdir = self.namespace.output
        else:
            self.newdir = 'output'
        if os.path.exists(self.newdir):
            shutil.rmtree(self.newdir)
        try:
            os.makedirs(self.newdir, exist_ok=True)
        except OSError:
            print('Невозможно создать каталог ' + self.newdir)
            sys.exit(-1)
        self.log_name = os.path.join(self.newdir, '{:s}'.format(self.namespace.output + '.log'))
        self.cls = ClusterPdb()
        self.open_file()
        self.parse_pdb()
        self.run()
        self.resi()
        self.save_state()
        self.graph()
        self.colormap()


    def log_append(self, line):
        print(line, end='')
        with open(self.log_name, 'at') as f:
            f.write(line)

    def run(self) -> None:
        """Основной алгоритм программы."""

        metric = self.namespace.score
        min_eps = self.namespace.emin
        max_eps = self.namespace.emax
        step_eps = self.namespace.estep
        min_min_samples = self.namespace.smin
        max_min_samples = self.namespace.smax
        nstep_eps = round((max_eps - min_eps) / step_eps)
        self.log_append(('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                         'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
            min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
        bar1 = progressbar.ProgressBar(maxval=(max_min_samples - min_min_samples + 1) * nstep_eps,
                                       redirect_stdout=True).start()
        try:
            for n, j, i in self.cls.auto_yield(min_eps=min_eps, max_eps=max_eps, nstep_eps=nstep_eps,
                                               min_min_samples=min_min_samples,
                                               max_min_samples=max_min_samples):
                self.log_append(('Step No {0:d}: EPS = {1:.2f} \u212B, '
                                 'min_samples = {2:d}, No clusters = {3:d}, '
                                 'Silhouette score = {4:.3f} Calinski score = {5:.3f}\n').format(
                    n, j, i, self.cls.n_clusters, self.cls.si_score, self.cls.calinski))

                bar1.update(n)
            eps, min_samples = self.cls.auto(metric=metric)
            self.log_append('Autoscan done... \n')
        except ValueError:
            self.log_append('Не загружен файл\nили ошибка кластерного анализа!\n')
            sys.exit(-1)
        else:
            self.log_append(('Estimated number of clusters: {0:d}\nSilhouette Coefficient: {1:.3f}\n'
                             'Calinski and Harabaz score: {4:.3f}\nEPS: {2:.1f} \u212B\n'
                             'MIN_SAMPLES: {3:d}\n').format(self.cls.n_clusters,
                                                            self.cls.si_score, eps, min_samples, self.cls.calinski))

    def graph(self):
        grid, legend = True, True
        sa = os.path.join(self.newdir, '{:s}'.format(self.namespace.output + '.png'))
        try:
            fig, ax = self.cls.graph(grid, legend)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Ошибка график не построен!\n')

    def open_file(self):
        if not self.namespace.input:
            self.log_append('Filename not defined')
            sys.exit(-1)
        if self.namespace.input.split('.')[-1].strip().lower() == 'pdb':
            pdb_f = self.namespace.input
            try:
                self.cls.open_pdb(pdb_f)
            except FileNotFoundError:
                self.log_append('File {:s} not found!\n'.format(pdb_f))
                sys.exit(-1)
            else:
                self.log_append('Файл {:s} прочитан!\n'.format(pdb_f))
                self.parse_pdb()
        elif self.namespace.input.split('.')[-1].strip().lower() == 'cif':
            cif_f = self.namespace.input
            try:
                self.cls.open_cif(cif_f)
            except ImportError:
                self.log_append('Ошибка импорта! BioPython недоступен!\nДля исправления установите biopython и mmtf!\n')
                sys.exit(-1)
            except FileNotFoundError:
                sys.exit(-1)
            except ValueError:
                self.log_append('Ошибка! Некорректный CIF файл: {0:s}!\n'.format(cif_f))
                sys.exit(-1)
            else:
                self.log_append('Файл {:s} прочитан!\n')
                self.parse_pdb()
        else:
            url = self.namespace.input
            try:
                self.cls.open_url(url)
            except ImportError:
                self.log_append('Ошибка импорта! BioPython недоступен!\nДля исправления установите biopython и mmtf!')
                return
            except HTTPError as e:
                self.log_append('Ошибка! {1:s}\nID PDB: {0:s} не найден или ссылается на некорректный файл!').format(
                    url, str(e))
            else:
                self.log_append('Файл загружен\n')
                self.parse_pdb()

    def parse_pdb(self):
        try:
            self.cls.parser()
        except ValueError:
            self.log_append('Ошибка! Неверный формат\nлибо файл не содержит гидрофоьных остатков!\n')
            return

    def save_state(self):
        st = os.path.join(self.newdir, '{:s}'.format(self.namespace.output + '.bin'))
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

    def colormap(self):
        sa = os.path.join(self.newdir, '{:s}'.format(self.namespace.output + '.cm.png'))
        try:
            fig = self.cls.colormap(grid_state=True)
            canvas = FigureCanvasAgg(fig)
            canvas.print_png(sa)
        except AttributeError:
            self.log_append('Ошибка график не построен!\n')
