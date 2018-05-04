#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import os
import sys
from urllib.error import HTTPError

import matplotlib

matplotlib.use('Agg')
try:
    from .pdbcluster import ClusterPdb
except ImportError:
    print('Ошибка! Библиотека scikit-learn не установлена!')
    sys.exit()
try:
    import progressbar2 as progressbar
except ImportError:
    import progressbar


class Cli():
    def __init__(self, namespace) -> None:
        self.namespace = namespace
        if self.namespace.output:
            self.newdir = self.namespace.output
        else:
            self.newdir = 'output'
        try:
            os.makedirs(self.newdir, exist_ok=True)
        except OSError:
            print('Невозможно создать каталог ' + self.newdir)
            sys.exit(-1)
        self.log_name = os.path.join(self.newdir, '{:s}'.format(self.output))
        self.fig = None
        self.cls = ClusterPdb()
        self.open_file()
        self.parse_pdb()
        self.run()
        self.graph()

    def log_append(self, line):
        print(line)
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
            self.graph()

    def graph(self):
        grid, legend = True, True
        try:
            self.fig, ax = self.cls.graph(grid, legend)
        except AttributeError:
            self.log_append()

    def open_file(self):
        if not self.namespace.input:
            print('File not defined')
            sys.exit(-1)
        if self.namespace.input.split('.')
            pdb = askopenfilename(**opt)
        if pdb:
            try:
                self.cls.open_pdb(pdb)
            except FileNotFoundError:
                return
            else:
                showinfo('Информация', 'Файл прочитан!')
                self.parse_pdb()
        else:
            return

    def open_url(self):
        url = askstring('Загрузить', 'ID PDB:')
        if url is not None:
            try:
                self.cls.open_url(url)
            except ImportError:
                showerror('Ошибка импорта',
                          'BioPython недоступен!'
                          '\nДля исправления установите biopython и mmtf!')
                return
            except HTTPError as e:
                showerror('Ошибка!', ('{1:s}\nID PDB: {0:s} не найден'
                                      ' или ссылается на некорректный файл!').format(url, str(e)))
            else:
                showinfo('Информация', 'Файл загружен')
                self.parse_pdb()

    def open_cif(self):
        opt = {'filetypes': [('Файлы mmCIF', ('.cif', '.CIF')), ('Все файлы', '.*')]}
        cif_f = askopenfilename(**opt)
        if cif_f:
            try:
                self.cls.open_cif(cif_f)
            except ImportError:
                showerror('Ошибка импорта',
                          'BioPython недоступен!'
                          '\nДля исправления установите biopython и mmtf!')
                return
            except FileNotFoundError:
                return
            except ValueError:
                showerror('Ошибка!', 'Некорректный CIF файл: {0:s}!'.format(cif_f))
                return
            else:
                showinfo('Информация', 'Файл прочитан!')
                self.parse_pdb()

    def parse_pdb(self):
        try:
            self.cls.parser()
        except ValueError:
            self.log_append('Ошибка! Неверный формат\nлибо файл не содержит гидрофоьных остатков!\n')
            return

    def clean_txt(self):
        self.tx.configure(state='normal')
        self.tx.delete('1.0', tk.END)
        self.tx.configure(state='disabled')

    def open_state(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет уже идёт!')
            return
        opt = {'filetypes': [('Файл данных', ('.dat', '.DAT')), ('Все файлы', '.*')],
               'title': 'Загрузить состояние'}
        state = askopenfilename(**opt)
        try:
            self.cls.loadstate(state)
        except FileNotFoundError:
            return
        except (ValueError, OSError):
            showerror("Ошибка!", "Файл неверного формата!")
            return
        else:
            self.run(auto=True)

    def save_state(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет уже идёт!')
            return
        opt = {'filetypes': [('Файл данных', ('.dat', '.DAT')), ('Все файлы', '.*')],
               'initialfile': 'myfile.dat',
               'title': 'Сохранить состояние'}
        state = asksaveasfilename(**opt)
        try:
            self.cls.savestate(state)
        except FileNotFoundError:
            return

    def save_log(self):
        opt = {'parent': self, 'filetypes': [('LOG', '.log'), ],
               'initialfile': 'myfile.log',
               'title': 'Сохранить LOG'}
        sa = asksaveasfilename(**opt)
        if sa:
            letter = self.tx.get(1.0, tk.END)
            try:
                with open(sa, 'w') as f:
                    f.write(letter)
            except FileNotFoundError:
                pass

    def save_graph(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет не закончен!')
            return
        if self.fig is None:
            showerror('Ошибка!', 'График недоступен!')
            return
        opt = {'parent': self,
               'filetypes': [('Все поддерживаесые форматы', ('.eps', '.jpeg', '.jpg', '.pdf', '.pgf', '.png', '.ps',
                                                             '.raw', '.rgba', '.svg', '.svgz', '.tif', '.tiff')), ],
               'initialfile': 'myfile.png',
               'title': 'Сохранить график'}
        sa = asksaveasfilename(**opt)
        if sa:
            try:
                self.fig.savefig(sa, dpi=600)
            except FileNotFoundError:
                return
            except AttributeError:
                showerror('Ошибка!', 'График недоступен!')
            except ValueError:
                showerror('Неподдерживаемый формат файла рисунка!',
                          'Поддреживаемые форматы: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.')

    def grid_set(self):
        self.grid = bool(askyesno('Cетка', 'Отобразить?'))
        if self.run_flag:
            return
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.graph()

    def legend_set(self):
        self.legend = bool(askyesno('Легенда', 'Отобразить?'))
        if self.run_flag:
            return
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.graph()

    def resi(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет уже идёт!')
            return
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
                self.tx.configure(state='normal')
                self.tx.insert(tk.END, '\nIn Core cluster No {:d} included: '.format(k + 1))
                self.tx.configure(state='disabled')
                for aac in core_aa_list:
                    self.tx.configure(state='normal')
                    self.tx.insert(tk.END, '{2:s}:{1:s}{0:d} '.format(*aac))
                    self.tx.configure(state='disabled')
                if uncore_aa_list:
                    self.tx.configure(state='normal')
                    self.tx.insert(tk.END, '\nIn UNcore cluster No {:d} included: '.format(k + 1))
                    self.tx.configure(state='disabled')
                    for aac in uncore_aa_list:
                        self.tx.configure(state='normal')
                        self.tx.insert(tk.END, '{2:s}:{1:s}{0:d} '.format(*aac))
                        self.tx.configure(state='disabled')
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, '\n\n')
        self.tx.configure(state='disabled')

    def colormap(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет не закончен!')
            return
        try:
            grid = self.grid
            fig = self.cls.colormap(grid)
        except ValueError:
            showinfo('Информация', 'Данные недоступны')
            return
        win_cls = tk.Toplevel(self)
        win_cls.title("ColorMaps")
        win_cls.minsize(width=600, height=600)
        win_cls.resizable(False, False)
        fra4 = ttk.Frame(win_cls)
        fra4.grid(row=0, column=0)
        canvas = FigureCanvasTkAgg(fig, master=fra4)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar = NavigationToolbar2TkAgg(canvas, fra4)
        toolbar.update()
        canvas._tkcanvas.pack(fill=tk.BOTH, side=tk.TOP, expand=1)
