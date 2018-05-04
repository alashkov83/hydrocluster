#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import sys
import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.messagebox import askyesno
from tkinter.messagebox import showerror
from tkinter.messagebox import showinfo
from tkinter.simpledialog import askstring
from urllib.error import HTTPError

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

try:
    from .pdbcluster import ClusterPdb
except ImportError:
    showerror('Ошибка!', 'Библиотека scikit-learn не установлена!')
    sys.exit()


class TkGui(tk.Tk):
    def __init__(self, namespace) -> None:
        super().__init__()
        self.title('HydroCluster')
        self.resizable(False, False)
        self.protocol('WM_DELETE_WINDOW', self.close_win)
        self.menu()
        fra1 = tk.Frame(self)
        fra1.grid(row=0, column=0)
        lab1 = tk.LabelFrame(fra1, text='EPS (\u212B)',
                             labelanchor='n', borderwidth=5)
        lab1.grid(row=0, column=0, pady=5, padx=5)
        self.sca1 = tk.Scale(lab1, length=300, from_=1.0, to=15.0,
                             showvalue=1, orient=tk.HORIZONTAL, resolution=0.1)
        self.sca1.pack()
        lab2 = tk.LabelFrame(fra1, text='MIN_SAMPLES',
                             labelanchor='n', borderwidth=5)
        lab2.grid(row=0, column=1, pady=5, padx=5)
        self.sca2 = tk.Scale(lab2, length=300, from_=1,
                             to=50, showvalue=1, orient=tk.HORIZONTAL)
        self.sca2.pack()
        but1 = tk.Button(fra1, text='Старт!',
                         command=lambda: self.run(auto=False))
        but1.grid(row=0, column=2, padx=10)
        fra2 = tk.Frame(self)
        fra2.grid(row=1, column=0)
        but2 = tk.Button(fra2, text='Авто',
                         command=lambda: self.run(auto=True))
        but2.grid(row=0, column=1, padx=10)
        lab3 = tk.LabelFrame(fra2, text='Метрика автоподбора',
                             labelanchor='n', borderwidth=5)
        lab3.grid(row=0, column=0, pady=5, padx=5)
        listbox_items = ['calinski', 'si_score']
        self.combox = ttk.Combobox(
            lab3, height=5, width=15, values=listbox_items)
        self.combox.pack()
        self.combox.set('calinski')
        lab4 = tk.LabelFrame(fra2, text='Progress: ',
                             labelanchor='n', borderwidth=5)
        lab4.grid(row=0, column=2, pady=5, padx=5)
        self.pb = ttk.Progressbar(
            lab4, orient='horizontal', mode='determinate', length=450)
        self.pb.pack()
        self.fra3 = tk.Frame(self, width=800, height=650)
        self.fra3.grid(row=2, column=0)
        self.fra3.grid_propagate(False)
        fra4 = tk.Frame(self)
        fra4.grid(row=3, column=0, pady=10)
        self.tx = tk.Text(fra4, width=100, height=10)
        scr = tk.Scrollbar(fra4, command=self.tx.yview)
        self.tx.configure(yscrollcommand=scr.set, state='disabled')
        self.tx.pack(side=tk.LEFT)
        scr.pack(side=tk.RIGHT, fill=tk.Y)
        self.tx.bind(
            '<Enter>', lambda e: self._bound_to_mousewheel(e, self.tx))
        self.tx.bind('<Leave>', self._unbound_to_mousewheel)
        self.run_flag = False
        self.s_array = None
        self.run_flag = False
        self.fig = None
        self.canvas = None
        self.toolbar = None
        self.grid = False
        self.legend = False
        self.namespace = namespace
        self.cls = ClusterPdb()

    def _bound_to_mousewheel(self, event, tx):
        self.bind_all('<MouseWheel>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Button-4>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Button-5>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Up>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Down>', lambda e: self._on_mousewheel(e, tx))

    def _unbound_to_mousewheel(self, event):
        self.unbind_all('<MouseWheel>')
        self.unbind_all('<Button-4>')
        self.unbind_all('<Button-5>')
        self.unbind_all('<Up>')
        self.unbind_all('<Down>')

    @staticmethod
    def _on_mousewheel(event, tx):
        if event.num == 4 or event.keysym == 'Up':
            tx.yview_scroll(-1, 'units')
        elif event.num == 5 or event.keysym == 'Down':
            tx.yview_scroll(1, 'units')
        else:
            tx.yview_scroll(int(-1 * (event.delta / 120)), 'units')

    @staticmethod
    def about():
        showinfo('Информация', 'Кластерный анализ гидрофобных областей макромолекул')

    def menu(self) -> None:
        """Метод инициалиции меню."""
        m = tk.Menu(self)  # создается объект Меню на главном окне
        self.config(menu=m)  # окно конфигурируется с указанием меню для него
        # создается пункт меню с размещением на основном меню (m)
        fm = tk.Menu(m)
        # пункту располагается на основном меню (m)
        m.add_cascade(label='Файл', menu=fm)
        # формируется список команд пункта меню
        fm.add_command(label='Открыть PDB', command=self.open_pdb)
        fm.add_command(label='Открыть CIF', command=self.open_cif)
        fm.add_command(label='Открыть ID PDB', command=self.open_url)
        fm.add_command(label='Открыть состояние', command=self.open_state)
        fm.add_command(label='Сохранить состояние', command=self.save_state)
        fm.add_command(label='Сохранить рисунок', command=self.save_graph)
        fm.add_command(label='Сохранить LOG', command=self.save_log)
        fm.add_command(label='Выход', command=self.close_win)
        # создается пункт меню с размещением на основном меню (m)
        om = tk.Menu(m)
        # пункту располагается на основном меню (m)
        m.add_cascade(label='Опции', menu=om)
        om.add_command(label='Сетка графика', command=self.grid_set)
        om.add_command(label='Легенда', command=self.legend_set)
        om.add_command(label='Состав гидрофобных ядер', command=self.resi)
        om.add_command(label='Цветовая карта автонастройки', command=self.colormap)
        om.add_command(label='Очистить LOG', command=self.clean_txt)
        m.add_command(label='Справка', command=self.about)

    def close_win(self) -> None:
        """Самоуничтожение с вопросом."""
        if askyesno('Выход', 'Вы точно хотите выйти?'):
            self.destroy()

    def run(self, auto: bool = False) -> None:
        """Основной алгоритм программы."""
        if self.run_flag:
            showerror('Ошибка!', 'Расчет уже идёт!')
            return
        self.run_flag = True
        self.pb['value'] = 0
        self.pb.update()
        if auto:
            metric = self.combox.get()
            if self.cls.states:
                eps, min_samples = self.cls.auto(metric=metric)
            else:
                min_eps = 3.0
                max_eps = 15.0
                step_eps = 0.1
                min_min_samples = 2
                max_min_samples = 50
                nstep_eps = round((max_eps - min_eps) / step_eps)
                self.pb['maximum'] = (max_min_samples - min_min_samples + 1) * nstep_eps
                self.tx.configure(state='normal')
                self.tx.insert(tk.END, ('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                                        'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
                    min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
                try:
                    for n, j, i in self.cls.auto_yield(min_eps=min_eps, max_eps=max_eps, nstep_eps=nstep_eps,
                                                       min_min_samples=min_min_samples,
                                                       max_min_samples=max_min_samples):
                        self.tx.insert(tk.END, ('Step No {0:d}: EPS = {1:.2f} \u212B, '
                                                'min_samples = {2:d}, No clusters = {3:d}, '
                                                'Silhouette score = {4:.3f} '
                                                'Calinski score = {5:.3f}\n').format(
                            n, j, i, self.cls.n_clusters, self.cls.si_score, self.cls.calinski))
                        self.pb['value'] = n
                        self.pb.update()
                    eps, min_samples = self.cls.auto(metric=metric)
                    self.tx.insert(tk.END, 'Autoscan done... \n')
                    self.tx.configure(state='disabled')
                except ValueError:
                    showerror('Ошибка!', 'Не загружен файл\nили ошибка кластерного анализа!')
                    self.tx.insert(tk.END, 'Ошибка! Не загружен файл или ошибка кластерного анализа!\n')
                    self.tx.configure(state='disabled')
                    self.run_flag = False
                    return
            self.sca1.set(eps)
            self.sca2.set(min_samples)
        else:
            eps = self.sca1.get()
            min_samples = self.sca2.get()
            self.pb['maximum'] = 1
            try:
                self.cls.cluster(eps, min_samples)
                self.pb['value'] = 1
                self.pb.update()
            except ValueError:
                showerror('Ошибка!', 'Не загружен файл\nили ошибка кластерного анализа!')
                self.run_flag = False
                return
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, ('Estimated number of clusters: {0:d}\nSilhouette Coefficient: {1:.3f}\n'
                                'Calinski and Harabaz score: {4:.3f}\nEPS: {2:.1f} \u212B\n'
                                'MIN_SAMPLES: {3:d}\n').format(self.cls.n_clusters,
                                                               self.cls.si_score, eps, min_samples, self.cls.calinski))
        self.tx.configure(state='disabled')
        self.graph()
        self.run_flag = False

    def graph(self):
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        grid, legend = self.grid, self.legend
        try:
            self.fig, ax = self.cls.graph(grid, legend)
        except AttributeError:
            return
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.fra3)
        ax.mouse_init()
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.fra3)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(fill=tk.BOTH, side=tk.TOP, expand=1)

    def open_pdb(self):
        if self.run_flag:
            showerror('Ошибка!', 'Расчет уже идёт!')
            return
        opt = {'filetypes': [('Файлы PDB', ('.pdb', '.PDB', '.ent')), ('Все файлы', '.*')]}
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
            showerror('Ошибка', 'Неверный формат\nлибо файл не содержит гидрофоьных остатков!')
            return
        else:
            showinfo('Информация', 'Файл распарсен!')
        self.pb['value'] = 0
        self.pb.update()
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.fig = None
        self.clean_txt()

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
