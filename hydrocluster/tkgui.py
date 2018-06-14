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
from tkinter.simpledialog import askstring, askfloat
from urllib.error import HTTPError

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

try:
    from .pdbcluster import ClusterPdb
except ImportError:
    showerror('Error!', 'Scikit-learn is not installed!')
    sys.exit()


class TkGui(tk.Tk):
    """

    """
    def __init__(self, namespace) -> None:
        super().__init__()
        self.title('HydroCluster')
        self.resizable(False, False)
        self.protocol('WM_DELETE_WINDOW', self.close_win)
        self.menu()
        fra1 = tk.Frame(self)
        fra1.grid(row=0, rowspan=2, column=0)
        lab1 = tk.LabelFrame(fra1, text='Parsing', labelanchor='n', borderwidth=5)
        lab1.grid(row=0, column=0, pady=5)
        lab11 = tk.LabelFrame(lab1, text='Property table', labelanchor='n', borderwidth=5)
        lab11.grid(row=0, column=0, pady=5, padx=5)
        listbox_items = ['hydropathy', 'nanodroplet', 'positive', 'negative']
        self.combox_p = ttk.Combobox(lab11, height=5, width=15, values=listbox_items)
        self.combox_p.pack()
        self.combox_p.set('hydropathy')
        but3 = tk.Button(lab1, text='Start', command=self.parse_pdb)
        but3.grid(row=1, column=0, pady=5)
        fra11 = tk.Frame(lab1)
        fra11.grid(row=2, column=0, pady=5, padx=5)
        l11 = tk.Label(fra11, text="No of residues:", anchor=tk.NW)
        l11.grid(row=0, column=0, pady=5, padx=5)
        self.l11 = tk.Label(fra11, text="{0:>5d}".format(0), anchor=tk.NW)
        self.l11.grid(row=0, column=1, pady=5, padx=5)
        l12 = tk.Label(fra11, text="Min distance (\u212B): ", anchor=tk.NW)
        l12.grid(row=1, column=0, pady=5, padx=5)
        self.l12 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.NW)
        self.l12.grid(row=1, column=1, pady=5, padx=5)
        l13 = tk.Label(fra11, text="Max distance (\u212B): ", anchor=tk.NW)
        l13.grid(row=2, column=0, pady=5, padx=5)
        self.l13 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.NW)
        self.l13.grid(row=2, column=1, pady=5, padx=5)
        l14 = tk.Label(fra11, text="Mean distance (\u212B): ", anchor=tk.NW)
        l14.grid(row=3, column=0, pady=5, padx=5)
        self.l14 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.NW)
        self.l14.grid(row=3, column=1, pady=5, padx=5)
        lab3 = tk.LabelFrame(fra1, text='Manual mode', labelanchor='n', borderwidth=5)
        lab3.grid(row=2, column=0, pady=5)
        lab31 = tk.LabelFrame(lab3, text='EPS (\u212B)', labelanchor='n', borderwidth=5)
        lab31.grid(row=0, column=0, pady=5, padx=5)
        self.sca1 = tk.Scale(lab31, length=200, from_=1.0, to=15.0, showvalue=1, orient=tk.HORIZONTAL, resolution=0.1)
        self.sca1.pack()
        lab32 = tk.LabelFrame(lab3, text='MIN_SAMPLES', labelanchor='n', borderwidth=5)
        lab32.grid(row=1, column=0, pady=5, padx=5)
        self.sca2 = tk.Scale(lab32, length=200, from_=1, to=50, showvalue=1, orient=tk.HORIZONTAL)
        self.sca2.pack()
        but1 = tk.Button(lab3, text='Start', command=lambda: self.run(auto=False))
        but1.grid(row=2, column=0, pady=5)
        lab2 = tk.LabelFrame(fra1, text='Auto mode', labelanchor='n', borderwidth=5)
        lab2.grid(row=1, column=0, pady=5)
        lab22 = tk.Frame(lab2)
        lab22.grid(row=1, column=0)
        l1 = tk.Label(lab22, text="Min EPS (\u212B):", anchor=tk.NW)
        l1.grid(row=0, column=0, pady=5, padx=5)
        self.ent_min_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_min_eps.delete(0, tk.END)
        self.ent_min_eps.insert(0, '3.0')
        self.ent_min_eps.grid(row=0, column=1, pady=5, padx=5)
        l2 = tk.Label(lab22, text="Max EPS (\u212B):", anchor=tk.NW)
        l2.grid(row=1, column=0, pady=5, padx=5)
        self.ent_max_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_max_eps.delete(0, tk.END)
        self.ent_max_eps.insert(0, '15.0')
        self.ent_max_eps.grid(row=1, column=1, pady=5, padx=5)
        l3 = tk.Label(lab22, text="Step EPS (\u212B):", anchor=tk.NW)
        l3.grid(row=2, column=0, pady=5, padx=5)
        self.ent_step_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_step_eps.delete(0, tk.END)
        self.ent_step_eps.insert(0, '0.1')
        self.ent_step_eps.grid(row=2, column=1, pady=5, padx=5)
        l4 = tk.Label(lab22, text="Min MIN_SAMPLES:", anchor=tk.NW)
        l4.grid(row=3, column=0, pady=5, padx=5)
        self.ent_min_min_samples = tk.Entry(lab22, width=4, bd=3)
        self.ent_min_min_samples.delete(0, tk.END)
        self.ent_min_min_samples.insert(0, '2')
        self.ent_min_min_samples.grid(row=3, column=1, pady=5, padx=5)
        l5 = tk.Label(lab22, text="Max MIN_SAMPLES:", anchor=tk.NW)
        l5.grid(row=4, column=0, pady=5, padx=5)
        self.ent_max_min_samples = tk.Entry(lab22, width=4, bd=3)
        self.ent_max_min_samples.delete(0, tk.END)
        self.ent_max_min_samples.insert(0, '50')
        self.ent_max_min_samples.grid(row=4, column=1, pady=5, padx=5)
        but2 = tk.Button(lab2, text='Start', command=lambda: self.run(auto=True))
        but2.grid(row=3, column=0, pady=5)
        lab21 = tk.LabelFrame(lab2, text='Metric shrink', labelanchor='n', borderwidth=5)
        lab21.grid(row=0, column=0, pady=5, padx=5)
        listbox_items = ['calinski', 'si_score']
        self.combox = ttk.Combobox(lab21, height=5, width=15, values=listbox_items)
        self.combox.pack()
        self.combox.set('calinski')
        lab23 = tk.LabelFrame(lab2, text='Progress: ', labelanchor='n', borderwidth=5)
        lab23.grid(row=4, column=0, pady=5, padx=5)
        self.pb = ttk.Progressbar(lab23, orient='horizontal', mode='determinate', length=200)
        self.pb.pack()
        self.fra3 = tk.Frame(self, width=800, height=720)
        self.fra3.grid_propagate(False)
        self.fra3.grid(row=0, column=1)
        fra4 = tk.Frame(self)
        fra4.grid(row=1, column=1, pady=10)
        self.tx = tk.Text(fra4, width=100, height=8, wrap=tk.WORD)
        scr = tk.Scrollbar(fra4, command=self.tx.yview)
        self.tx.configure(yscrollcommand=scr.set, state='disabled')
        self.tx.pack(side=tk.LEFT)
        scr.pack(side=tk.RIGHT, fill=tk.Y)
        self.tx.bind('<Enter>', lambda e: self._bound_to_mousewheel(e, self.tx))
        self.tx.bind('<Leave>', self._unbound_to_mousewheel)
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
        """

        """
        showinfo('About', 'Cluster analysis of hydrophobic or charged regions of macromolecules')

    def menu(self) -> None:
        """The method of initialize menu."""
        m = tk.Menu(self)  # creates a Menu object for main window
        self.config(menu=m)  # the window is configured with the indication of the menu for it
        # creates a menu item with the placement on the main menu (m)
        fm = tk.Menu(m)
        # item is located on the main menu (m)
        m.add_cascade(label='File', menu=fm)
        # a list of commands of a menu item
        fm.add_command(label='Open PDB', command=self.open_pdb)
        fm.add_command(label='Open CIF', command=self.open_cif)
        fm.add_command(label='Open ID PDB', command=self.open_url)
        fm.add_command(label='Load state', command=self.open_state)
        fm.add_command(label='Save PyMOL script', command=self.save_pymol_script)
        fm.add_command(label='Save state', command=self.save_state)
        fm.add_command(label='Save picture', command=self.save_graph)
        fm.add_command(label='Save LOG', command=self.save_log)
        fm.add_command(label='Quit', command=self.close_win)
        # creates a menu item with the placement on the main menu (m)
        om = tk.Menu(m)
        # item is located on the main menu (m)
        m.add_cascade(label='Options', menu=om)
        om.add_command(label='Plot grid', command=self.grid_set)
        om.add_command(label='Plot legend', command=self.legend_set)
        om.add_command(label='Cores content', command=self.resi)
        om.add_command(label='Autotune colormap', command=self.colormap)
        om.add_command(label='Clear LOG', command=self.clean_txt)
        m.add_command(label='About', command=self.about)

    def close_win(self) -> None:
        """Self-destruct with the ask."""
        if askyesno('Quit', 'Are your sure?'):
            self.destroy()

    def run(self, auto: bool = False) -> None:
        """The main algorithm of the program."""
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        self.run_flag = True
        self.pb['value'] = 0
        self.pb.update()
        if auto:
            metric = self.combox.get()
            try:
                min_eps = float(self.ent_min_eps.get())
                if min_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = True
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Non correct value for Min EPS")
                return
            try:
                max_eps = float(self.ent_max_eps.get())
                if max_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = True
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Non correct value for Max EPS")
                return
            try:
                step_eps = float(self.ent_step_eps.get())
                if step_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = True
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Non correct value for Step EPS")
                return
            try:
                min_min_samples = int(self.ent_min_min_samples.get())
                if min_min_samples < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = True
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Non correct value for Min MIN_SAMPLES")
                return
            try:
                max_min_samples = int(self.ent_max_min_samples.get())
                if max_min_samples < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = True
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Non correct value for Max MIN_SAMPLES")
                return
            if self.cls.states and self.cls.auto_params == (min_eps, max_eps, step_eps,
                                                            min_min_samples, max_min_samples):
                eps, min_samples = self.cls.auto(metric=metric)
            else:
                self.pb['maximum'] = self.cls.init_cycles(min_eps, max_eps, step_eps, min_min_samples, max_min_samples)
                self.tx.configure(state='normal')
                self.tx.insert(tk.END, ('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                                        'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
                    min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
                try:
                    for n, j, i in self.cls.auto_yield():
                        self.tx.insert(tk.END, 'Step No {0:d}: EPS = {1:.2f} \u212B, min_s = {2:d}, '
                                               'No cls = {3:d}, Si sc = {4:.3f}, Calinski_sc = {5:.3f}\n'
                                       .format(n, j, i, self.cls.n_clusters,
                                               self.cls.si_score, self.cls.calinski))
                        self.pb['value'] = n
                        self.pb.update()
                    eps, min_samples = self.cls.auto(metric=metric)
                    self.tx.insert(tk.END, 'Autoscan done... \n')
                    self.tx.configure(state='disabled')
                except ValueError:
                    showerror('Error!', 'File was not parse or clustering was fail')
                    self.tx.insert(tk.END, 'Error! File was not parse or clustering was fail\n')
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
                showerror('Error!', 'File was not parse or clustering was fail')
                self.run_flag = False
                return
        showinfo('Autoscan done', ('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                                   'Calinski-Harabaz score = {4:.3f}\nEPS = {2:.1f} \u212B\nMIN_SAMPLES = {3:d}\n'
                                   ).format(self.cls.n_clusters, self.cls.si_score, eps, min_samples,
                                            self.cls.calinski))
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, ('Number of clusters = {0:d}\nSilhouette Coefficient = {1:.3f}\n'
                                'Calinski-Harabaz score = {4:.3f}\nEPS = {2:.1f} \u212B\nMIN_SAMPLES = {3:d}\n'
                                ).format(self.cls.n_clusters, self.cls.si_score, eps, min_samples, self.cls.calinski))
        self.tx.configure(state='disabled')
        self.graph()
        self.run_flag = False

    def graph(self):
        """

        :return:
        """
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
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('Files PDB', ('.pdb', '.PDB', '.ent')), ('All files', '.*')]}
        pdb = askopenfilename(**opt)
        if pdb:
            try:
                self.cls.open_pdb(pdb)
            except FileNotFoundError:
                return
            else:
                showinfo('Info', 'File was read!')
                self.parse_pdb()
        else:
            return

    def open_url(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        url = askstring('Download', 'ID PDB:')
        if url is not None:
            try:
                self.cls.open_url(url)
            except ImportError:
                showerror('Import error', 'Bio Python unavailable!\nPlease install biopython and mmtf!')
                return
            except HTTPError as e:
                showerror('Error!!', ('{1:s}\nID PDB: {0:s} not found'
                                      ' or refers to an incorrect file!').format(url, str(e)))
            else:
                showinfo('Info', 'File ID PDB: {0:s} was downloaded'.format(url))
                self.parse_pdb()

    def open_cif(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('Files mmCIF', ('.cif', '.CIF')), ('Все файлы', '.*')]}
        cif_f = askopenfilename(**opt)
        if cif_f:
            try:
                self.cls.open_cif(cif_f)
            except ImportError:
                showerror('Import error', 'Bio Python is not available!\nPlease install biopython and mmtf!')
                return
            except FileNotFoundError:
                return
            except ValueError:
                showerror('Error', 'Incorrect CIF file: {0:s}!'.format(cif_f))
                return
            else:
                showinfo('Infor', 'File {0:s} was read!'.format(cif_f))
                self.parse_pdb()

    def parse_pdb(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        htable = self.combox_p.get()
        try:
            if htable == 'positive' or htable == 'negative':
                pH = askfloat('Your pH', 'pH value:', initialvalue=7.0, minvalue=0.0, maxvalue=14.0)
                parse_results = self.cls.parser(htable=htable, pH=pH)
            else:
                parse_results = self.cls.parser(htable=htable)
        except ValueError:
            showerror('Error!', 'Invalid file format\nor file does not {:s} contain resides\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'nanodroplet')
                else 'negative' if htable == 'negative' else 'positive'))
            return
        else:
            showinfo('Info', 'File was parsed!\nPTable: {:s}\n'.format(htable) +
                     "No of residues: {:d}\nMinimum distance = {:.3f} \u212B\n"
                     "Maximum distance = {:.3f} \u212B\nMean distance = {:.3f} \u212B\n".format(*parse_results))
            self.l11.configure(text="{0:>5d}".format(parse_results[0]))
            self.l12.configure(text="{0:>5.3f}".format(parse_results[1]))
            self.l13.configure(text="{0:>5.3f}".format(parse_results[2]))
            self.l14.configure(text="{0:>5.3f}".format(parse_results[3]))
        self.pb['value'] = 0
        self.pb.update()
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.fig = None
        self.clean_txt()
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, 'PTable: {:s}\n'.format(htable) +
                       "No of residues: {:d}\nMinimum distance = {:.3f} \u212B\n"
                       "Maximum distance = {:.3f} \u212B\nMean distance = {:.3f} \u212B\n\n".format(*parse_results))
        self.tx.configure(state='disabled')

    def clean_txt(self):
        """

        """
        self.tx.configure(state='normal')
        self.tx.delete('1.0', tk.END)
        self.tx.configure(state='disabled')

    def open_state(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('Data file', ('.dat', '.DAT')), ('All files', '.*')], 'title': 'Load state'}
        state = askopenfilename(**opt)
        try:
            self.cls.loadstate(state)
        except FileNotFoundError:
            return
        except (ValueError, OSError):
            showerror("Error!", "Invalid file format!")
            return
        else:
            htable = self.cls.htable
            self.combox_p.set(htable)
            nor, mind, maxd, meand = self.cls.parse_results
            self.l11.configure(text="{0:>5d}".format(nor))
            self.l12.configure(text="{0:>5.3f}".format(mind))
            self.l13.configure(text="{0:>5.3f}".format(maxd))
            self.l14.configure(text="{0:>5.3f}".format(meand))
            min_eps, max_eps, step_eps, min_min_samples, max_min_samples = self.cls.auto_params
            self.ent_min_eps.delete(0, tk.END)
            self.ent_min_eps.insert(0, '{:.1f}'.format(min_eps))
            self.ent_max_eps.delete(0, tk.END)
            self.ent_max_eps.insert(0, '{:.1f}'.format(max_eps))
            self.ent_step_eps.delete(0, tk.END)
            self.ent_step_eps.insert(0, '{:.1f}'.format(step_eps))
            self.ent_min_min_samples.delete(0, tk.END)
            self.ent_min_min_samples.insert(0, '{:d}'.format(min_min_samples))
            self.ent_max_min_samples.delete(0, tk.END)
            self.ent_max_min_samples.insert(0, '{:d}'.format(max_min_samples))
            self.run(auto=True)

    def save_state(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('Data file', ('.dat', '.DAT')), ('All files', '.*')], 'initialfile': 'myfile.dat',
               'title': 'Save state'}
        state = asksaveasfilename(**opt)
        try:
            self.cls.savestate(state)
        except FileNotFoundError:
            return

    def save_pymol_script(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('PyMOL script', ('.py',)), ('All files', '.*')], 'initialfile': 'myfile.py',
               'title': 'Save PyMOL script'}
        pmlf = asksaveasfilename(**opt)
        try:
            self.cls.save_pymol_script(pmlf)
        except FileNotFoundError:
            return
        except AttributeError:
            showerror('Error!', 'Script unavailable!')

    def save_log(self):
        """

        """
        opt = {'parent': self, 'filetypes': [('LOG', '.log'), ], 'initialfile': 'myfile.log', 'title': 'Save LOG'}
        sa = asksaveasfilename(**opt)
        if sa:
            letter = self.tx.get(1.0, tk.END)
            try:
                with open(sa, 'w', encoding='utf-8') as f:
                    f.write(letter)
            except FileNotFoundError:
                pass

    def save_graph(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        if self.fig is None:
            showerror('Error!', 'Plot unavailable!')
            return
        opt = {'parent': self,
               'filetypes': [('All support formats', ('.eps', '.jpeg', '.jpg', '.pdf', '.pgf', '.png', '.ps',
                                                      '.raw', '.rgba', '.svg', '.svgz', '.tif', '.tiff')), ],
               'initialfile': 'myfile.png',
               'title': 'Save plot'}
        sa = asksaveasfilename(**opt)
        if sa:
            try:
                self.fig.savefig(sa, dpi=600)
            except FileNotFoundError:
                return
            except AttributeError:
                showerror('Error!', 'Plot unavailable!')
            except ValueError:
                showerror('Unsupported file format!',
                          'Supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.')

    def grid_set(self):
        """

        :return:
        """
        self.grid = bool(askyesno('Plot grid', 'Display?'))
        if self.run_flag:
            return
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.graph()

    def legend_set(self):
        """

        :return:
        """
        self.legend = bool(askyesno('Plot legend', 'Display?'))
        if self.run_flag:
            return
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.graph()

    def resi(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        dict_aa = self.cls.get_dict_aa()
        if not dict_aa:
            return
        for k, aa_list in dict_aa.items():
            self.tx.configure(state='normal')
            self.tx.insert(tk.END, '\nIn {:s} cluster No {:d} included: {:s}'.format(
                ("Core" if k[0] else "Uncore"), k[1], ", ".join(['{2:s}:{1:s}{0:d}'.format(*aac) for aac in aa_list])))
            self.tx.configure(state='disabled')
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, '\n\n')
        self.tx.configure(state='disabled')

    def colormap(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        try:
            grid = self.grid
            fig = self.cls.colormap(grid)
        except ValueError:
            showinfo('Information', 'Data unavailable')
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
