#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import math
import os.path
import sys
import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter.messagebox import askyesno, showerror, showinfo, showwarning
from tkinter.simpledialog import askstring, askfloat
from urllib.error import HTTPError

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from .. import __version__, __license__, __skversion__, newversioncheck

try:
    from ..core.pdbcluster import ClusterPdb
except ImportError:
    showerror('Error!', 'Scikit-learn not installed!')
    sys.exit(-1)


class FigureTk(tk.Toplevel):
    def __init__(self, master, title, fig, ax=None, modal=True):
        super().__init__(master=master)
        self.title(title)
        self.fig = fig
        x = master.winfo_x()
        y = master.winfo_y()
        self.menu()
        fra1 = ttk.Frame(self)
        fra1.grid(row=0, column=0)
        canvas = FigureCanvasTkAgg(fig, master=fra1)
        if ax:
            ax.mouse_init()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.wm_geometry("+{:d}+{:d}".format(x, y))
        self.resizable(False, False)
        if modal:
            self.transient(master)
            self.grab_set()
        self.focus_force()

    def save_graph(self):
        """
        """
        opt = {'parent': self,
               'filetypes': [('All supported formats', ('.eps', '.jpeg', '.jpg', '.pdf', '.pgf', '.png', '.ps',
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
                showerror('Error!', 'Failed to plot!')
            except ValueError:
                showerror('Unsupported file format!',
                          'Supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.')

    def menu(self):
        m = tk.Menu(self)
        self.config(menu=m)
        fm = tk.Menu(m)
        m.add_cascade(label='File', menu=fm)
        fm.add_command(label='Save picture', command=self.save_graph)
        fm.add_command(label='Quit', command=self.destroy)


class DialogOK(tk.Toplevel):
    def __init__(self, master, title, message):
        super().__init__(master=master)
        self.title(title)
        tk.Label(self, text=message, anchor=tk.W, justify=tk.LEFT).pack(pady=10, padx=10)
        x = master.winfo_x() + master.winfo_width() // 2
        y = master.winfo_y() + master.winfo_height() // 2
        self.wm_geometry("+{:d}+{:d}".format(x, y))
        self.resizable(False, False)
        self.transient(master)
        self.grab_set()
        self.focus_force()
        b = tk.Button(self, text="OK", command=self.ok)
        b.pack(pady=5, padx=5)

    def ok(self):
        self.destroy()


class TkGui(tk.Tk):
    """

    """

    def __init__(self) -> None:
        super().__init__()
        self.title('HydroCluster')
        self.resizable(False, False)
        self.protocol('WM_DELETE_WINDOW', self.close_win)
        self.grid = tk.BooleanVar()
        self.legend = tk.BooleanVar()
        self.grid.set(False)
        self.legend.set(False)
        self.menu()
        fra1 = tk.Frame(self)
        fra1.grid(row=0, rowspan=2, column=0)
        lab1 = tk.LabelFrame(fra1, text='Parsing', labelanchor='n', borderwidth=5)
        lab1.pack(expand=1, fill=tk.X, pady=5, padx=5)
        lab11 = tk.LabelFrame(lab1, text='Property table', labelanchor='n', borderwidth=5)
        lab11.grid(row=0, column=0, pady=5, padx=5)
        listbox_items = ['hydropathy', 'menv', 'fuzzyoildrop', 'rekkergroup', 'nanodroplet',
                         'aliphatic_core', 'hydrophilic', 'positive', 'negative', 'pgroup', 'ngroup']
        self.combox_p = ttk.Combobox(lab11, height=5, width=15, values=listbox_items, state='readonly')
        self.combox_p.pack()
        self.combox_p.set('hydropathy')
        but3 = tk.Button(lab1, text='Start', command=self.parse_pdb)
        but3.grid(row=1, column=0, pady=5)
        fra11 = tk.Frame(lab1)
        fra11.grid(row=2, column=0, pady=5, padx=5)
        l11 = tk.Label(fra11, text="No. of res.(groups):", anchor=tk.W)
        l11.grid(row=0, column=0, pady=5, padx=5, sticky="W")
        self.l11 = tk.Label(fra11, text="{0:<5d}".format(0), anchor=tk.W)
        self.l11.grid(row=0, column=1, pady=5, padx=5, sticky="W")
        l12 = tk.Label(fra11, text="Min distance (\u212B): ", anchor=tk.W)
        l12.grid(row=1, column=0, pady=5, padx=5, sticky="W")
        self.l12 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.W)
        self.l12.grid(row=1, column=1, pady=5, padx=5, sticky="W")
        l13 = tk.Label(fra11, text="Max distance (\u212B): ", anchor=tk.W)
        l13.grid(row=2, column=0, pady=5, padx=5, sticky="W")
        self.l13 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.W)
        self.l13.grid(row=2, column=1, pady=5, padx=5, sticky="W")
        l14 = tk.Label(fra11, text="Mean distance (\u212B): ", anchor=tk.W)
        l14.grid(row=3, column=0, pady=5, padx=5, sticky="W")
        self.l14 = tk.Label(fra11, text="{0:>5.3f}".format(0), anchor=tk.W)
        self.l14.grid(row=3, column=1, pady=5, padx=5, sticky="W")
        lab4 = tk.LabelFrame(fra1, text="Options", labelanchor='n', borderwidth=5)
        lab4.pack(expand=1, fill=tk.X, pady=5, padx=5)
        l41 = tk.Label(lab4, text="Metric: ", anchor=tk.W)
        l41.grid(row=0, column=0, pady=5, padx=5, sticky="W")
        listbox_items = ['calinski',
                         'si_score',
                         's_dbw',
                         ]
        self.combox = ttk.Combobox(lab4, height=5, width=10, values=listbox_items, state='readonly')
        self.combox.grid(row=0, column=1, pady=5, padx=5, sticky="W")
        self.combox.set('calinski')
        l42 = tk.Label(lab4, text="Noise filter: ", anchor=tk.W)
        l42.grid(row=1, column=0, pady=5, padx=5, sticky="W")
        listbox_items = ['comb', 'filter', 'sep', 'bind']
        self.combox_n = ttk.Combobox(lab4, height=5, width=10, values=listbox_items, state='readonly')
        self.combox_n.grid(row=1, column=1, pady=5, padx=5, sticky="W")
        self.combox_n.set('comb')
        lab2 = tk.LabelFrame(fra1, text='Auto mode', labelanchor='n', borderwidth=5)
        lab2.pack(expand=1, fill=tk.X, pady=5, padx=5)
        lab3 = tk.LabelFrame(fra1, text='Manual mode', labelanchor='n', borderwidth=5)
        lab3.pack(expand=1, fill=tk.X, pady=5, padx=5)
        lab31 = tk.LabelFrame(lab3, text='EPS (\u212B)', labelanchor='n', borderwidth=5)
        lab31.grid(row=0, column=0, pady=5, padx=5)
        self.sca1 = tk.Scale(lab31, length=210, from_=1.0, to=15.0, showvalue=1,
                             orient=tk.HORIZONTAL, resolution=0.1)
        self.sca1.set(3.0)
        self.sca1.pack()
        lab32 = tk.LabelFrame(lab3, text='MIN_SAMPLES', labelanchor='n', borderwidth=5)
        lab32.grid(row=1, column=0, pady=5, padx=5)
        self.sca2 = tk.Scale(lab32, length=210, from_=1, to=50, showvalue=1,
                             orient=tk.HORIZONTAL)
        self.sca2.set(3)
        self.sca2.pack()
        but1 = tk.Button(lab3, text='Start', command=self.run)
        but1.grid(row=2, column=0, pady=5)
        lab22 = tk.Frame(lab2)
        lab22.grid(row=1, column=0)
        l1 = tk.Label(lab22, text="Min EPS (\u212B):", anchor=tk.W)
        l1.grid(row=0, column=0, pady=5, padx=5, sticky="W")
        self.ent_min_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_min_eps.delete(0, tk.END)
        self.ent_min_eps.insert(0, '3.0')
        self.ent_min_eps.grid(row=0, column=1, pady=5, padx=5, sticky="W")
        l2 = tk.Label(lab22, text="Max EPS (\u212B):", anchor=tk.W)
        l2.grid(row=1, column=0, pady=5, padx=5, sticky="W")
        self.ent_max_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_max_eps.delete(0, tk.END)
        self.ent_max_eps.insert(0, '15.0')
        self.ent_max_eps.grid(row=1, column=1, pady=5, padx=5, sticky="W")
        l3 = tk.Label(lab22, text="Step EPS (\u212B):", anchor=tk.W)
        l3.grid(row=2, column=0, pady=5, padx=5, sticky="W")
        self.ent_step_eps = tk.Entry(lab22, width=4, bd=3)
        self.ent_step_eps.delete(0, tk.END)
        self.ent_step_eps.insert(0, '0.1')
        self.ent_step_eps.grid(row=2, column=1, pady=5, padx=5, sticky="W")
        l4 = tk.Label(lab22, text="Min MIN_SAMPLES:", anchor=tk.W)
        l4.grid(row=3, column=0, pady=5, padx=5, sticky="W")
        self.ent_min_min_samples = tk.Entry(lab22, width=4, bd=3)
        self.ent_min_min_samples.delete(0, tk.END)
        self.ent_min_min_samples.insert(0, '3')
        self.ent_min_min_samples.grid(row=3, column=1, pady=5, padx=5, sticky="W")
        l5 = tk.Label(lab22, text="Max MIN_SAMPLES:", anchor=tk.W)
        l5.grid(row=4, column=0, pady=5, padx=5, sticky="W")
        self.ent_max_min_samples = tk.Entry(lab22, width=4, bd=3)
        self.ent_max_min_samples.delete(0, tk.END)
        self.ent_max_min_samples.insert(0, '50')
        self.ent_max_min_samples.grid(row=4, column=1, pady=5, padx=5, sticky="W")
        but2 = tk.Button(lab2, text='Start', command=lambda: self.run(auto=True))
        but2.grid(row=3, column=0, pady=5)
        self.pb = ttk.Progressbar(lab2, orient='horizontal', mode='determinate', length=220)
        self.pb.grid(row=4, column=0, pady=5, padx=5)
        self.fra3 = tk.Frame(self, width=800, height=710)
        self.fra3.grid_propagate(False)
        self.fra3.grid(row=0, column=1)
        fra4 = tk.Frame(self)
        fra4.grid(row=1, column=1, pady=10)
        self.tx = tk.Text(fra4, width=100, height=9, wrap=tk.WORD)
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
        self.cls = ClusterPdb()
        self.eval('tk::PlaceWindow {:s} center'.format(self.winfo_pathname(self.winfo_id())))  # Center on screen
        self.tk.eval('::msgcat::mclocale en')  # Set the English language for standard tkinter dialog

    def _bound_to_mousewheel(self, event, tx: tk.Text):
        _ = event
        self.bind_all('<MouseWheel>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Button-4>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Button-5>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Up>', lambda e: self._on_mousewheel(e, tx))
        self.bind_all('<Down>', lambda e: self._on_mousewheel(e, tx))

    def _unbound_to_mousewheel(self, event):
        _ = event
        self.unbind_all('<MouseWheel>')
        self.unbind_all('<Button-4>')
        self.unbind_all('<Button-5>')
        self.unbind_all('<Up>')
        self.unbind_all('<Down>')

    @staticmethod
    def _on_mousewheel(event, tx: tk.Text):
        if event.num == 4 or event.keysym == 'Up':
            tx.yview_scroll(-1, 'units')
        elif event.num == 5 or event.keysym == 'Down':
            tx.yview_scroll(1, 'units')
        else:
            tx.yview_scroll(int(-1 * (event.delta / 120)), 'units')

    @staticmethod
    def about():
        """

        :return:
        """
        showinfo('About', 'Cluster analysis of hydrophobic or charged regions of macromolecules\n\n'
                          'Version: {:s}\n'
                          'Latest version on PyPi: {:s}\n'
                          'License: {:s}\n\n'
                          'Python version: {:s}\n'
                          'Platform: {:s}\n'
                          'Scikit-learn version: {:s}'.format(__version__, newversioncheck(), __license__, sys.version,
                                                              sys.platform, __skversion__))

    @staticmethod
    def readme():
        """

        :return:
        """
        import webbrowser
        webbrowser.open_new('https://alashkov83.github.io/hydrocluster/')

    def menu(self) -> None:
        """The method of initialize menu."""
        m = tk.Menu(self)  # creates a Menu object for main window
        self.config(menu=m)  # the window is configured with the indication of the menu for it
        # creates a menu item with the placement on the main menu (m)
        fm = tk.Menu(m)
        # item is located on the main menu (m)
        m.add_cascade(label='File', menu=fm)
        # a list of commands of a menu item
        fm.add_command(label='Open File', command=self.open_file)
        fm.add_command(label='Open ID PDB', command=self.open_url)
        fm.add_command(label='Load State', command=self.open_state)
        fm.add_separator()
        fm.add_command(label='Save PyMOL script', command=self.save_pymol_script)
        fm.add_command(label='Save State', command=self.save_state)
        fm.add_command(label='Save Picture', command=self.save_graph)
        fm.add_command(label='Save LOG', command=self.save_log)
        fm.add_separator()
        fm.add_command(label='Quit', command=self.close_win)
        # creates a menu item with the placement on the main menu (m)
        om = tk.Menu(m)
        # item is located on the main menu (m)
        m.add_cascade(label='Options', menu=om)
        oms = tk.Menu(m)
        om.add_cascade(label='Select clustering solution', menu=oms)
        oms.add_command(label='By local max (min)', command=lambda: self.select_sol(True))
        oms.add_command(label='By max (min) values', command=self.select_sol)
        oma = tk.Menu(m)
        # item is located on the main menu (m)
        om.add_cascade(label='Solution analysis', menu=oma)
        oma.add_command(label='Autotune colormap', command=self.colormap)
        oma.add_command(label='Autotune 3D-map', command=self.colormap3d)
        oma.add_command(label='Scan by parameter', command=self.scan_param)
        om.add_separator()
        om.add_command(label='Open PyMol', command=self.open_pymol)
        om.add_command(label='About Protein', command=self.about_protein)
        om.add_separator()
        omg = tk.Menu(m)
        # item is located on the main menu (m)
        om.add_cascade(label='Plot settings', menu=omg)
        omg.add_checkbutton(label='Plot grid', onvalue=True, offvalue=False,
                            variable=self.grid, command=self.plot_set)
        omg.add_checkbutton(label='Plot legend', onvalue=True, offvalue=False,
                            variable=self.legend, command=self.plot_set)
        om.add_command(label='Clear LOG', command=self.clean_txt)
        hm = tk.Menu(m)
        m.add_cascade(label='Help', menu=hm)
        hm.add_command(label='Readme', command=self.readme)
        hm.add_command(label='About', command=self.about)

    def close_win(self) -> None:
        """Self-destruct with the ask."""
        if askyesno('Quit', 'Are your sure?'):
            self.destroy()

    def run(self, auto: bool = False, load_state: bool = False, eps=None, min_samples=None) -> None:
        """The main algorithm of the program."""
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        if self.cls.X is None:
            return
        self.run_flag = True
        self.pb['value'] = 0
        self.pb.update()
        noise_filter = self.combox_n.get()
        metric = self.combox.get()
        if metric == 's_dbw':
            if not askyesno('Warning!', "Very slow function!\nA you sure?"):
                self.run_flag = False
                return
        if auto and not load_state:
            try:
                min_eps = float(self.ent_min_eps.get())
                if min_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = False
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Not correct value for Min EPS")
                return
            try:
                max_eps = float(self.ent_max_eps.get())
                if max_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = False
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Not correct value for Max EPS")
                return
            try:
                step_eps = float(self.ent_step_eps.get())
                if step_eps < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = False
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Not correct value for Step EPS")
                return
            try:
                min_min_samples = int(self.ent_min_min_samples.get())
                if min_min_samples < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = False
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Not correct value for Min MIN_SAMPLES")
                return
            try:
                max_min_samples = int(self.ent_max_min_samples.get())
                if max_min_samples < 0:
                    raise ValueError
            except ValueError:
                self.run_flag = False
                self.pb['value'] = 0
                self.pb.update()
                showerror("Error", "Not correct value for Max MIN_SAMPLES")
                return
            if load_state:
                self.graph()
                self.run_flag = False
                return
            if self.cls.states and self.cls.auto_params == (min_eps, max_eps,
                                                            step_eps,
                                                            min_min_samples,
                                                            max_min_samples, metric) and \
                    self.cls.noise_filter == noise_filter:
                self.run_flag = False
                return
            else:
                self.cls.noise_filter = noise_filter
                self.pb['maximum'] = self.cls.init_cycles(min_eps, max_eps, step_eps,
                                                          min_min_samples, max_min_samples, metric=metric)
                self.tx.configure(state='normal')
                self.tx.insert(tk.END, ('Starting Autoscan (range EPS: {0:.2f} - {1:.2f} \u212B,'
                                        'step EPS = {2:.2f} \u212B, range min_samples: {3:d} - {4:d}...\n').format(
                    min_eps, max_eps, step_eps, min_min_samples, max_min_samples))
                self.tx.see(tk.END)
                try:
                    for n, j, i, n_clusters, score in self.cls.auto_yield():
                        self.tx.insert(tk.END, 'Step No {0:d}: EPS = {1:.2f} \u212B, min_s = {2:d}, '
                                               'No cls = {3:d}, {4:s} = {5:.3f}\n'
                                       .format(n, j, i, n_clusters,
                                               self.cls.metrics_name[self.cls.metric], score))
                        self.tx.see(tk.END)
                        self.pb['value'] = n
                        self.pb.update()
                    eps, min_samples = self.cls.auto()
                    self.tx.insert(tk.END, 'Autoscan done... \n')
                    self.tx.configure(state='disabled')
                except ValueError:
                    showerror('Error!', 'Could not parse file or clustering failed')
                    self.tx.insert(tk.END, 'Error! Could not parse file or clustering failed\n')
                    self.tx.see(tk.END)
                    self.tx.configure(state='disabled')
                    self.run_flag = False
                    return
            showinfo('Autoscan done', ('Number of clusters = {0:d}\n{4:s} = {1:.3f}\nEPS = {2:.1f} \u212B'
                                       '\nMIN_SAMPLES = {3:d}\n').format(
                self.cls.n_clusters, self.cls.score, eps, min_samples, self.cls.metrics_name[self.cls.metric]))
            self.sca1.set(eps)
            self.sca2.set(min_samples)
        elif not auto and not load_state:
            self.cls.noise_filter = noise_filter
            if eps is None:
                eps = self.sca1.get()
            if min_samples is None:
                min_samples = self.sca2.get()
            self.pb['maximum'] = 1
            try:
                self.cls.cluster(eps, min_samples, metric)
                self.pb['value'] = 1
                self.pb.update()
            except ValueError:
                showerror('Error!', 'Could not parse file or clustering failed')
                self.run_flag = False
                return
        if self.cls.noise_percent() > 30:
            showwarning('Warning!', 'Percent of noise = {:.2f} %'.format(self.cls.noise_percent()))
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, ('Number of clusters = {0:d}\n{1:s} = {4:.3f}\nEPS = {2:.1f} \u212B\n'
                                'MIN_SAMPLES = {3:d}\nPercent of noise = {5:.2f} %{6:s}\n').format(
            self.cls.n_clusters, self.cls.metrics_name[self.cls.metric],
            self.cls.eps, self.cls.min_samples, self.cls.score,
            self.cls.noise_percent(), (' !!WARNING!!!' if self.cls.noise_percent() > 30 else '')))
        self.tx.see(tk.END)
        dict_aa = self.cls.get_dict_aa()
        if dict_aa:
            for k, aa_list in dict_aa.items():
                self.tx.configure(state='normal')
                self.tx.insert(tk.END, '\n{:s} cluster No. {:d} contains: {:s}'.format(
                    ("Core" if k[0] else "Uncore"), k[1],
                    ", ".join(['{2:s}:{1:s}{0:d}{3:s}'.format(aac[0], aac[1], aac[2],
                                                              ((':' + '-'.join(aac[3])) if (
                                                                      len(aac) > 3 and aac[3]) else ''))
                               for aac in aa_list])))
                self.tx.insert(tk.END, '\n')
        self.tx.configure(state='disabled')
        self.graph()
        self.run_flag = False

    def graph(self):
        """

        :return:
        """
        try:
            self.canvas.get_tk_widget().destroy()
            self.fig = None
        except AttributeError:
            pass
        grid, legend = self.grid.get(), self.legend.get()
        try:
            self.fig, ax = self.cls.graph(grid, legend)
        except AttributeError:
            return
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.fra3)
        ax.mouse_init()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def open_file(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        opt = {'filetypes': [('Files PDB', ('.pdb', '.PDB', '.ent')), ('Files mmCIF', ('.cif', '.CIF'))]}
        fname = askopenfilename(**opt)
        if fname:
            try:
                if os.path.splitext(fname)[1].lower() in ('.pdb', '.ent'):
                    self.cls.open_pdb(fname)
                elif os.path.splitext(fname)[1].lower() in ('.cif',):
                    self.cls.open_cif(fname)
                else:
                    showerror('Error', 'Incorrect file extension: {0:s}!'.format(os.path.splitext(fname)[1]))
                    return
            except ImportError:
                showerror('Import error', 'BioPython is not available!\nPlease install biopython and mmtf!')
                return
            except (ValueError, TypeError):
                showerror('Error', 'Incorrect CIF file: {0:s}!'.format(os.path.basename(fname)))
                return
            except FileNotFoundError:
                return
            else:
                showinfo('Infor', 'File {0:s} successfully read!'.format(os.path.basename(fname)))
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
        if url:
            try:
                self.cls.open_url(url)
            except ImportError:
                showerror('Import error', 'BioPython unavailable!\nPlease install biopython and mmtf!')
                return
            except HTTPError as e:
                showerror('Error!!', ('{1:s}\nID PDB: {0:s} not found'
                                      ' or refers to an incorrect file!').format(url, str(e)))
            else:
                showinfo('Info', 'File ID PDB: {0:s} downloaded'.format(url))
                self.parse_pdb()

    def parse_pdb(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        win = tk.Toplevel(self)
        x = self.winfo_x() + self.winfo_width() // 2
        y = self.winfo_y() + self.winfo_height() // 2
        win.wm_geometry("+{:d}+{:d}".format(x, y))
        win.title("Choice chains")
        self.resizable(False, False)
        win.transient(self)
        win.grab_set()
        win.focus_force()
        fra1 = tk.Frame(win)
        fra1.grid(row=0, column=0)
        checkVars = []
        for chain in self.cls.preparser():
            checkVar = tk.StringVar()
            CheckBox = tk.Checkbutton(fra1, text="Chain {:s}".format(chain), variable=checkVar, onvalue=chain,
                                      offvalue='', anchor=tk.NW)
            CheckBox.select()
            CheckBox.pack(expand=1, fill=tk.X, pady=5, padx=5)
            checkVars.append(checkVar)
        fra2 = tk.Frame(win)
        fra2.grid(row=1, column=0)
        butChoice = tk.Button(fra2, text='Choice', command=lambda: (self.parse_pdb_main(
            [ch for ch in (ch.get() for ch in checkVars) if ch], win)))
        butChoice.grid(row=0, column=0, pady=5)
        butAll = tk.Button(fra2, text='All', command=lambda: self.parse_pdb_main(win=win))
        butAll.grid(row=0, column=1, pady=5)
        butRes = tk.Button(fra2, text='Select residues', command=lambda: self.choise_residues(win=win))
        butRes.grid(row=1, column=0, columnspan=2, pady=5)

    def choise_residues(self, win: tk.Toplevel = None):
        """

        :param win:
        """
        if win:
            win.destroy()
        win = tk.Toplevel(self)
        win.title("Choice chains")
        self.resizable(False, False)
        win.transient(self)
        win.grab_set()
        win.focus_force()
        x = self.winfo_x() + self.winfo_width() // 2
        y = self.winfo_y() + self.winfo_height() // 2
        win.wm_geometry("+{:d}+{:d}".format(x, y))
        win.title("Choice residues")
        fra1 = tk.Frame(win)
        fra1.grid(row=0, column=0)
        tx = tk.Text(fra1, width=30, height=8, wrap=tk.WORD)
        scr = tk.Scrollbar(fra1, command=self.tx.yview)
        tx.configure(yscrollcommand=scr.set)
        tx.pack(side=tk.LEFT)
        scr.pack(side=tk.RIGHT, fill=tk.Y)
        tx.bind('<Enter>', lambda e: self._bound_to_mousewheel(e, tx))
        tx.bind('<Leave>', self._unbound_to_mousewheel)
        fra2 = tk.Frame(win)
        fra2.grid(row=1, column=0)
        butChoice = tk.Button(fra2, text='Choice',
                              command=lambda: self.parse_pdb_main(residues=tx.get('1.0', tk.END), win=win))
        butChoice.grid(row=0, column=0, pady=5)
        butCancel = tk.Button(fra2, text='Cancel', command=lambda: self.parse_pdb_main(win=win))
        butCancel.grid(row=0, column=1, pady=5)

    def parse_pdb_main(self, chains=None, win=None, residues=None):
        """

        :param residues:
        :param win:
        :param chains:
        :return:
        """
        if win:
            win.destroy()
        htable = self.combox_p.get()
        try:
            if htable in ('positive', 'negative', 'pgroup', 'ngroup'):
                pH = askfloat('Your pH', 'pH value:', initialvalue=7.0, minvalue=0.0, maxvalue=14.0)
                parse_results = self.cls.parser(htable=htable, pH=pH, selectChains=chains, res=residues)
            else:
                parse_results = self.cls.parser(htable=htable, selectChains=chains, res=residues)
        except ValueError:
            showerror('Error!', 'Invalid file format\nor file does not {:s} contain {:s}\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'menv', 'nanodroplet', 'fuzzyoildrop', 'rekkergroup')
                else 'negative' if htable in ('ngroup', 'negative') else
                'positive' if htable in ('pgroup', 'positive') else 'aliphatic',
                'groups' if htable in ('rekkergroup', 'pgroup', 'ngroup') else 'residues'))
            return
        else:
            showinfo('Info', 'File successfully parsed!\n'
                             'Property table: {0:s}\n'
                             "No. of {5:s} {6:s}: {1:d}\n"
                             "Minimum distance = {2:.3f} \u212B\n"
                             "Maximum distance = {3:.3f} \u212B\n"
                             "Mean distance = {4:.3f} \u212B\n".format(htable, parse_results[0], parse_results[1],
                                                                       parse_results[2], parse_results[3], 'hydrophobic'
                                                                       if htable in (
                    'hydropathy', 'menv', 'nanodroplet', 'fuzzyoildrop', 'rekkergroup')
                                                                       else 'negative' if htable in (
                    'ngroup', 'negative') else
                'positive' if htable in ('pgroup', 'positive') else 'aliphatic',
                                                                       'groups' if htable in ('rekkergroup', 'pgroup',
                                                                                              'ngroup') else 'residues'))
            self.l11.configure(text="{0:>5d}".format(parse_results[0]))
            self.l12.configure(text="{0:>5.3f}".format(parse_results[1]))
            self.l13.configure(text="{0:>5.3f}".format(parse_results[2]))
            self.l14.configure(text="{0:>5.3f}".format(parse_results[3]))
        self.pb['value'] = 0
        self.pb.update()
        try:
            self.canvas.get_tk_widget().destroy()
        except AttributeError:
            pass
        self.fig = None
        self.clean_txt()
        self.tx.configure(state='normal')
        self.tx.insert(tk.END, 'Property table: {0:s}\n'
                               "No. of {5:s} {6:s}: {1:d}\n"
                               "Minimum distance = {2:.3f} \u212B\n"
                               "Maximum distance = {3:.3f} \u212B\n"
                               "Mean distance = {4:.3f} \u212B\n\n".format(htable, parse_results[0], parse_results[1],
                                                                           parse_results[2], parse_results[3],
                                                                           'hydrophobic' if htable in (
                                                                               'hydropathy', 'menv', 'nanodroplet',
                                                                               'fuzzyoildrop', 'rekkergroup')
                                                                           else 'negative' if htable in (
                                                                               'ngroup', 'negative') else
                                                                           'positive' if htable in (
                                                                               'pgroup', 'positive') else 'aliphatic',
                                                                           'groups' if htable in (
                                                                               'rekkergroup', 'pgroup',
                                                                               'ngroup') else 'residues'))
        self.tx.see(tk.END)
        self.tx.configure(state='disabled')
        self.sca1.set(parse_results[1])
        self.sca2.set(round(math.log(parse_results[0], math.e)))
        self.ent_min_eps.delete(0, tk.END)
        self.ent_min_eps.insert(0, '{:.1f}'.format(parse_results[1]))

    def clean_txt(self):
        """

        """
        self.tx.configure(state='normal')
        self.tx.delete('1.0', tk.END)
        self.tx.see(tk.END)
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
        if state:
            try:
                self.cls.loadstate(state)
            except (FileNotFoundError, TypeError):
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
                min_eps, max_eps, step_eps, min_min_samples, max_min_samples, metric = self.cls.auto_params
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
                _, _, _, _, eps, min_samples = self.cls.states[0]
                self.clean_txt()
                self.sca1.set(eps)
                self.sca2.set(min_samples)
                self.combox_n.set(self.cls.noise_filter)
                self.combox.set(metric)
                self.run(auto=True, load_state=True)

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
        if state:
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
        if pmlf:
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
            showerror('Error!', 'Failed to plot!')
            return
        opt = {'parent': self,
               'filetypes': [('All supported formats', ('.eps', '.jpeg', '.jpg', '.pdf', '.pgf', '.png', '.ps',
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
                showerror('Error!', 'Failed to plot!')
            except ValueError:
                showerror('Unsupported file format!',
                          'Supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.')

    def plot_set(self):
        """

        :return:
        """
        if self.run_flag:
            return
        try:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        except AttributeError:
            pass
        self.graph()

    def colormap(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        try:
            grid = self.grid.get()
            fig = self.cls.colormap(grid)
        except ValueError:
            showinfo('Information', 'Data unavailable')
            return
        FigureTk(self, 'Colormap', fig)

    def colormap3d(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        try:
            grid = self.grid.get()
            ax, fig = self.cls.colormap3d(grid)
        except ValueError:
            showinfo('Information', 'Data unavailable')
            return
        FigureTk(self, '3DColormap', fig, ax)

    def scan_param(self):
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        if not self.cls.states:
            showerror('Error!', 'States is not available!')
            return
        min_eps, max_eps, step_eps, min_min_samples, max_min_samples, _ = self.cls.auto_params
        win = tk.Toplevel(self)
        x = self.winfo_x() + self.winfo_width() // 2
        y = self.winfo_y() + self.winfo_height() // 2
        win.wm_geometry("+{:d}+{:d}".format(x, y))
        win.title("Clustering solutions")
        win.resizable(False, False)
        win.focus_force()
        fra1 = tk.Frame(win)
        fra1.grid(row=0, column=0)
        lab31 = tk.LabelFrame(fra1, text='EPS (\u212B)', labelanchor='n', borderwidth=5)
        lab31.grid(row=0, column=0, pady=5, padx=5)
        sca1 = tk.Scale(lab31, length=400, from_=min_eps, to=max_eps, showvalue=1,
                        orient=tk.HORIZONTAL, resolution=step_eps)
        sca1.pack(fill=tk.X)
        lab32 = tk.LabelFrame(fra1, text='MIN_SAMPLES', labelanchor='n', borderwidth=5)
        lab32.grid(row=1, column=0, pady=5, padx=5)
        sca2 = tk.Scale(lab32, length=400, from_=min_min_samples, to=max_min_samples, showvalue=1, orient=tk.HORIZONTAL)
        sca2.pack(fill=tk.X)
        sca1.set(self.cls.eps)
        sca2.set(self.cls.min_samples)
        fra2 = tk.Frame(win)
        fra2.grid(row=1, column=0)
        b_eps = tk.Button(fra2, text="Scan EPS",
                          command=lambda: FigureTk(self, "Scan EPS",
                                                   self.cls.fig_scan_param('eps', int(sca2.get())), modal=False))
        b_eps.grid(row=0, column=0, pady=5, padx=5)
        b_ms = tk.Button(fra2, text="Scan MIN_SAMPLES",
                         command=lambda: FigureTk(self, "Scan MIN_SAMPLES",
                                                  self.cls.fig_scan_param('min_samples',
                                                                          float(sca1.get())), modal=False))
        b_ms.grid(row=0, column=1, pady=5, padx=5)
        b_c = tk.Button(fra2, text="Cancel", command=lambda: win.destroy())
        b_c.grid(row=0, column=2, pady=5, padx=5)

    def open_pymol(self):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        try:
            self.cls.open_pymol()
        except FileNotFoundError:
            showerror("Erros!", "PyMol not found!")

    def select_sol(self, ext=False):
        """

        :return:
        """
        if self.run_flag:
            showerror('Error!', 'The calculation is still running!')
            return
        if ext:
            sol_list = self.cls.get_nsol_ext(5)
        else:
            sol_list = self.cls.get_nsol(5)
        if not sol_list:
            return
        win = tk.Toplevel(self)
        x = self.winfo_x() + self.winfo_width() // 2
        y = self.winfo_y() + self.winfo_height() // 2
        win.wm_geometry("+{:d}+{:d}".format(x, y))
        win.title("Clustering solutions")
        self.resizable(False, False)
        win.transient(self)
        win.grab_set()
        win.focus_force()
        fra1 = tk.Frame(win)
        fra1.grid(row=0, column=0)
        for sol in sol_list:
            CheckBox = tk.Button(fra1, text="{:d}: score: {:.3f}, nclustes: {:d}, eps: {:.2f} \u212B, "
                                            "min_samples: {:d}".format(sol[0], sol[1], sol[2], sol[3], sol[4]),
                                 command=lambda x=sol[3], y=sol[4]: self.set_sol(win, x, y),
                                 anchor=tk.W)
            CheckBox.flash()
            CheckBox.pack(expand=1, fill=tk.X, pady=5, padx=5)

    def set_sol(self, win, eps, min_samples):
        if win:
            win.destroy()
        self.run(eps=eps, min_samples=min_samples)
        self.sca1.set(self.cls.eps)
        self.sca2.set(self.cls.min_samples)

    def about_protein(self):
        if not self.cls.s_array:
            showerror('Error!', 'The structure is not avaible!')
            return
        name, head, method, res, ncomp, nchain, ec, nres, mmass = self.cls.get_protein_info()
        DialogOK(self, "About protein", """Name: {:s}
HINFO: {:s}
Method: {:s}
Resolution: {:.2f} \u212B
Number of compounds: {:d}
Number of peptide chains: {:d}
EC: {:s}
Total number of residues: {:d}
Total molecular weight: {:d} Da""".format(name, head, method, res, ncomp, nchain, ec, nres, mmass))
