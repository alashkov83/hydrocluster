#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import gzip
import io
import pickle
import random
import warnings
from collections import OrderedDict
from multiprocessing import Queue, Process
from urllib.error import HTTPError

import matplotlib.cm as cm
import numpy as np
import psutil
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d

try:
    from sklearn.cluster import DBSCAN
    from sklearn.metrics import silhouette_score
    from sklearn.metrics import calinski_harabaz_score
    from sklearn.metrics.pairwise import euclidean_distances
except ImportError:
    raise ImportError

from .DBCV import DBCV

warnings.filterwarnings("ignore")


class ClusterPdb:
    """

    """

    def __init__(self) -> None:
        self.X = None
        self.pdist = None
        self.labels = None
        self.core_samples_mask = []
        self.n_clusters = 0
        self.si_score = -1
        self.calinski = 0
        self.dbcv = -1
        self.s_array = []
        self.htable = 'hydropathy'
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0)
        self.weight_array = []
        self.aa_list = []
        self.states = []
        self.clusterThreads = []
        self.queue = Queue()

    def clean(self) -> None:
        """

        """
        self.X = None
        self.pdist = None
        self.labels = None
        self.htable = 'hydropathy'
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0)
        self.core_samples_mask = []
        self.n_clusters = 0
        self.si_score = -1
        self.calinski = 0
        self.dbcv = -1
        self.weight_array.clear()
        self.aa_list.clear()
        self.states.clear()
        self.clusterThreads.clear()
        self.queue = Queue()

    @staticmethod
    def clusterDBSCAN(X: np.ndarray, pdist: np.ndarray, weight_array, eps: float, min_samples: int):
        """

        :param X:
        :param pdist:
        :param weight_array:
        :param eps:
        :param min_samples:
        """
        db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1, metric='precomputed'
                    ).fit(pdist, sample_weight=weight_array)
        # The DBSCAN algorithm considers clusters as areas of high density separated by areas of low density.
        # Due to this rather generic view, clusters found by DBSCAN can be of any shape,
        # as opposed to k-means which assumes that clusters are convex shaped.
        # The result of the DBSCAN follows the concept of core samples, namely the samples that are located in areas
        # of high density. A cluster is therefore a set of core samples,
        # each close to each other (measured by some distance measure) and a set of non-core samples that are close
        # to a core sample (but are not core samples themselves).
        # There algorithm has two parameters: min_samples and eps,
        #  which define formally what we mean when we say dense.
        # Higher min_samples or lower eps indicate higher density needed to form a cluster.
        # For more info see:
        # “A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise”
        # Ester, M., H. P. Kriegel, J. Sander, and X. Xu,
        # In Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining,
        # Portland, OR, AAAI Press, pp. 226–231. 1996
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        try:
            si_score = silhouette_score(pdist, labels, metric='precomputed')
        # The Silhouette Coefficient is calculated using the mean intra-cluster distance (a)
        # and the mean nearest-cluster distance (b) for each sample.
        # The Silhouette Coefficient for a sample is (b - a) / max(a, b).
        # To clarify, b is the distance between a sample and a nearest neighbor cluster (not containing this sample).
        # Note that Silhouette Coefficient is only defined if number of labels is 2 <= n_labels <= n_samples - 1.
        # For more info see:
        # Peter J. Rousseeuw (1987).
        # “Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis”.
        # Computational and Applied Mathematics 20: 53-65.
        except ValueError:
            si_score = -1
        try:
            calinski = calinski_harabaz_score(X, labels)
            # The score is defined as ratio between the within-cluster dispersion and the between-cluster dispersion.
            # For more info see:
            # T.Calinski and J.Harabasz, 1974. “A dendrite method for cluster analysis”.Communications in Statistics
        except ValueError:
            calinski = 0
        try:
            dbcv = DBCV(X, labels)
        except ValueError:
            dbcv = -1
        return labels, n_clusters, core_samples_mask, si_score, calinski, dbcv

    def cluster(self, eps: float, min_samples: int):
        """

        :param eps:
        :param min_samples:
        """
        if self.X is None or self.pdist is None:
            raise ValueError
        self.labels, self.n_clusters, self.core_samples_mask, self.si_score, self.calinski, self.dbcv = self.clusterDBSCAN(
            self.X, self.pdist, self.weight_array, eps, min_samples)

    def clusterThread(self, subParams):
        """

        :param subParams:
        """
        for eps, min_samples in subParams:
            labels, n_clusters, core_samples_mask, si_score, calinski, dbcv = self.clusterDBSCAN(
                self.X, self.pdist, self.weight_array, eps, min_samples)
            clusterResults = labels, core_samples_mask, n_clusters, si_score, calinski, dbcv, eps, min_samples
            self.queue.put(clusterResults)
        self.queue.put(None)

    @staticmethod
    def chunkIt(seq, num):
        out = [[] for x in range(num)]
        seq = seq.copy()
        n = 0
        while seq:
            element = random.choice(seq)
            seq.remove(element)
            out[n % num].append(element)
            n += 1
        return out

    def init_cycles(self, min_eps: float, max_eps: float, step_eps: float,
                    min_min_samples: int, max_min_samples: int, n_jobs=0) -> tuple:
        """

        :param n_jobs:
        :param min_eps:
        :param max_eps:
        :param step_eps:
        :param min_min_samples:
        :param max_min_samples:
        :return:
        """
        hyperParams = []
        for eps in np.arange(min_eps, max_eps + step_eps, step_eps):
            for min_samples in range(min_min_samples, max_min_samples + 1):
                hyperParams.append((eps, min_samples))
        if n_jobs == 0:
            if psutil.cpu_count(logical=True) == 1:
                n_jobs = 1
            else:
                n_jobs = psutil.cpu_count(logical=True) - 1
        hyperParams = self.chunkIt(hyperParams, n_jobs)
        self.clusterThreads.clear()
        for subParams in hyperParams:
            p = Process(target=self.clusterThread, args=(subParams,))
            p.start()
            self.clusterThreads.append(p)
        self.auto_params = min_eps, max_eps, step_eps, min_min_samples, max_min_samples
        return (max_min_samples - min_min_samples + 1) * np.arange(min_eps, max_eps + step_eps, step_eps).size

    def auto_yield(self) -> iter:
        """
        """
        n = 1
        self.states.clear()
        k = len(self.clusterThreads)
        while k:
            clusterResults = self.queue.get()
            if clusterResults is None:
                k -= 1
                continue
            self.states.append(clusterResults)
            yield n, clusterResults[5], clusterResults[6]
            n += 1
        for p in self.clusterThreads:
            p.join()
        self.clusterThreads.clear()

    def auto(self, metric: str = 'si_score') -> tuple:
        """

        :param metric:
        :return:
        """
        if metric == 'si_score':
            self.states.sort(key=lambda x: x[3], reverse=True)
        elif metric == 'calinski':
            self.states.sort(key=lambda x: x[4], reverse=True)
        if metric == 'dbcv':
            self.states.sort(key=lambda x: x[5], reverse=True)
        state = self.states[0]
        self.labels = state[0]
        self.core_samples_mask = state[1]
        self.n_clusters = state[2]
        self.si_score = state[3]
        self.calinski = state[4]
        self.dbcv = state[5]
        return state[6], state[7]

    def open_pdb(self, pdb: str) -> None:
        """

        :return:
        """
        try:
            with open(pdb) as f:
                self.s_array = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError

    def open_url(self, url: str) -> None:
        """

        :return:
        """
        try:
            import Bio.PDB as PDB
            from Bio.PDB.mmtf import MMTFParser
        except ImportError:
            raise ImportError
        try:
            structure = MMTFParser.get_structure_from_url(url)
        except HTTPError as e:
            raise e
        else:
            with io.StringIO() as f:
                iopdb = PDB.PDBIO()
                iopdb.set_structure(structure)
                iopdb.save(f)
                f.flush()
                f.seek(0, 0)
                self.s_array = f.readlines()

    def open_cif(self, cif_f: str):
        """

        :return:
        """
        try:
            import Bio.PDB as PDB
        except ImportError:
            raise ImportError
        parser = PDB.MMCIFParser()
        try:
            structure = parser.get_structure('X', cif_f)
        except FileNotFoundError:
            raise FileNotFoundError
        except (KeyError, ValueError, AssertionError):
            raise ValueError
        else:
            with io.StringIO() as f:
                iopdb = PDB.PDBIO()
                iopdb.set_structure(structure)
                iopdb.save(f)
                f.flush()
                f.seek(0, 0)
                self.s_array = f.readlines()

    @staticmethod
    def calc_abs_charge(res_type: str, pH: float) -> dict:
        """

        :param res_type:
        :param pH:
        :return:
        """
        pKa_dict = {'ARG': 12.5, 'ASP': 3.9, 'GLU': 4.35, 'HIS': 6.5, 'LIS': 10.35, 'TYR': 9.9, 'CYS': 8.3}
        # DEXTER S MOORE Amino Acid and Peptide Net Charges: A Simple Calculational Procedure
        # BIOCHEMICAL EDUCATION 13(1) 1985
        shrink_value = 0.1
        if res_type == 'positive':
            return {res: 1 / (1 + 10 ** (pH - pKa_dict[res])) for res in ['HIS', 'LIS', 'ARG']
                    if (1 / (1 + 10 ** (pH - pKa_dict[res]))) > shrink_value}
        elif res_type == 'negative':
            return {res: 1 / (1 + 10 ** (pKa_dict[res] - pH)) for res in ['ASP', 'GLU', 'TYR', 'CYS']
                    if (1 / (1 + 10 ** (pKa_dict[res] - pH))) > shrink_value}

    def parser(self, htable: str = 'hydropathy', pH: float = 7.0) -> tuple:
        """

        :param htable:
        :param pH:
        :return:
        """
        self.clean()
        self.htable = htable
        xyz_array = []
        # www.pnas.org/cgi/doi/10.1073/pnas.1616138113 # 1.0 - -7.55 kj/mol A;; residues with delta mu < 0
        nanodroplet = {'ALA': 1.269, 'VAL': 1.094, 'PRO': 1.0, 'LEU': 1.147, 'ILE': 1.289, 'PHE': 1.223, 'MET': 1.013,
                       'TRP': 1.142, 'CYS': 0.746, 'GLY': 0.605, 'THR': 0.538, 'SER': 0.472}
        # Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32. 1.0 - 1.8 residues with kdHydrophobicity > 0
        hydropathy = {'ALA': 1.0, 'VAL': 2.333, 'LEU': 2.111, 'ILE': 2.5, 'PHE': 1.556, 'MET': 1.056, 'CYS': 1.389}
        if htable == 'hydropathy':
            hydrfob = hydropathy
        elif htable == 'nanodroplet':
            hydrfob = nanodroplet
        elif htable == 'positive' or htable == 'negative':
            hydrfob = self.calc_abs_charge(htable, pH)
        # OLD CA-BASED PARSER
        # for s in strarr:
        #     if s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob) and s[12:16] == ' CA ':
        #         xyz = [float(s[30:38]), float(s[38:46]), float(s[46:54])]
        #         xyz_array = np.hstack((xyz_array, xyz))
        #         self.weight_array.append(hydrfob[s[17:20]])
        xyzm_array = []
        current_resn = None
        current_chainn = None
        current_resname = None
        mass = {
            ' H': 1.0,
            ' C': 12.0,
            ' N': 14.0,
            ' O': 16.0,
            ' P': 31.0,
            ' S': 32.0,
            ' F': 19.0}
        for s in self.s_array:
            if s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob) and ((current_resn is None or current_chainn is None) or (
                    current_resn == int(s[22:26]) and current_chainn == s[21])):
                current_resn = int(s[22:26])
                current_chainn = s[21]
                current_resname = s[17:20]
                xyzm = [float(s[30:38]), float(s[38:46]), float(s[46:54]), mass[s[76:78]]]
                xyzm_array = np.hstack((xyzm_array, xyzm))
            elif s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob):
                self.aa_list.append(
                    (current_resn, current_chainn, current_resname))
                self.weight_array.append(hydrfob[current_resname])
                try:
                    xyzm_array.shape = (-1, 4)
                except AttributeError:
                    raise ValueError
                xyz_array = np.hstack((xyz_array, self.cmass(xyzm_array)))
                xyzm_array = []
                current_resn = int(s[22:26])
                current_chainn = s[21]
                current_resname = s[17:20]
                xyz = [float(s[30:38]), float(s[38:46]), float(s[46:54]), mass[s[76:78]]]
                xyzm_array = np.hstack((xyzm_array, xyz))
        try:
            xyz_array.shape = (-1, 3)
        except AttributeError:
            raise ValueError
        self.X = xyz_array
        self.pdist = euclidean_distances(self.X)
        self.parse_results = len(self.aa_list), np.min(self.pdist[np.nonzero(
            self.pdist)]), np.max(self.pdist[np.nonzero(self.pdist)]), np.mean(self.pdist[np.nonzero(self.pdist)])
        return self.parse_results

    def graph(self, grid_state: bool, legend_state: bool) -> tuple:
        """

        :return:
        """
        fig = Figure(figsize=(8, 6))
        try:
            unique_labels = sorted(set(self.labels))
            xyz_all = self.X
        except AttributeError:
            raise AttributeError
        ax = axes3d.Axes3D(fig)
        colors = [cm.get_cmap('rainbow')(each) for each in np.linspace(0, 1, len(unique_labels))]
        for k, col in zip(unique_labels, colors):
            # Noise is marked in black.
            if k == -1:
                xyz_noise = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k])
                if xyz_noise.any():
                    ax.scatter(xyz_noise[:, 0], xyz_noise[:, 1], xyz_noise[:, 2], c='k', s=12, label='Noise')
            else:
                xyz_core = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k
                                     and self.core_samples_mask[i]])
                if xyz_core.any():
                    ax.scatter(xyz_core[:, 0], xyz_core[:, 1], xyz_core[:, 2], c=tuple(col), s=32,
                               label='Core Cluster No {:d}'.format(k + 1))
                xyz_uncore = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k
                                       and not self.core_samples_mask[i]])
                if xyz_uncore.any():
                    ax.scatter(xyz_uncore[:, 0], xyz_uncore[:, 1], xyz_uncore[:, 2], c=tuple(col), s=12,
                               label='Non-core Cluster No {:d}'.format(k + 1))
        fig.suptitle('Cluster analysis\n')
        ax.set_ylabel(r'$y\ \AA$')
        ax.set_xlabel(r'$x\ \AA$')
        ax.set_zlabel(r'$z\ \AA$')
        ax.grid(grid_state)
        ax.set_aspect('equal')
        max_range = np.array([self.X[:, 0].max() - self.X[:, 0].min(),
                              self.X[:, 1].max() - self.X[:, 1].min(),
                              self.X[:, 2].max() - self.X[:, 2].min()]).max() / (2.0 * 0.8)
        mid_x = (self.X[:, 0].max() + self.X[:, 0].min()) * 0.5
        mid_y = (self.X[:, 1].max() + self.X[:, 1].min()) * 0.5
        mid_z = (self.X[:, 2].max() + self.X[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        if legend_state:
            ax.legend(loc='best', frameon=False)
        return fig, ax

    def colormap(self, grid_state: bool) -> object:
        """

        :return:
        """
        if not self.states:
            raise ValueError
        colormap_data = [(state[7], state[6], state[4], state[3], state[5]) for state in self.states]
        colormap_data.sort(key=lambda i: i[0])
        colormap_data.sort(key=lambda i: i[1])
        x = np.array(sorted(list({data[0] for data in colormap_data})), ndmin=1)
        y = np.array(sorted(list({data[1] for data in colormap_data})), ndmin=1)
        z = np.array([data[2] for data in colormap_data])
        z.shape = (y.size, x.size)
        fig = Figure(figsize=(7, 8))
        ax1 = fig.add_subplot(211)
        ax1.set_title('Calinski-Harabaz_score')
        ax1.set_ylabel('EPS, \u212B')
        ax1.grid(grid_state)
        pc1 = ax1.pcolor(x, y, z, cmap='gnuplot')
        fig.colorbar(pc1, ax=ax1, extend='max', extendfrac=0.1)
        ax2 = fig.add_subplot(212)
        ax2.set_title('Silhouette score')
        ax2.set_xlabel('MIN SAMPLES')
        ax2.set_ylabel('EPS, \u212B')
        ax2.grid(grid_state)
        z = np.array([data[3] for data in colormap_data])
        try:
            z_min = min([x for x in z.flat if x > -1.0])
        except ValueError:
            z_min = -1.0
        z.shape = (y.size, x.size)
        pc2 = ax2.pcolor(x, y, z, cmap='gnuplot', vmin=z_min)
        fig.colorbar(pc2, ax=ax2, extend='max', extendfrac=0.1)
        ax3 = fig.add_subplot(232)
        ax3.set_title('DBCV score')
        ax3.set_xlabel('MIN SAMPLES')
        ax3.set_ylabel('EPS, \u212B')
        ax3.grid(grid_state)
        z = np.array([data[4] for data in colormap_data])
        try:
            z_min = min([x for x in z.flat if x > -1.0])
        except ValueError:
            z_min = -1.0
        z.shape = (y.size, x.size)
        pc3 = ax3.pcolor(x, y, z, cmap='gnuplot', vmin=z_min)
        fig.colorbar(pc3, ax=ax3, extend='max', extendfrac=0.1)
        return fig

    @staticmethod
    def cmass(str_nparray: np.ndarray) -> list:
        """Calculate the position of the center of mass."""
        mass_sum = float(str_nparray[:, 3].sum())
        mx = (str_nparray[:, 3]) * (str_nparray[:, 0])
        my = (str_nparray[:, 3]) * (str_nparray[:, 1])
        mz = (str_nparray[:, 3]) * (str_nparray[:, 2])
        c_mass_x = float(mx.sum()) / mass_sum
        c_mass_y = float(my.sum()) / mass_sum
        c_mass_z = float(mz.sum()) / mass_sum
        return [c_mass_x, c_mass_y, c_mass_z]

    def get_dict_aa(self):
        """

        :return:
        """
        if (not self.aa_list) or self.labels is None:
            return None
        aa_list = list(zip(self.aa_list, self.labels, self.core_samples_mask))
        dict_aa = OrderedDict()
        for k in sorted(set(self.labels)):
            if k != -1:
                core_aa_list = []
                uncore_aa_list = []
                for aa in aa_list:
                    if aa[1] == k and aa[2]:
                        core_aa_list.append(aa[0])
                    elif aa[1] == k and not aa[2]:
                        uncore_aa_list.append(aa[0])
                dict_aa[(1, k + 1)] = core_aa_list
                dict_aa[(0, k + 1)] = uncore_aa_list
        return dict_aa

    def save_pymol_script(self, filename):
        """

        :param filename:
        """
        s = "from pymol import cmd\n\ncmd.set('label_color','white')\ncmd.delete ('sele')\ncmd.hide ('everything')\n" \
            "cmd.show_as('sticks', 'all')\n"
        dict_aa = self.get_dict_aa()
        for k, aa_list in dict_aa.items():
            s += "cmd.select('{:s}_cluster_{:d}', '{:s}')\n".format(
                ("Core" if k[0] else "Uncore"), k[1], "+".join(
                    ['(chain {1:s} and resi {0:d})'.format(*aac) for aac in aa_list]))
            s += "cmd.color('{:s}', '{:s}_cluster_{:d}')\n".format(
                ('forest' if k[0] else 'olive'), ("Core" if k[0] else "Uncore"), k[1])
            s += "cmd.show_as('spheres', '{:s}_cluster_{:d}')\n".format(("Core" if k[0] else "Uncore"), k[1])
        with open(filename, 'wt') as f:
            f.write(s)

    def savestate(self, file: str):
        """

        :param file:
        """
        glob_state = {
            'X': self.X,
            'pdist': self.pdist,
            'labels': self.labels,
            'core_samples_mask': self.core_samples_mask,
            'n_clusters': self.n_clusters,
            'si_score': self.si_score,
            'calinski': self.calinski,
            'dbcv': self.dbcv,
            'weight_array': self.weight_array,
            'aa_list': self.aa_list,
            's_array': self.s_array,
            'htable': self.htable,
            'parse_results': self.parse_results,
            'auto_params': self.auto_params,
            'states': self.states}
        with gzip.open(file, 'wb') as f:
            pickle.dump(glob_state, f)

    def loadstate(self, file: str):
        """

        :param file:
        """
        with gzip.open(file) as f:
            global_state = pickle.load(f)
            self.clean()
        self.X = global_state['X']
        self.pdist = global_state['pdist']
        self.labels = global_state['labels']
        self.core_samples_mask = global_state['core_samples_mask']
        self.n_clusters = global_state['n_clusters']
        self.s_array = global_state['s_array']
        self.htable = global_state['htable'],
        self.parse_results = global_state['parse_results']
        self.auto_params = global_state['auto_params']
        self.si_score = global_state['si_score']
        self.calinski = global_state['calinski']
        self.dbcv = global_state['dbcv']
        self.weight_array = global_state['weight_array']
        self.aa_list = global_state['aa_list']
        self.states = global_state['states']
