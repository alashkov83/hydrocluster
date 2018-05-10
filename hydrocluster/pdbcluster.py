#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import gzip
import io
import pickle
import warnings
from urllib.error import HTTPError

import matplotlib.cm as cm
import numpy as np
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d

warnings.filterwarnings("ignore")

try:
    from sklearn.cluster import DBSCAN
    from sklearn.metrics import silhouette_score
    from sklearn.metrics import calinski_harabaz_score
    from sklearn.metrics.pairwise import euclidean_distances
except ImportError:
    raise ImportError


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
        self.s_array = []
        self.htable = 'hydropathy'
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0)
        self.weight_array = []
        self.aa_list = []
        self.states = []

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
        self.weight_array.clear()
        self.aa_list.clear()
        self.states.clear()

    def cluster(self, eps, min_samples):
        """

        :param eps:
        :param min_samples:
        """
        if self.X is None or self.pdist is None:
            raise ValueError
        db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1, metric='precomputed'
                    ).fit(self.pdist, sample_weight=self.weight_array)
        # The DBSCAN algorithm views clusters as areas of high density separated by areas of low density.
        # Due to this rather generic view, clusters found by DBSCAN can be any shape,
        # as opposed to k-means which assumes that clusters are convex shaped.
        # The central component to the DBSCAN is the concept of core samples, which are samples that are in areas
        # of high density. A cluster is therefore a set of core samples,
        # each close to each other (measured by some distance measure) and a set of non-core samples that are close
        # to a core sample (but are not themselves core samples).
        # There are two parameters to the algorithm, min_samples and eps,
        #  which define formally what we mean when we say dense.
        # Higher min_samples or lower eps indicate higher density necessary to form a cluster.
        # Cite:
        # “A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise”
        # Ester, M., H. P. Kriegel, J. Sander, and X. Xu,
        # In Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining,
        # Portland, OR, AAAI Press, pp. 226–231. 1996
        self.core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        self.core_samples_mask[db.core_sample_indices_] = True
        self.labels = db.labels_
        self.n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        try:
            self.si_score = silhouette_score(self.X, self.labels)
        # The Silhouette Coefficient is calculated using the mean intra-cluster distance (a)
        # and the mean nearest-cluster distance (b) for each sample.
        # The Silhouette Coefficient for a sample is (b - a) / max(a, b).
        # To clarify, b is the distance between a sample and the nearest cluster that the sample is not a part of.
        # Note that Silhouette Coefficient is only defined if number of labels is 2 <= n_labels <= n_samples - 1.
        # Cite:
        # Peter J. Rousseeuw (1987).
        # “Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis”.
        # Computational and Applied Mathematics 20: 53-65.
        except ValueError:
            self.si_score = -1
        try:
            self.calinski = calinski_harabaz_score(self.X, self.labels)
            # The score is defined as ratio between the within-cluster dispersion and the between-cluster dispersion.
            # Cite:
            # T.Calinski and J.Harabasz, 1974. “A dendrite method for cluster analysis”.Communications in Statistics
        except ValueError:
            self.calinski = 0

    def init_cycles(self, min_eps, max_eps, step_eps, min_min_samples, max_min_samples):
        self.auto_params = min_eps, max_eps, step_eps, min_min_samples, max_min_samples
        return (max_min_samples - min_min_samples + 1) * np.arange(min_eps, max_eps + step_eps, step_eps).size

    def auto_yield(self):
        """
        """
        min_eps, max_eps, step_eps, min_min_samples, max_min_samples = self.auto_params
        n = 1
        for j in np.arange(min_eps, max_eps + step_eps, step_eps):
            for i in range(min_min_samples, max_min_samples + 1):
                self.cluster(eps=j, min_samples=i)
                self.states.append(
                    (self.labels, self.core_samples_mask, self.n_clusters, self.si_score, self.calinski, j, i))
                yield n, j, i
                n += 1

    def auto(self, metric: str = 'si_score') -> tuple:
        """

        :param metric:
        :return:
        """
        if metric == 'si_score':
            self.states.sort(key=lambda x: x[3], reverse=True)
        elif metric == 'calinski':
            self.states.sort(key=lambda x: x[4], reverse=True)
        state = self.states[0]
        self.labels = state[0]
        self.core_samples_mask = state[1]
        self.n_clusters = state[2]
        self.si_score = state[3]
        self.calinski = state[4]
        return state[5], state[6]

    def open_pdb(self, pdb) -> None:
        """

        :return:
        """
        try:
            with open(pdb) as f:
                self.s_array = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError

    def open_url(self, url) -> None:
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

    def open_cif(self, cif_f):
        """

        :return:
        """
        try:
            import Bio.PDB as PDB
            from Bio.PDB.mmtf import MMTFParser
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

    def parser(self, htable='hydropathy'):
        self.clean()
        self.htable = htable
        xyz_array = []
        # www.pnas.org/cgi/doi/10.1073/pnas.1616138113 # 1.0 - -7.55 kj/mol A;; residues with delta mu < 0
        nanodroplet = {'ALA': 1.269, 'VAL': 1.094, 'PRO': 1.0, 'LEU': 1.147, 'ILE': 1.289, 'PHE': 1.223, 'MET': 1.013,
                       'TRP': 1.142, 'CYS': 0.746, 'GLY': 0.605, 'THR': 0.472}
        # Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32. 1.0 - 1.8 residues with kdHydrophobicity > 0
        hydropathy = {'ALA': 1.0, 'VAL': 2.333, 'LEU': 2.111, 'ILE': 2.5, 'PHE': 1.556, 'MET': 1.056, 'CYS': 1.389}
        if htable == 'hydropathy':
            hydrfob = hydropathy
        elif htable == 'nanodroplet':
            hydrfob = nanodroplet
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
                xyzm = [float(s[30:38]), float(s[38:46]),
                        float(s[46:54]), mass[s[76:78]]]
                xyzm_array = np.hstack((xyzm_array, xyzm))
            elif s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob):
                self.aa_list.append(
                    (current_resn, current_chainn, current_resname))
                self.weight_array.append(hydrfob[current_resname])
                try:
                    xyzm_array.shape = (-1, 4)
                except AttributeError:
                    raise ValueError
                xyz_array = np.hstack((xyz_array, self._cmass(xyzm_array)))
                xyzm_array = []
                current_resn = int(s[22:26])
                current_chainn = s[21]
                current_resname = s[17:20]
                xyz = [float(s[30:38]), float(s[38:46]),
                       float(s[46:54]), mass[s[76:78]]]
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

    def graph(self, grid_state, legend_state):
        """

        :return:
        """
        fig = Figure(figsize=(8, 6))
        try:
            unique_labels = set(self.labels)
        except AttributeError:
            raise AttributeError
        ax = axes3d.Axes3D(fig)
        colors = [cm.get_cmap('rainbow')(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        for k, col in zip(unique_labels, colors):
            class_member_mask = (self.labels == k)
            if k == -1:
                # Black used for noise.
                xyz = self.X[class_member_mask & ~self.core_samples_mask]
                ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c='k', s=12, label='Noise')
            else:
                xyz = self.X[class_member_mask & self.core_samples_mask]
                ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=tuple(
                    col), s=32, label='Core Cluster No {:d}'.format(k + 1))
                xyz = self.X[class_member_mask & ~self.core_samples_mask]
                ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=tuple(
                    col), s=12, label='Added Cluster No {:d}'.format(k + 1))
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

    def colormap(self, grid_state):
        """

        :return:
        """
        if not self.states:
            raise ValueError
        colormap_data = [(state[6], state[5], state[4], state[3]) for state in self.states]
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
        z.shape = (y.size, x.size)
        pc2 = ax2.pcolor(x, y, z, cmap='gnuplot')
        fig.colorbar(pc2, ax=ax2, extend='max', extendfrac=0.1)
        return fig

    @staticmethod
    def _cmass(str_nparray: np.ndarray) -> list:
        """Calculate the position of the center of mass."""
        mass_sum = float(str_nparray[:, 3].sum())
        mx = (str_nparray[:, 3]) * (str_nparray[:, 0])
        my = (str_nparray[:, 3]) * (str_nparray[:, 1])
        mz = (str_nparray[:, 3]) * (str_nparray[:, 2])
        c_mass_x = float(mx.sum()) / mass_sum
        c_mass_y = float(my.sum()) / mass_sum
        c_mass_z = float(mz.sum()) / mass_sum
        return [c_mass_x, c_mass_y, c_mass_z]

    def savestate(self, file):
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
            'weight_array': self.weight_array,
            'aa_list': self.aa_list,
            's_array': self.s_array,
            'hteble': self.htable,
            'parse_result': self.parse_results,
            'auto_params': self.auto_params,
            'states': self.states}
        with gzip.open(file, 'wb') as f:
            pickle.dump(glob_state, f)

    def loadstate(self, file):
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
        self.weight_array = global_state['weight_array']
        self.aa_list = global_state['aa_list']
        self.states = global_state['states']
