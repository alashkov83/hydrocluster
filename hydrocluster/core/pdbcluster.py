#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import bz2
import io
import pickle
import random
import re
import time
import warnings
from collections import OrderedDict
from multiprocessing import Queue, Process
from subprocess import Popen
from urllib.error import HTTPError
from xmlrpc.client import ServerProxy

import matplotlib.cm as cm
import numpy as np
import psutil
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d

try:
    from sklearn.cluster import dbscan
    from sklearn.linear_model import LinearRegression, RANSACRegressor
    from sklearn.metrics import silhouette_score, calinski_harabaz_score
    from sklearn.metrics.pairwise import euclidean_distances
    from sklearn.neighbors import NearestNeighbors
except ImportError:
    raise ImportError

from .dbcv import DBCV

warnings.filterwarnings("ignore")

Vcoeff = 4 * np.pi / 3

n_usecpus = psutil.cpu_count()
if n_usecpus > 2:
    n_usecpus -= 1


def filterXYZandRData(Label, XYZ, Dist):
    """

    :param Label:
    :param XYZ:
    :param Dist:
    :return:
    """
    filterLabel = Label[Label != -1]
    filterXYZ = XYZ[Label != -1]
    filterR = Dist[Label != -1, :][:, Label != -1]
    return filterLabel, filterXYZ, filterR


def clusterDBSCAN(X: np.ndarray, pdist: np.ndarray, sparse_n, weight_array, eps: float, min_samples: int,
                  metric: str = 'calinski', noise_filter: bool = False) -> tuple:
    """

    :param sparse_n:
    :param noise_filter:
    :param X:
    :param pdist:
    :param weight_array:
    :param eps:
    :param min_samples:
    :param metric:
    :param: noise_filter:
    """
    core_sample_indices, labels = dbscan(sparse_n, sample_weight=weight_array,
                                         eps=eps, min_samples=min_samples, n_jobs=-1, metric='precomputed')
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
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[core_sample_indices] = True
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    if noise_filter:
        filterLabel, filterXYZ, filterR = filterXYZandRData(labels, X, pdist)
    else:
        filterLabel, filterXYZ, filterR = labels, X, pdist
    if metric == 'si_score':
        try:
            score = silhouette_score(filterR, filterLabel, metric='precomputed')
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
            score = -1
    elif metric == 'calinski':
        try:
            score = calinski_harabaz_score(filterXYZ, filterLabel)
            # The score is defined as ratio between the within-cluster dispersion and the between-cluster dispersion.
            # For more info see:
            # T.Calinski and J.Harabasz, 1974. “A dendrite method for cluster analysis”.Communications in Statistics
        except ValueError:
            score = 0
    elif metric == 'dbcv':
        # from hdbscan import validity_index as DBCV # TODO: Поэкспериментировать - работает быстрее, но странно
        # score = DBCV(filterR, filterLabel, d=3)
        try:
            score = DBCV(filterXYZ, filterLabel)  # TODO: Провести оценку с уже реализованными функциями.
            if np.isnan(score):
                raise ValueError
        except ValueError:
            score = -1
    return labels, n_clusters, core_samples_mask, score


def chunkIt(seq: list, num: int) -> list:
    """

    :param seq:
    :param num:
    :return:
    """
    out = [[] for _ in range(num)]
    seq = seq.copy()
    n = 0
    while seq:
        element = random.choice(seq)
        seq.remove(element)
        out[n % num].append(element)
        n += 1
    return out


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


def lineaRegressor(X: np.ndarray, Y: np.ndarray) -> tuple:
    """

    :param X:
    :param Y:
    :return:
    """
    modelLinear = LinearRegression()
    modelLinear.fit(X, Y)
    YfitLinear = modelLinear.predict(X)
    return YfitLinear, modelLinear.coef_[0][0], modelLinear.intercept_[0], modelLinear.score(X, Y)


def ransacRegressor(X: np.ndarray, Y: np.ndarray) -> tuple:
    """

    :param X:
    :param Y:
    :return:
    """
    modelRAMSAC = RANSACRegressor()
    modelRAMSAC.fit(X, Y)
    YfitRANSAC = modelRAMSAC.predict(X)
    return (YfitRANSAC, modelRAMSAC.estimator_.coef_[0][0], modelRAMSAC.estimator_.intercept_[0],
            modelRAMSAC.score(X[modelRAMSAC.inlier_mask_], Y[modelRAMSAC.inlier_mask_]))


def regr_cube(x: np.ndarray, y: np.ndarray, z: np.ndarray, z_correct):
    """

    :param z_correct:
    :param x:
    :param y:
    :param z:
    :return:
    """
    z = np.copy(z)
    z_uncorrect = ~z_correct
    z[z_uncorrect] = np.nan
    x = x[~np.all(z_uncorrect, axis=1)]
    z = z[~np.all(z_uncorrect, axis=1), :]
    X = np.array((Vcoeff * x ** 3), ndmin=2).T
    k = np.nanargmax(z, axis=1)
    Y = np.array(y[k], ndmin=2).T
    return X, Y


def regr_cube_alt(x: np.ndarray, y: np.ndarray, z: np.ndarray, z_correct):
    """

    :param z_correct:
    :param x:
    :param y:
    :param z:
    :return:
    """
    z = np.copy(z)
    z_uncorrect = ~z_correct
    z[z_uncorrect] = np.nan
    Y = np.array(y[~np.all(z_uncorrect, axis=0)], ndmin=2).T
    z = z[:, ~np.all(z_uncorrect, axis=0)]
    X = np.array((Vcoeff * x ** 3), ndmin=2).T
    k = np.nanargmax(z, axis=0)
    X = np.array(X[k], ndmin=2)
    return X, Y


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


class ClusterPdb:
    """

    """

    def __init__(self) -> None:
        self.metrics_name = {'calinski': 'Calinski-Harabaz score', 'si_score': 'Silhouette score', 'dbcv': 'DBCV score'}
        self.X = None
        self.pdist = None
        self.sparse_n = None
        self.labels = None
        self.noise_filter = False
        self.core_samples_mask = []
        self.n_clusters = 0
        self.score = 0
        self.min_samples = 3
        self.eps = 3.0
        self.metric = 'calinski'
        self.s_array = []
        self.htable = 'hydropathy'
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0, 'si_score')
        self.weight_array = []
        self.aa_list = []
        self.states = []
        self.figs = {'colormap': None, 'linear': None, 'ransac': None}
        self.clusterThreads = []
        self.queue = Queue()

    def clean(self) -> None:
        """

        """
        self.X = None
        self.pdist = None
        self.sparse_n = None
        self.labels = None
        self.htable = 'hydropathy'
        self.noise_filter = False
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0, 'calinski')
        self.core_samples_mask = []
        self.n_clusters = 0
        self.score = 0
        self.min_samples = 3
        self.eps = 3.0
        self.metric = 'calinski'
        self.weight_array.clear()
        self.aa_list.clear()
        self.states.clear()
        self.figs = {'colormap': None, 'linear': None, 'ransac': None}
        self.clusterThreads.clear()
        self.queue = Queue()

    def cluster(self, eps: float, min_samples: int, metric):
        """

        :param metric:
        :param eps:
        :param min_samples:
        """
        if self.X is None or self.pdist is None:
            raise ValueError
        self.min_samples = min_samples
        self.eps = eps
        self.metric = metric
        self.labels, self.n_clusters, self.core_samples_mask, self.score = clusterDBSCAN(
            self.X, self.pdist, self.sparse_n, self.weight_array, eps, min_samples, metric, self.noise_filter)

    def clusterThread(self, subParams):
        """

        :param subParams:
        """
        for eps, min_samples in subParams:
            labels, n_clusters, core_samples_mask, score = clusterDBSCAN(
                self.X, self.pdist, self.sparse_n, self.weight_array, eps, min_samples, self.metric, self.noise_filter)
            clusterResults = labels, core_samples_mask, n_clusters, score, eps, min_samples
            self.queue.put(clusterResults)
        self.queue.put(None)

    def init_cycles(self, min_eps: float, max_eps: float, step_eps: float,
                    min_min_samples: int, max_min_samples: int, n_jobs=0, metric: str = 'si_score') -> int:
        """

        :param n_jobs:
        :param min_eps:
        :param max_eps:
        :param step_eps:
        :param min_min_samples:
        :param max_min_samples:
        :param metric:
        :return:
        """
        self.metric = metric
        hyperParams = []
        for eps in np.arange(min_eps, max_eps + step_eps, step_eps):
            for min_samples in range(min_min_samples, max_min_samples + 1):
                hyperParams.append((eps, min_samples))
        n_cycles = len(hyperParams)
        if n_jobs == 0:
            n_jobs = n_usecpus
        hyperParams = chunkIt(hyperParams, n_jobs)
        self.clusterThreads.clear()
        for subParams in hyperParams:
            p = Process(target=self.clusterThread, args=(subParams,))
            p.start()
            self.clusterThreads.append(p)
        self.auto_params = min_eps, max_eps, step_eps, min_min_samples, max_min_samples, metric
        return n_cycles

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
            yield n, clusterResults[4], clusterResults[5], clusterResults[2], clusterResults[3]
            n += 1
        for p in self.clusterThreads:
            p.join()
        self.clusterThreads.clear()

    def auto(self) -> tuple:
        """

        :return:
        """
        self.states.sort(key=lambda l: l[3], reverse=True)
        state = self.states[0]
        self.labels = state[0]
        self.core_samples_mask = state[1]
        self.n_clusters = state[2]
        self.score = state[3]
        self.eps, self.min_samples = state[4], state[5]
        colormap_data = [(state[5], state[4], state[3], state[0]) for state in self.states]
        colormap_data.sort(key=lambda i: i[0])
        colormap_data.sort(key=lambda i: i[1])
        x = np.array(sorted(list({data[0] for data in colormap_data})), ndmin=1)
        y = np.array(sorted(list({data[1] for data in colormap_data})), ndmin=1)
        z = np.array([data[2] for data in colormap_data])
        z.shape = (y.size, x.size)
        if self.metric == 'si_score' or self.metric == 'dbcv':
            try:
                z_min = min([x for x in z.flat if x > -1.0])
            except ValueError:
                z_min = -1.0
        else:
            z_min = 0
        self.figs['colormap'] = (y, x, z.copy().T, z_min)
        z_correct = np.array([(True if (len(set(data[3]))) > 1 else False) for data in colormap_data], dtype=bool)
        z_correct.shape = (y.size, x.size)
        V, N, = regr_cube(y, x, z, z_correct)
        V_alt, N_alt = regr_cube_alt(y, x, z, z_correct)
        Nfit, C, B, R2 = lineaRegressor(V, N)
        Nfit_alt, C_alt, B_alt, R2_alt = lineaRegressor(V_alt, N_alt)
        if R2_alt > R2:
            V, N, Nfit, C, B, R2 = V_alt, N_alt, Nfit_alt, C_alt, B_alt, R2_alt
        self.figs['linear'] = (V, N, Nfit, C, B, R2)
        try:
            self.figs['ransac'] = ransacRegressor(V, N)
        except ValueError:
            self.figs['ransac'] = None
        return self.eps, self.min_samples

    def get_conc(self):
        """

        :return:
        """
        if self.figs is None:
            return 0.0, 0.0, 0.0, 0.0
        cl, r2l = self.figs['linear'][3], self.figs['linear'][5]
        if self.figs['ransac'] is None:
            cr, r2r = 0.0, 0.0
        else:
            cr, r2r = self.figs['ransac'][1], self.figs['ransac'][3]
        return cl, cr, r2l, r2r

    def noise_percent(self):
        """

        :return:
        """
        if self.labels is not None:
            labflat = self.labels.flatten()
            n = len(labflat)
            noise_n = len([x for x in labflat if x == -1])
            return noise_n * 100 / n
        else:
            return 0

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

    def preparser(self) -> list:
        """

        :return:
        """
        return list(sorted(set((s[21] for s in self.s_array if s[0:6] == 'ATOM  '))))

    def parser(self, selectChains: list = None, htable: str = 'hydropathy', pH: float = 7.0, res: str = '') -> tuple:
        """
        :param res:
        :param selectChains:
        :param htable:
        :param pH:
        :return:
        """
        self.clean()
        self.htable = htable
        xyz_array = []
        # www.pnas.org/cgi/doi/10.1073/pnas.1616138113 # 1.0 - -7.55 kj/mol A; residues with delta mu < 0
        nanodroplet = {'ALA': 1.269, 'VAL': 1.094, 'PRO': 1.0, 'LEU': 1.147, 'ILE': 1.289, 'PHE': 1.223, 'MET': 1.013,
                       'TRP': 1.142, 'CYS': 0.746, 'GLY': 0.605, 'THR': 0.538, 'SER': 0.472}
        # Quantitative expression of protein heterogeneity:
        # Response of amino acid side chains to their local environment.
        # Bandyopadhyay D, Mehler EL.
        # Proteins. 2008 Aug;72(2):646-59. doi: 10.1002/prot.21958. # 1.0 - 0.33; residues with values > 0
        menv = {'ALA': 1.0, 'VAL': 2.52, 'LEU': 2.64, 'ILE': 2.94, 'PHE': 2.58, 'MET': 1.64,
                'TRP': 2.03, 'CYS': 3.48, 'THR': 1.82}
        # Brylinski, M., Konieczny, L. and Roterman, I. (2006b) ‘Fuzzy-oil-drop hydrophobic force
        # field – a model to represent late-stage folding (In Silico) of lysozyme’, J. Biomol. Struct. Dyn .,
        # Vol. 23, pp.519–528.
        # Is the protein folding an aim-oriented process? Human haemoglobin as example.
        # Brylinski M, Konieczny L, Roterman I.
        # Int J Bioinform Res Appl. 2007;3(2):234-60.
        fuzzyoildrop = {'ALA': 1.0, 'VAL': 1.418, 'LEU': 1.369, 'ILE': 1.544, 'PHE': 1.583, 'MET': 1.448,
                        'TRP': 1.497, 'CYS': 1.748, 'THR': 0.538}
        # Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32. 1.0 - 1.8 residues with kdHydrophobicity > 0
        hydropathy = {'ALA': 1.0, 'VAL': 2.333, 'LEU': 2.111, 'ILE': 2.5, 'PHE': 1.556, 'MET': 1.056, 'CYS': 1.389}
        # Ikai, A.J. (1980) Thermostability and aliphatic index of globular proteins. J. Biochem. 88, 1895-1898
        aliphatic_core = {'ALA': 1.0, 'VAL': 2.9, 'LEU': 3.9, 'ILE': 3.9}
        # Reversed hydropathy table
        hydropathy_h2o = {'GLY': 1.0, 'THR': 1.75, 'TRP': 2.25, 'SER': 2., 'TYR': 3.25, 'PRO': 4., 'HIS': 8.,
                          'GLU': 8.75, 'GLN': 8.75, 'ASP': 8.75, 'ASN': 8.75, 'LYS': 9.75, 'ARG': 11.25}
        if htable == 'hydropathy':
            hydrfob = hydropathy
        elif htable == 'menv':
            hydrfob = menv
        elif htable == 'fuzzyoildrop':
            hydrfob = fuzzyoildrop
        elif htable == 'nanodroplet':
            hydrfob = nanodroplet
        elif htable == 'aliphatic_core':  # TODO: Понять биологический смысл кластеризации по этой таблице
            hydrfob = aliphatic_core
        elif htable == 'hydrophilic':  # TODO: Понять биологический смысл кластеризации по этой таблице
            hydrfob = hydropathy_h2o
        elif htable == 'positive' or htable == 'negative':
            hydrfob = calc_abs_charge(htable, pH)
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
            ' F': 19.0,
            'SE': 79.0,
            ' D': 2.0}
        if selectChains is None:
            selectChains = self.preparser()
        if res:
            list_res = list(map(lambda x: (re.findall(r'[A-Z]', x)[0], int(re.findall(r'\d{1,4}', x)[0])),
                                re.findall(r'[A-Z]\d{1,4}', res.strip().upper())))
            s_array = [s for s in self.s_array if (s[0:6] == 'ATOM  ' and (s[21], int(s[22:26])) in list_res)]
        else:
            s_array = self.s_array
        for s in (s for s in s_array if (len(s) > 21 and s[21] in selectChains)):
            if s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob) and ((current_resn is None or current_chainn is None) or (
                    current_resn == int(s[22:26]) and current_chainn == s[21])):
                current_resn = int(s[22:26])
                current_chainn = s[21]
                current_resname = s[17:20]
                xyzm = [float(s[30:38]), float(s[38:46]), float(s[46:54]), mass[s[76:78]]]
                xyzm_array = np.hstack((xyzm_array, xyzm))
            elif s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob):
                self.aa_list.append((current_resn, current_chainn, current_resname))
                self.weight_array.append(hydrfob[current_resname])
                try:
                    xyzm_array.shape = (-1, 4)
                except AttributeError:
                    raise ValueError
                xyz_array = np.concatenate((xyz_array, cmass(xyzm_array)))
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
        pdist = euclidean_distances(xyz_array)
        parse_results = len(self.aa_list), np.min(pdist[np.nonzero(
            pdist)]), np.max(pdist[np.nonzero(pdist)]), np.mean(pdist[np.nonzero(pdist)])
        sparse_n = NearestNeighbors(radius=parse_results[2], algorithm='brute', n_jobs=-1,
                                    ).fit(xyz_array).radius_neighbors_graph(xyz_array, mode='distance')
        self.X, self.pdist, self.parse_results, self.sparse_n = xyz_array, pdist, parse_results, sparse_n
        return parse_results

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
            col = np.array(col, ndmin=2)
            if k == -1:
                xyz_noise = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k])
                if xyz_noise.any():
                    ax.scatter(xyz_noise[:, 0], xyz_noise[:, 1], xyz_noise[:, 2], c='k', s=12, label='Noise')
            else:
                xyz_core = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k
                                     and self.core_samples_mask[i]])
                if xyz_core.any():
                    ax.scatter(xyz_core[:, 0], xyz_core[:, 1], xyz_core[:, 2], c=col, s=32,
                               label='Core Cluster No {:d}'.format(k + 1))
                xyz_uncore = np.array([x for i, x in enumerate(xyz_all) if self.labels[i] == k
                                       and not self.core_samples_mask[i]])
                if xyz_uncore.any():
                    ax.scatter(xyz_uncore[:, 0], xyz_uncore[:, 1], xyz_uncore[:, 2], c=col, s=12,
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
        if self.figs is not None:
            es, ms, z, z_min = self.figs['colormap']
        else:
            raise ValueError
        fig = Figure(figsize=(12, 6))
        ax1 = fig.add_subplot(121)
        ax1.set_title(self.metrics_name[self.metric])
        ax1.set_xlabel('EPS, \u212B')
        ax1.set_ylabel('MIN SAMPLES')
        ax1.grid(grid_state)
        pc1 = ax1.pcolor(es, ms, z, cmap='gnuplot', vmin=z_min)
        fig.colorbar(pc1, ax=ax1, extend='max', extendfrac=0.1)
        V, N, Nfit, C, B, R2 = self.figs['linear']
        ax11 = fig.add_subplot(122)
        ax11.set_ylabel('MIN SAMPLES')
        ax11.set_xlabel(r'$V,\ \AA^3$')
        ax11.scatter(V, N, c='k')
        ax11.scatter(Vcoeff * self.eps ** 3, self.min_samples, c='r', label='Current')
        texLINEAR = 'Linear:\n' + r'$C_h\ =\ ' + '{:.4f}'.format(C) + r'\ \AA^{-3}$' + "\n" + r'$N_0\ =\ ' + \
                    '{:.1f}$'.format(B) + '\n' + r'$R^2\ =\ ' + '{:.4f}'.format(R2) + r'$'
        ax11.plot(V, Nfit, c='r', label=texLINEAR)
        if self.figs['ransac'] is not None:
            NRfit, CR, BR, SR = self.figs['ransac']
            texRANSAC = 'RANSAC:\n' + r'$C_h\ =\ ' + '{:.4f}'.format(CR) + r'\ \AA^{-3}$' + "\n" + r'$N_0\ =\ ' + \
                        '{:.1f}$'.format(BR) + '\n' + r'$R^2\ =\ ' + '{:.4f}'.format(SR) + r'$'
            ax11.plot(V, NRfit, c='g', label=texRANSAC)
        ax11.legend(loc='best', fontsize=6)
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        return fig

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

    def open_pymol(self):
        """

        :return:
        """
        for pid in psutil.pids():
            if psutil.Process(pid).name() == 'pymol':
                cmd = psutil.Process(pid).cmdline()
                if (len(cmd) == 2 and cmd[1] == '-R') or (len(cmd) == 3 and cmd[2] == '-R'):
                    break
        else:
            Popen(args=('pymol', '-R'))
        n = 0
        while n < 10:
            try:
                pymol = ServerProxy(uri="http://localhost:9123/RPC2")
                pymol.read_pdbstr(''.join(self.s_array), 'obj')
            except ConnectionRefusedError:
                time.sleep(1)
                n += 1
            else:
                break
        else:
            raise ConnectionRefusedError
        pymol.set('label_color', 'white')
        pymol.delete('sele')
        pymol.hide('everything')
        pymol.show_as('sticks', 'all')
        dict_aa = self.get_dict_aa()
        if dict_aa is None:
            return
        colors = [cm.get_cmap('tab20b')(each) for each in np.linspace(0, 1, len(dict_aa))]
        color_names = []
        for n, colindex in enumerate(colors):
            colindex = [float(col) for n, col in enumerate(colindex) if n < 3]
            pymol.set_color('col_{:d}'.format(n), colindex)
            color_names.append("col_{:d}".format(n))
        for (k, aa_list), color in zip(dict_aa.items(), color_names):
            if aa_list:
                pymol.select('{:s}_cluster_{:d}'.format(("Core" if k[0] else "Uncore"), k[1]),
                             '{:s}'.format(
                                 "+".join(['(chain {1:s} and resi {0:d})'.format(*aac) for aac in aa_list])))
                pymol.color(color, '{:s}_cluster_{:d}'.format(("Core" if k[0] else "Uncore"), k[1]))
                pymol.show_as('spheres', '{:s}_cluster_{:d}'.format(("Core" if k[0] else "Uncore"), k[1]))
        pymol.deselect()

    def save_pymol_script(self, filename):
        """

        :param filename:
        """
        s = "from pymol import cmd\n\ncmd.set('label_color','white')\ncmd.delete ('sele')\ncmd.hide ('everything')\n" \
            "cmd.show_as('sticks', 'all')\n"
        dict_aa = self.get_dict_aa()
        colors = (cm.get_cmap('tab20b')(each) for each in np.linspace(0, 1, len(dict_aa)))
        color_names = []
        for n, colindex in enumerate(colors):
            s += "cmd.set_color('col_{:d}', [{:f}, {:f}, {:f}])\n".format(n, *colindex)
            color_names.append("col_{:d}".format(n))
        for (k, aa_list), color in zip(dict_aa.items(), color_names):
            if aa_list:
                s += "cmd.select('{:s}_cluster_{:d}', '{:s}')\n".format(
                    ("Core" if k[0] else "Uncore"), k[1], "+".join(
                        ['(chain {1:s} and resi {0:d})'.format(*aac) for aac in aa_list]))
                s += "cmd.color('{:s}', '{:s}_cluster_{:d}')\n".format(
                    color, ("Core" if k[0] else "Uncore"), k[1])
                s += "cmd.show_as('spheres', '{:s}_cluster_{:d}')\n".format(("Core" if k[0] else "Uncore"), k[1])
        s += "cmd.deselect()\n"
        with open(filename, 'wt') as f:
            f.write(s)

    def savestate(self, file: str):
        """

        :param file:
        """
        glob_state = {
            'X': self.X,
            'pdist': self.pdist,
            'sparse_n': self.sparse_n,
            'labels': self.labels,
            'noise_filter': self.noise_filter,
            'core_samples_mask': self.core_samples_mask,
            'n_clusters': self.n_clusters,
            'score': self.score,
            'eps': self.eps,
            'min_samples': self.min_samples,
            'metric': self.metric,
            'weight_array': self.weight_array,
            'aa_list': self.aa_list,
            's_array': self.s_array,
            'htable': self.htable,
            'parse_results': self.parse_results,
            'auto_params': self.auto_params,
            'states': self.states,
            'figs': self.figs}
        with bz2.open(file, 'wb') as f:
            pickle.dump(glob_state, f, protocol=4)

    def loadstate(self, file: str):
        """

        :param file:
        """
        with bz2.open(file) as f:
            global_state = pickle.load(f)
            self.clean()
        self.X = global_state['X']
        self.pdist = global_state['pdist']
        self.labels = global_state['labels']
        self.sparse_n = global_state['sparse_n']
        self.noise_filter = global_state['noise_filter']
        self.core_samples_mask = global_state['core_samples_mask']
        self.n_clusters = global_state['n_clusters']
        self.s_array = global_state['s_array']
        self.htable = global_state['htable']
        self.parse_results = global_state['parse_results']
        self.auto_params = global_state['auto_params']
        self.score = global_state['score']
        self.eps = global_state['eps']
        self.min_samples = global_state['min_samples']
        self.metric = global_state['metric']
        self.weight_array = global_state['weight_array']
        self.aa_list = global_state['aa_list']
        self.states = global_state['states']
        self.figs = global_state['figs']
