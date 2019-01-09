#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 04.05.18"""

import bz2
import gc
import io
import pickle
import random
import re
import time
import warnings
from collections import OrderedDict
from multiprocessing import Queue, Process, Lock
from subprocess import Popen
from urllib.error import HTTPError
from xmlrpc.client import ServerProxy

import matplotlib.cm as cm
import numpy as np
import psutil
import scipy.ndimage.filters as filters
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d
from scipy.spatial.distance import cdist

try:
    from sklearn.cluster import dbscan
    from sklearn.linear_model import LinearRegression, RANSACRegressor
    from sklearn.metrics import silhouette_score
    from sklearn.metrics.pairwise import euclidean_distances
except ImportError:
    raise ImportError

try:
    from sklearn.metrics import calinski_harabaz_score
except ImportError:
    from sklearn.metrics import calinski_harabasz_score as calinski_harabaz_score  # Compat for 0.23 and upper

from .. import __version__
from .s_dbw import S_Dbw

warnings.filterwarnings("ignore")

Vcoeff = 4 * np.pi / 3

n_usecpus = psutil.cpu_count()
if n_usecpus > 2:
    n_usecpus -= 1


def detect_local_extrema(arr: np.ndarray, x: np.ndarray, y: np.ndarray, ext: str = 'max', bg=0, n: int = 5) -> list:
    """

    :param arr:
    :param x:
    :param y:
    :param ext:
    :param bg:
    :param n:
    :return:
    """
    epsilon = 0.00000001
    arr = arr.copy()
    shape = arr.shape
    rand_matrix = np.random.uniform(-epsilon, epsilon, size=shape)
    rand_matrix[arr == bg] = 0
    arr_noise = arr + rand_matrix
    i = 1
    j = 1
    while j < n:
        neighborhood = np.ones((shape[0] // i, shape[1] // i))
        if ext == 'min':
            local_ext_val = filters.minimum_filter(arr_noise, footprint=neighborhood)
        else:
            local_ext_val = filters.maximum_filter(arr_noise, footprint=neighborhood)
        local_ext = (local_ext_val == arr_noise) & (local_ext_val != bg)
        j = len(np.nonzero(local_ext)[0])
        i += 1
    xind, yind = np.nonzero(local_ext)[0], np.nonzero(local_ext)[1]
    xext, yext = x[xind], y[yind]
    val_ext = arr[xind, yind]
    sols = list(zip(val_ext, xext, yext, xind, yind))
    if ext == 'min':
        sols.sort(key=lambda a: a[0])
    else:
        sols.sort(key=lambda a: a[0], reverse=True)
    return sols


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


def bind_noise_lab(X, labels):
    """

    :param X:
    :param labels:
    :return:
    """
    labels = labels.copy()
    if -1 not in set(labels):
        return labels
    if len(set(labels)) == 1 and -1 in set(labels):
        return labels
    label_id = []
    label_new = []
    for i in range(len(labels)):
        if labels[i] == -1:
            point = np.array([X[i]])
            dist = cdist(X[labels != -1], point)
            lid = np.where(np.all(X == X[labels != -1][np.argmin(dist), :], axis=1))[0][0]
            label_id.append(i)
            label_new.append(labels[lid])
    labels[np.array(label_id)] = np.array(label_new)
    return labels


def sep_noise_lab(labels):
    """

    :param labels:
    :return:
    """
    labels = labels.copy()
    max_label = np.max(labels)
    j = max_label + 1
    for i in range(len(labels)):
        if labels[i] == -1:
            labels[i] = j
            j += 1
    return labels


def clusterDBSCAN(X: np.ndarray, pdist: np.ndarray, weight_array, eps: float, min_samples: int,
                  metric: str = 'calinski', noise_filter: str = 'comb') -> tuple:
    """

    :param noise_filter:
    :param X:
    :param pdist:
    :param weight_array:
    :param eps:
    :param min_samples:
    :param metric:
    :param: noise_filter:
    """
    # TODO: Separate clustering and evaluation functions
    core_sample_indices, labels = dbscan(pdist, sample_weight=weight_array,
                                         eps=eps, min_samples=min_samples,
                                         algorithm='brute', n_jobs=-1, metric='precomputed')
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
    # For precomputed metric use only 'brute' algorythm
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[core_sample_indices] = True
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    if noise_filter == 'filter':
        filterLabel, filterXYZ, filterR = filterXYZandRData(labels, X, pdist)
    elif noise_filter == 'sep':
        filterLabel = sep_noise_lab(labels)
        filterXYZ, filterR = X, pdist
    elif noise_filter == 'bind':
        filterLabel = bind_noise_lab(X, labels)
        filterXYZ, filterR = X, pdist
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
    elif metric == 's_dbw':
        # M. Halkidi and M. Vazirgiannis, “Clustering validity assess-
        # ment: Finding the optimal partitioning of a data set,” in
        # ICDM, Washington, DC, USA, 2001, pp. 187–194.
        try:
            score = S_Dbw(filterXYZ, filterLabel)
        except ValueError:
            score = np.inf
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
    pKa_dict = {'ARG': 12.48, 'ASP': 3.86, 'GLU': 4.35, 'HIS': 6.04, 'LYS': 10.35, 'TYR': 9.04, 'CYS': 8.3}
    # DEXTER S MOORE Amino Acid and Peptide Net Charges: A Simple Calculational Procedure
    # BIOCHEMICAL EDUCATION 13(1) 1985
    shrink_value = 0.1
    if res_type == 'positive':
        return {res: 1 / (1 + 10 ** (pH - pKa_dict[res])) for res in ['HIS', 'LYS', 'ARG']
                if (1 / (1 + 10 ** (pH - pKa_dict[res]))) > shrink_value}
    elif res_type == 'negative':
        return {res: 1 / (1 + 10 ** (pKa_dict[res] - pH)) for res in ['ASP', 'GLU', 'TYR', 'CYS']
                if (1 / (1 + 10 ** (pKa_dict[res] - pH))) > shrink_value}


def calc_group_charge(table_type: str, pH: float) -> dict:
    """

    :param table_type:
    :param pH:
    :return:
    """
    group_dict = {'ARG': (('NE',), ('CZ',), ('NH1',), ('NH2',)),
                  'ASP': (('CG',), ('OD1',), ('OD2',)),
                  'CYS': (('SG',),),
                  'GLU': (('CD',), ('OE1',), ('OE2',)),
                  'HIS': (('CG',), ('CD2',), ('CE1',), ('ND1',), ('NE2',)),
                  'LYS': (('NZ',),),
                  'TYR': (('OH',),),
                  }
    dict_aa = {}
    if table_type == 'ngroup':
        dict_aa = calc_abs_charge('negative', pH)
    elif table_type == 'pgroup':
        dict_aa = calc_abs_charge('positive', pH)
    return {key: ((group_dict[key] + (value,)),) for (key, value) in dict_aa.items()}


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


def regr_cube(x: np.ndarray, y: np.ndarray, z: np.ndarray, z_correct, rev: bool = False):
    """

    :param rev:
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
    if rev:
        k = np.nanargmax(z, axis=1)
    else:
        k = np.nanargmax(z, axis=1)
    Y = np.array(y[k], ndmin=2).T
    return X, Y


def regr_cube_alt(x: np.ndarray, y: np.ndarray, z: np.ndarray, z_correct, rev: bool = False):
    """

    :param rev:
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
    if rev:
        k = np.nanargmin(z, axis=0)
    else:
        k = np.nanargmax(z, axis=0)
    X = np.array(X[k], ndmin=2)
    return X, Y


def cmass(str_nparray: np.ndarray) -> list:
    """Calculate the position of the center of mass."""
    center = np.average(str_nparray[:, 0:3], axis=0, weights=str_nparray[:, 3])
    return center


def create_group(group_table: dict, res_name: str, xyzm_array: np.ndarray, atom_list: list) -> tuple:
    """

    :param group_table:
    :param res_name:
    :param xyzm_array:
    :param atom_list:
    :return:
    """
    if res_name in group_table:
        real_groups = []
        weights = []
        groups_xyz_array = []
        for group in group_table[res_name]:
            real_atoms = []
            xyzm_group = []
            for atoms in group[:-1]:
                for alt_atom in atoms:
                    for index_atom, real_atom in enumerate(atom_list):
                        if real_atom == alt_atom:
                            real_atoms.append(real_atom)
                            xyzm_group = np.hstack((xyzm_group, xyzm_array[index_atom]))
            if real_atoms:
                xyzm_group.shape = (-1, 4)
                groups_xyz_array = np.concatenate((groups_xyz_array, cmass(xyzm_group)))
                weights.append(group[-1])
                real_groups.append(tuple(real_atoms))
        return groups_xyz_array, weights, real_groups


def draw_scan_param(x, y, y1, y2, y3, y4, htable, metric, xparametr, const_str):
    """

    :param x:
    :param y:
    :param y1:
    :param y2:
    :param y3:
    :param y4:
    :param htable:
    :param metric:
    :param xparametr:
    :param const_str:
    :return:
    """
    fig = Figure(figsize=(8, 8))
    ax1 = fig.add_subplot(311)
    ax1.set_title(metric + ' vs ' + xparametr + '\nhtable: ' + htable + ", " + const_str)
    ax1.set_xlabel(xparametr)
    ax1.set_ylabel(metric)
    ax1.grid(True)
    if xparametr == 'EPS (\u212B)':
        ax1.plot(x, y)
    elif xparametr == 'MIN_SAMPLES':
        ax1.bar(x, y)
    ax2 = fig.add_subplot(312)
    ax2.set_title('No. of clusters vs ' + xparametr)
    ax2.set_xlabel(xparametr)
    ax2.set_ylabel('No. of clusters')
    ax2.grid(True)
    if xparametr == 'EPS (\u212B)':
        ax2.plot(x, y1)
    elif xparametr == 'MIN_SAMPLES':
        ax2.bar(x, y1)
    ax3 = fig.add_subplot(313)
    ax3.set_title('%  points in clusters vs ' + xparametr)
    ax3.set_xlabel(xparametr)
    ax3.set_ylabel('%  points in clusters')
    ax3.grid(True)
    if xparametr == 'EPS (\u212B)':
        ax3.plot(x, y2, color='b', label='all')
        ax3.plot(x, y3, color='r', label='core')
        ax3.plot(x, y4, color='g', label='uncore')
    elif xparametr == 'MIN_SAMPLES':
        ax3.bar(x, y2, color='b', label='all')
        ax3.bar(x, y3, color='r', label='core')
        ax3.bar(x, y4, color='g', label='uncore')
    ax3.grid(True)
    ax3.legend(loc='best')
    fig.tight_layout(pad=0.4, h_pad=0.5)
    return fig


def notnoise_percent(labels: np.ndarray) -> float:
    """

    :param labels:
    :return:
    """
    labflat = labels.copy().flatten()
    n = len(labflat)
    not_noise_n = len([x for x in labflat if x != -1])
    return not_noise_n * 100 / n


def notnoise_percent_core(labels: np.ndarray, core_mask: np.ndarray) -> float:
    """

    :param labels:
    :param core_mask:
    :return:
    """
    labflat = labels.copy().flatten()
    coreflat = core_mask.copy().flatten()
    n = len(labflat)
    not_noise_core = len([x for x in zip(labflat, coreflat) if x[0] != -1 and x[1]])
    return not_noise_core * 100 / n


def calculate_scan(states: list, param, xparm: str) -> tuple:
    """

    :param states:
    :param param:
    :param xparm:
    :return:
    """
    epsilon = 0.00000001
    if xparm == 'eps':
        x = [state[5] for state in states if abs(state[4] - param) <= epsilon]
        y = [state[3] for state in states if abs(state[4] - param) <= epsilon]
        y1 = [state[2] for state in states if abs(state[4] - param) <= epsilon]
        y2 = [notnoise_percent(state[0]) for state in states if abs(state[4] - param) <= epsilon]
        y3 = [notnoise_percent_core(state[0], state[1]) for state in states if abs(state[4] - param) <= epsilon]
    elif xparm == 'min_samples':
        x = [state[4] for state in states if state[5] == param]
        y = [state[3] for state in states if state[5] == param]
        y1 = [state[2] for state in states if state[5] == param]
        y2 = [notnoise_percent(state[0]) for state in states if state[5] == param]
        y3 = [notnoise_percent_core(state[0], state[1]) for state in states if state[5] == param]
    y4 = list(map(lambda a: a[0] - a[1], zip(y2, y3)))
    x, y, y1, y2, y3, y4 = list(zip(*sorted((zip(x, y, y1, y2, y3, y4)), key=lambda tup: tup[0])))
    return x, y, y1, y2, y3, y4


class ClusterPdb:
    """

    """

    def __init__(self) -> None:
        self.metrics_name = {'calinski': 'Calinski-Harabasz score',
                             'si_score': 'Silhouette score',
                             's_dbw'   : 'S_Dbw'
                             }
        self.X = None
        self.pdist = None
        self.labels = None
        self.noise_filter = 'comb'
        self.core_samples_mask = []
        self.n_clusters = 0
        self.score = 0
        self.min_samples = 3
        self.eps = 3.0
        self.metric = 'calinski'
        self.s_array = []
        self.htable = 'hydropathy'
        self.parse_results = (0, 0.0, 0.0, 0.0)
        self.auto_params = (0.0, 0.0, 0.0, 0, 0, 'calinski')
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
        self.labels = None
        self.htable = 'hydropathy'
        self.noise_filter = ''
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
            self.X, self.pdist, self.weight_array, eps=eps, min_samples=min_samples,
            metric=metric, noise_filter=self.noise_filter)

    def clusterThread_old(self, subParams: list):  # Old process-worker
        """

        :param subParams:
        """
        for eps, min_samples in subParams:
            labels, n_clusters, core_samples_mask, score = clusterDBSCAN(
                self.X, self.pdist, self.weight_array,
                eps=eps, min_samples=min_samples, metric=self.metric, noise_filter=self.noise_filter)
            clusterResults = labels, core_samples_mask, n_clusters, score, eps, min_samples
            self.queue.put(clusterResults)
        self.queue.put(None)
        gc.collect(2)

    def clusterThread(self, queue_in: Queue, lock_in: Lock):
        """

        :param queue_in:
        :param lock_in
        """
        while lock_in.acquire():
            if queue_in.empty():  # Because of multithreading/multiprocessing semantics, this is not reliable.
                lock_in.release()  # Additional synchronised primitive was added for safety
                break
            eps, min_samples = queue_in.get()
            lock_in.release()
            labels, n_clusters, core_samples_mask, score = clusterDBSCAN(
                self.X, self.pdist, self.weight_array,
                eps=eps, min_samples=min_samples, metric=self.metric, noise_filter=self.noise_filter)
            clusterResults = labels, core_samples_mask, n_clusters, score, eps, min_samples
            self.queue.put(clusterResults)
        self.queue.put(None)
        gc.collect(2)

    def init_cycles_old(self, min_eps: float, max_eps: float, step_eps: float,  # Old process manager
                        min_min_samples: int, max_min_samples: int, n_jobs=0, metric: str = 'calinski') -> int:
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
        # TODO: can find results in self.states first?
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
            p = Process(target=self.clusterThread_old, args=(subParams,))
            p.start()
            self.clusterThreads.append(p)
        self.auto_params = min_eps, max_eps, step_eps, min_min_samples, max_min_samples, metric
        return n_cycles

    def init_cycles(self, min_eps: float, max_eps: float, step_eps: float,
                    min_min_samples: int, max_min_samples: int, n_jobs=0, metric: str = 'calinski') -> int:
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
        # TODO: can find results in self.states first?
        self.metric = metric
        hyperParams = []
        for eps in np.arange(min_eps, max_eps + step_eps, step_eps):
            for min_samples in range(min_min_samples, max_min_samples + 1):
                hyperParams.append((eps, min_samples))
        n_cycles = len(hyperParams)
        if n_jobs == 0:
            n_jobs = n_usecpus
        queue_in = Queue()
        lock_in = Lock()
        random.shuffle(hyperParams)
        for param in hyperParams:
            queue_in.put(param)
        self.clusterThreads.clear()
        for _ in range(n_jobs):
            p = Process(target=self.clusterThread, args=(queue_in, lock_in))
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
        gc.collect(2)

    def auto(self) -> tuple:
        """

        :return:
        """
        if self.auto_params[5] == 's_dbw':
            self.states.sort(key=lambda l: l[3])
        else:
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
        if self.metric == 'si_score':
            try:
                z_min = min([x for x in z.flat if x > -1.0])
            except ValueError:
                z_min = -1.0
        else:
            z_min = _min = min(z.flat)
        z_max = max([x for x in z.flat if x < np.inf])
        self.figs['colormap'] = (y, x, z.copy().T, z_min, z_max)
        z_correct = np.array([(True if (len(set(data[3]))) > 1 else False) for data in colormap_data], dtype=bool)
        z_correct.shape = (y.size, x.size)
        if self.auto_params[5] == 's_dbw':
            V, N, = regr_cube(y, x, z, z_correct, rev=True)
        else:
            V, N, = regr_cube(y, x, z, z_correct)
        if self.auto_params[5] == 's_dbw':
            V_alt, N_alt = regr_cube_alt(y, x, z, z_correct, rev=True)
        else:
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

    def get_nsol(self, nsol: int) -> list:
        """

        :param nsol:
        :return:
        """
        if self.states and len(self.states) >= nsol:
            sols = []
            for n in range(nsol):
                state = self.states[n]
                sols.append((n, state[3], state[2], state[4], state[5]))
        else:
            sols = None
        return sols

    def get_nsol_ext(self, nsol: int) -> list:
        """

        :param nsol:
        :return:
        """
        if self.figs:
            y, x, z = self.figs['colormap'][0], self.figs['colormap'][1], self.figs['colormap'][2]
            if self.auto_params[5] == 'si_score':
                ext = 'max'
                bg = -1
                z_cut = 0
            elif self.auto_params[5] == 's_dbw':
                ext = 'min'
                bg = np.inf
                z_cut = 0
            else:
                ext = 'max'
                bg = 0
                z_min = z.min()
                z_max = z.max()
                z_cut = (z_max + z_min) / 2
            ext_sols = detect_local_extrema(z, x, y, ext=ext, bg=bg, n=nsol)
            sols = []
            n_cls = [(state[5], state[4], state[2]) for state in self.states]
            for n, (val, min_samples, eps, xind, yind) in enumerate(ext_sols):
                if val >= z_cut:
                    ncl = 0
                    for st in n_cls:
                        if st[0] == min_samples and st[1] == eps:
                            ncl = st[2]
                    sols.append((n, val, ncl, eps, min_samples))
        else:
            sols = None
        return sols

    def noise_percent(self) -> float:
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

        :param url:
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

    def get_protein_info(self):
        """

        :return:
        """
        try:
            from Bio.PDB import PDBParser, PDBList
            from Bio.PDB.Polypeptide import PPBuilder
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
        except ImportError:
            raise ImportError('Biopython is not installed!')
        parser = PDBParser()
        with io.StringIO() as f:
            f.writelines(self.s_array)
            f.flush()
            f.seek(0, 0)
            structure = parser.get_structure('X', f)
        name = parser.header.get('name', '')
        head = parser.header.get('head', '')
        method = parser.header.get('structure_method', '')
        res = parser.header.get('resolution', '')
        ncomp = 0
        eclist = []
        for values in parser.header['compound'].values():
            ncomp += 1
            eclist.append(values.get('ec', '') or values.get('ec_number', ''))
        ec = ", ".join([ec for ec in eclist if ec])
        nres = 0
        mmass = 0
        nchain = 0
        ppb = PPBuilder()
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()
            nres += len(seq)
            seqan = ProteinAnalysis(str(seq))
            mmass += int(seqan.molecular_weight())
            nchain += 1
        return name, head, method, res, ncomp, nchain, ec, nres, mmass

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
        # www.pnas.org/cgi/doi/10.1073/pnas.1616138113 # 1.0 - -9.58 kj/mol (ALA) A; residues with delta mu < 0
        nanodroplet = {'ALA': 1.0, 'VAL': 0.862, 'PRO': 0.788, 'LEU': 0.904, 'ILE': 1.016, 'PHE': 0.963, 'MET': 0.799,
                       'TRP': 0.900, 'CYS': 0.588, 'GLY': 0.477, 'THR': 0.424, 'SER': 0.372}
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
        table_group_flag = False
        if htable == 'hydropathy':
            hydrfob = hydropathy
        elif htable == 'menv':
            hydrfob = menv
        elif htable == 'fuzzyoildrop':
            hydrfob = fuzzyoildrop
        elif htable == 'nanodroplet':
            hydrfob = nanodroplet
        elif htable == 'aliphatic_core':
            hydrfob = aliphatic_core
        elif htable == 'hydrophilic':
            hydrfob = hydropathy_h2o
        elif htable == 'positive' or htable == 'negative':
            hydrfob = calc_abs_charge(htable, pH)
        elif htable == 'ngroup' or htable == 'pgroup':
            table_group_flag = True
            group_table = calc_group_charge(htable, pH)
            hydrfob = set(group_table.keys())
        elif htable == 'rekkergroup':
            # Raimund Mannhold, Roelof F. Rekker
            # 'The hydrophobic fragmental constant approach for calculating log P in octanol/water and aliphatic
            # hydrocarbon/water systems' Perspectives in Drug Discovery and Design, 18: 1–18, 2000.
            table_group_flag = True
            group_table = {'ALA': ((('CA',), 0.315),
                                   (('CB',), 0.519)),
                           'ARG': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519),
                                   (('CD',), 0.519)),
                           'ASP': ((('CA',), 0.315),
                                   (('CB',), 0.519)),
                           'ASN': ((('CA',), 0.315),
                                   (('CB',), 0.519)),
                           'CYS': ((('CA',), 0.315),
                                   (('CB',), 0.519)),
                           'GLU': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519)),
                           'GLN': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519)),
                           'HIS': ((('CA',), 0.315),
                                   (('CB',), 0.519),),
                           'ILE': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG1',), 0.724),
                                   (('CG2',), 0.519),
                                   (('CD1',), 0.724)),
                           'LEU': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.315),
                                   (('CD1',), 0.724),
                                   (('CD2',), 0.724)),
                           'LYS': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519),
                                   (('CE',), 0.519),
                                   (('CD',), 0.519)),
                           'MET': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519),
                                   (('CE',), 0.519)),
                           'PRO': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), 0.519),
                                   (('CD',), 0.519)),
                           'SER': ((('CA',), 0.315),
                                   (('CB',), 0.519)),
                           'THR': ((('CA',), 0.315),
                                   (('CB',), 0.110),
                                   (('CG2',), 0.724)),
                           'PHE': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), ('CD1',), ('CD2',), ('CE1',), ('CE2',), ('CZ',), 1.903)),
                           'TRP': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), ('CD1',), ('CD2',), ('CE2',),
                                    ('CE3',), ('CZ2',), ('CZ3',), ('CH2',), ('NE1',), 1.902)),
                           'TYR': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG',), ('CD1',), ('CD2',), ('CE1',), ('CE2',), ('CZ',), 1.903)),
                           'VAL': ((('CA',), 0.315),
                                   (('CB',), 0.519),
                                   (('CG1',), 0.724),
                                   (('CG2',), 0.724))
                           }
            hydrfob = set(group_table.keys())
        xyzm_array = []
        current_resn = None
        current_chainn = None
        current_resname = None
        atom_list = []
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
                atom_list.append(s[12:16].strip().upper())
                xyzm_array = np.hstack((xyzm_array, xyzm))
            elif s[0:6] == 'ATOM  ' and (s[17:20] in hydrfob):
                try:
                    xyzm_array.shape = (-1, 4)
                except AttributeError:
                    raise ValueError
                if table_group_flag:
                    if create_group(group_table, current_resname, xyzm_array, atom_list):
                        groups_xyz_array, weights, real_groups = create_group(group_table, current_resname, xyzm_array,
                                                                              atom_list)
                        for group in real_groups:
                            self.aa_list.append((current_resn, current_chainn, current_resname, group))
                        self.weight_array.extend(weights)
                        xyz_array = np.concatenate((xyz_array, groups_xyz_array))
                else:
                    self.aa_list.append((current_resn, current_chainn, current_resname, ()))
                    self.weight_array.append(hydrfob[current_resname])
                    xyz_array = np.concatenate((xyz_array, cmass(xyzm_array)))
                xyzm_array = []
                atom_list.clear()
                current_resn = int(s[22:26])
                current_chainn = s[21]
                current_resname = s[17:20]
                xyz = [float(s[30:38]), float(s[38:46]), float(s[46:54]), mass[s[76:78]]]
                xyzm_array = np.hstack((xyzm_array, xyz))
                atom_list.append(s[12:16].strip().upper())
        try:
            xyz_array.shape = (-1, 3)
        except AttributeError:
            raise ValueError
        pdist = euclidean_distances(xyz_array)
        parse_results = len(self.aa_list), np.min(pdist[np.nonzero(
            pdist)]), np.max(pdist[np.nonzero(pdist)]), np.mean(pdist[np.nonzero(pdist)])
        # Attension! This code makes problems and it removed. 1) In this case size of sparse matrix  and size of
        # pdist matrix are equivalents. 2) In dbscan (scikit-learn) founded bug if input array this precomputed
        # sparse neighbors graph. See https://github.com/scikit-learn/scikit-learn/pull/12105
        # sparse_n = NearestNeighbors(radius=parse_results[2], algorithm='brute', n_jobs=-1
        #                             ).fit(xyz_array).radius_neighbors_graph(xyz_array, mode='distance')
        self.X, self.pdist, self.parse_results = xyz_array, pdist, parse_results
        return parse_results

    def graph(self, grid_state: bool, legend_state: bool) -> tuple:
        """

        :return:
        """
        fig = Figure(figsize=(8, 6))
        try:
            unique_labels = sorted(set(self.labels))
            xyz_all = self.X
        except TypeError:
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
            es, ms, z, z_min, z_max = self.figs['colormap']
        else:
            raise ValueError
        fig = Figure(figsize=(12, 6))
        ax1 = fig.add_subplot(121)
        ax1.set_title(self.metrics_name[self.metric])
        ax1.set_xlabel('EPS, \u212B')
        ax1.set_ylabel('MIN SAMPLES')
        ax1.grid(grid_state)
        if self.auto_params[5] == 's_dbw':
            ax1.patch.set_facecolor("black")
            pc1 = ax1.pcolor(es, ms, z, cmap='gnuplot_r', vmin=z_min, vmax=z_max)
        else:
            pc1 = ax1.pcolor(es, ms, z, cmap='gnuplot', vmin=z_min, vmax=z_max)
        fig.colorbar(pc1, ax=ax1, extend='max', extendfrac=0.1)
        V, N, Nfit, C, B, R2 = self.figs['linear']
        ax11 = fig.add_subplot(122)
        ax11.set_ylabel('MIN SAMPLES')
        ax11.set_xlabel(r'$V,\ \AA^3$')
        ax11.scatter(V, N, c='k')
        vs = [Vcoeff * state[3] ** 3 for state in self.get_nsol(5)]
        ns = [state[4] for state in self.get_nsol(5)]
        ax11.scatter(vs, ns, c='g', label='Suboptimal values')
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

    def colormap3d(self, grid_state: bool = True):
        """

        :param grid_state:
        :return:
        """
        if not self.states:
            raise ValueError
        if self.figs is not None:
            es, ms, z, z_min, z_max = self.figs['colormap']
        else:
            raise ValueError
        es, ms = np.meshgrid(es, ms)
        fig = Figure(figsize=(8, 6))
        ax = axes3d.Axes3D(fig)
        surf = ax.plot_surface(es, ms, z, antialiased=True, cmap='coolwarm', vmin=z_min, vmax=z_max)
        ax.set_title(self.metrics_name[self.metric])
        ax.set_xlabel('EPS, \u212B')
        ax.set_ylabel('MIN SAMPLES')
        ax.set_zlabel(self.metrics_name[self.metric])
        ax.grid(grid_state)
        fig.colorbar(surf, shrink=0.5, aspect=5)
        return ax, fig

    def fig_scan_param(self, mode: str, value):
        """

        :param mode:
        :param value:
        :return:
        """
        if not self.states:
            raise ValueError
        states = self.states
        metric = self.metrics_name[self.auto_params[5]]
        htable = self.htable
        if mode == 'min_samples':
            x, y, y1, y2, y3, y4 = calculate_scan(states, value, 'eps')
            xparametr = 'MIN_SAMPLES'
            const_str = 'EPS: {:.2f} \u212B'.format(value)
        else:
            x, y, y1, y2, y3, y4 = calculate_scan(states, value, 'min_samples')
            xparametr = 'EPS (\u212B)'
            const_str = 'MIN_SAMPLES: {:d}'.format(value)
        fig = draw_scan_param(x, y, y1, y2, y3, y4, htable, metric, xparametr, const_str)
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
                                 "+".join(['(chain {1:s} and resi {0:d}){2:s}'.format(
                                     aac[0], aac[1],
                                     ' and name {:s}'.format('+'.join(aac[3])) if (len(aac) > 3 and aac[3]) else '')
                                     for aac in aa_list])))
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
                        ['(chain {1:s} and resi {0:d}{2:s})'.format(aac[0], aac[1],
                                                                    ' and name {:s}'.format('+'.join(aac[3])) if (
                                                                            len(aac) > 3 and aac[3]) else '')
                         for aac in aa_list]))
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
            'X'                : self.X,
            'pdist'            : self.pdist,
            'labels'           : self.labels,
            'noise_filter'     : self.noise_filter,
            'core_samples_mask': self.core_samples_mask,
            'n_clusters'       : self.n_clusters,
            'score'            : self.score,
            'eps'              : self.eps,
            'min_samples'      : self.min_samples,
            'metric'           : self.metric,
            'weight_array'     : self.weight_array,
            'aa_list'          : self.aa_list,
            's_array'          : self.s_array,
            'htable'           : self.htable,
            'parse_results'    : self.parse_results,
            'auto_params'      : self.auto_params,
            'states'           : self.states,
            'figs'             : self.figs,
            'version'          : __version__}
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
        noise_filter = global_state['noise_filter']
        if noise_filter is True:  # For 0.1 version saves compatibility
            self.noise_filter = 'filter'
        elif noise_filter is False:
            self.noise_filter = 'comb'
        else:
            self.noise_filter = noise_filter
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
        return global_state.get('version', '0.2.0d50')
