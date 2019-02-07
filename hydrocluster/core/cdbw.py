# -*- coding: utf-8 -*-

import importlib
import math
from collections import defaultdict

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import check_X_y


def gen_dist_func(metric):
    """
    Obtain the distances function from scipy.spatial.distance package

    Parameters
    ----------
    metric : str,
        The distance metric, can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
        ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’,
        ‘yule’.

    Returns
    -------
    func : function (array_like, array_like)
        Function for calculation distance between two points in n-dimensional space.
    """
    mod = importlib.import_module("scipy.spatial.distance")
    func = getattr(mod, metric)
    return func


def prep(X, labels):
    """
    Calculation necessary parameters

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.

    Returns
    -------
    n_clusters : int,
        Number of clusters.
    stdev : float,
        Average standard deviation of the considered clusters.
    dimension : int,
        Dimension of Data Set.
    """
    dimension = X.shape[1]
    n_clusters = labels.max() + 1
    stdev = 0
    for i in range(n_clusters):
        std_matrix_i = np.std(X[labels == i], axis=0)
        stdev += math.sqrt(np.dot(std_matrix_i.T, std_matrix_i))
    stdev = math.sqrt(stdev) / n_clusters
    return n_clusters, stdev, dimension


def rep(X, labels, n_clusters, dimension):
    """
    Select of representative points for each clusters

    Parameters
    ----------

    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.
    n_clusters : int,
        Number of clusters.
    dimension : int,
        Dimension of Data Set.


    Returns
    -------
    rep_dic : dict {i: indexes}
        Indexes of representative points for each clusters, i - No. of cluster, indexes - array of indexes.
    mean_arr : array_like shape (n_clusters, dimension)
        Coordinates of the centroid of each cluster.
    n_rep : array_like shape (n_clusters,)
        Number of representative points in each cluster.
    n_points_in_cl : array_like shape (n_clusters,)
        Number of points in each cluster
    """
    rep_dic = {}
    mean_arr = np.zeros(shape=(n_clusters, dimension), dtype=float)
    n_rep = np.zeros(shape=(n_clusters,), dtype=int)
    labels_in_cluster = {}
    n_points_in_cl = np.zeros(shape=(n_clusters,), dtype=int)
    for i in range(n_clusters):
        labels_in_cluster[i] = np.where(labels == i)
        if labels_in_cluster[i][0].size >= 4:
            ch = ConvexHull(X[labels == i])
            rep_dic[i] = labels_in_cluster[i][0][ch.vertices]
        else:
            rep_dic[i] = labels_in_cluster[i][0]
        mean_arr[i] = np.mean(X[labels == i], axis=0)
        n_rep[i] = rep_dic[i].size
        n_points_in_cl[i] = labels_in_cluster[i][0].size
    return rep_dic, mean_arr, n_rep, n_points_in_cl


def closest_rep(X, n_clusters, rep_dic, n_rep, metric, distvec):
    """
    Select of the closest representative points for two clusters

    Parameters
    ----------

    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    n_clusters : int,
        Number of clusters.
    rep_dic : dict {number: indexes}
        Indexes of representative points for each clusters, number - No. of cluster, indexes - array of indexes.
    n_rep : array_like shape (n_clusters,)
        Number of representative points in each cluster.
    metric : str,
        The distance metric, can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
        ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’,
        ‘yule’.
    distvec : function (array_like, array_like)
        Function for calculation distance between two points in n-dimensional space.

    Returns
    -------
    dist_min : defaultdict {tuple (i, j) : list [float]}
        List of distances between of of closest representative points for each pair of clusters,
        i - No. of 1 cluster, j - No. of 2 cluster.
    middle_point : defaultdict {tuple (i, j) : list [array_like (dimension)]}
        List of coordinates of middle points for two closest representative points of each pair of clusters,
        i - No. of 1 cluster, j - No. of 2 cluster.
    n_cl_rep : dict {tuple (i, j): indexes}
        Indices of closest representative points for each two clusters, i - No. of 1 cluster, j - No. of 2 cluster
        indexes - array of indexes.

    """
    b1 = {}
    b2 = {}
    dist_arr = {}
    min_value1 = {}
    min_value0 = {}
    min_index0 = {}
    min_index1 = {}
    min_index2 = {}
    min_index_r = {}
    min_index_c = {}
    cl_r = {}
    s1 = []
    s2 = []
    s2_t = []
    t1 = []
    t2 = []
    v1 = []
    v2 = []
    dist_min = defaultdict(list)
    middle_point = defaultdict(list)
    n_cl_rep = {}
    for i in range(n_clusters):
        for j in range(n_clusters):
            if i > j:
                dist_arr[(i, j)] = cdist(X[rep_dic[i]], X[rep_dic[j]], metric=metric)
                min_value1[(i, j)] = dist_arr[(i, j)].min(axis=1)
                min_value0[(i, j)] = dist_arr[(i, j)].min(axis=0)
                min_index1[(i, j)] = dist_arr[(i, j)].argmin(axis=1)
                min_index0[(i, j)] = dist_arr[(i, j)].argmin(axis=0)
                min_index_r[(i, j)] = np.add(np.arange(0, n_rep[i] * n_rep[j], n_rep[j]), min_index1[(i, j)])
                min_index_c[(i, j)] = np.add(np.arange(0, n_rep[i] * n_rep[j], n_rep[i]), min_index0[(i, j)])
                t1 += [n_rep[i]]
                t2 += [n_rep[j]]
                for k in range(n_rep[i]):
                    s1.append(np.unravel_index(min_index_r[(i, j)][k], (n_rep[i], n_rep[j])))
                for n in range(n_rep[j]):
                    s2.append(np.unravel_index(min_index_c[(i, j)][n], (n_rep[j], n_rep[i])))
                    s2_t = [(x[1], x[0]) for x in s2]
    p = 0
    for m in range(len(t1)):
        p += t1[m]
        v1.append(p)
    p = 0
    for m in range(len(t2)):
        p += t2[m]
        v2.append(p)
    min_index1[(1, 0)] = s1[0:v1[0]]
    min_index2[(1, 0)] = s2_t[0:v2[0]]
    l = 0
    for i in range(2, n_clusters):
        for j in range(n_clusters):
            if i > j:
                min_index1[(i, j)] = s1[v1[l]:v1[l + 1]]
                min_index2[(i, j)] = s2_t[v2[l]:v2[l + 1]]
                l += 1
    for i in range(n_clusters):
        for j in range(n_clusters):
            if i > j:
                b1[(i, j)] = set(min_index1[(i, j)])
                b2[(i, j)] = set(min_index2[(i, j)])
                cl_r[(i, j)] = list(b1[(i, j)] & b2[(i, j)])
                n_cl_rep[(i, j)] = len(cl_r[i, j])
                for u, v in cl_r[(i, j)]:
                    middle_point[(i, j)].append((X[rep_dic[i][u]] + X[rep_dic[j][v]]) / 2)
                    dist_min[(i, j)].append(distvec(X[rep_dic[i][u]], X[rep_dic[j][v]]))
    return middle_point, dist_min, n_cl_rep


def art_rep(X, n_clusters, rep_dic, n_rep, mean_arr, s):
    """
    Calculate of the art representative points

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    n_clusters : int,
        Number of clusters.
    rep_dic : dict {i: indexes}
        Indexes of representative points for each clusters, i - No. of cluster, indexes - array of indexes.
    n_rep : array_like shape (n_clusters,)
        Number of representative points in each cluster.
    mean_arr : array_like shape (n_clusters, dimension)
        Coordinates of the centroid of each cluster.
    s : int,
        Number of art representative points. (>2)

    Returns
    -------
    a_rep_shell : defaultdict {tuple (i, k) : list [array_like (dimension)]}
        List of n-dimensional coordinates of art representative points for pair i, k where i - No. of cluster,
        k - No. of shell.
    """
    a_rep_shell = defaultdict(list)
    for i in range(n_clusters):
        if n_rep[i] == 1:
            raise ValueError('Cluster No. {:d} obtain only 1 point'.format(i))
        for x in rep_dic[i]:
            for k in range(1, s + 1):
                a_rep_shell[i, k - 1].append((1 - (k / s)) * X[x] + (k / s) * mean_arr[i])
    return a_rep_shell


def compactness(X, labels, n_clusters, stdev, a_rep_shell, n_rep, n_points_in_cl, distvec, s, ms):
    """
    Clusters compactness and cohesion evaluation

    Parameters
    ----------

    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.
    n_clusters : integer
        Number of clusters.
    stdev : float
        Average standard deviation of the considered clusters.
    a_rep_shell : defaultdict {tuple (i, k) : list [array_like (dimension)]}
        List of n-dimensional coordinates of art representative points for pair i, k where i - No. of cluster,
        k - No. of shell.
    n_rep : array_like shape (n_clusters,)
        Number of representative points in each cluster.
    n_points_in_cl : array_like shape (n_clusters,)
        Number of points in each cluster
    distvec : function (array_like, array_like)
        Function for calculation distance between two points in n-dimensional space.
    s : integer
        Number of art representative points. (>2)

    Returns
    -------
    compact : float
        Compactness of clusters.
    cohesion : float
        Cohesion of clusters.

    """
    card = defaultdict(lambda: 0)
    a_rep_shell1 = {}
    intra_dens_shell = np.zeros(shape=(n_clusters, s), dtype=float)
    for i in range(n_clusters):
        for x in X[labels == i]:
            for k in range(s):
                a_rep_shell1[i, k] = np.array(a_rep_shell[i, k])
                for p in a_rep_shell1[i, k]:
                    dist = distvec(x, p)
                    if dist < stdev:
                        card[i, k] += 1
    for i in range(n_clusters):
        for k in range(s):
            intra_dens_shell[i, k] = card[i, k] / (n_rep[i] * n_points_in_cl[i])
    intra_dens = np.sum(intra_dens_shell, axis=0) / (stdev * n_clusters)
    compact = np.sum(intra_dens) / s
    intra_change = 0
    for l in range(s - 1):
        intra_change += abs(intra_dens[l + 1] - intra_dens[l]) / (s - 1)
    cohesion = compact / (1 + intra_change)
    return compact, cohesion


def separation(X, labels, n_clusters, stdev, middle_point, dist_min, n_cl_rep, n_points_in_cl, distvec):
    """
    Clusters separation evaluation

    Parameters
    ----------

    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.
    n_clusters : int,
        Number of clusters.
    stdev : float.
        Average standard deviation of the considered clusters.
    middle_point : defaultdict {tuple (i, j) : list [array_like (dimension)]},
        List of coordinates of middle points for two closest representative points of each pair of clusters,
        i - No. of 1 cluster, j - No. of 2 cluster.
    dist_min : defaultdict {tuple (i, j) : list [float]}
        List of distances between of of closest representative points for each pair of clusters,
        i - No. of 1 cluster, j - No. of 2 cluster.
    n_cl_rep : dict {tuple (i, j): indexes}
        Indices of closest representative points for each two clusters, i - No. of 1 cluster, j - No. of 2 cluster
        indexes - array of indexes.
    n_points_in_cl : array_like shape (n_clusters,)
        Number of points in each cluster
    distvec : function (array_like, array_like)
        Function for calculation distance between two points in n-dimensional space.

    Returns
    -------
    sep : float,
        Separation of clusters.
    """
    dens_mean = np.zeros((n_clusters, n_clusters))
    dist_mm = np.zeros((n_clusters, n_clusters))
    dist_mm[np.diag_indices_from(dist_mm)] = np.inf
    card1 = {k: [0 for _ in range(n)] for k, n in n_cl_rep.items()}
    for i in range(n_clusters):
        for j in range(n_clusters):
            if i > j:
                for s in range(n_cl_rep[i, j]):
                    for x in np.array(middle_point[(i, j)][s]):
                        for p in np.array(np.vstack([X[labels == i], X[labels == j]])):
                            dist1 = distvec(x, p)
                            if dist1 < stdev:
                                card1[i, j][s] += 1
    for i in range(n_clusters):
        for j in range(n_clusters):
            if i > j:
                dens_mean[i, j] = np.mean(np.array(dist_min[i, j]) * np.array(card1[i, j]))
                dist_mm[i, j] = np.mean(dist_min[i, j])
            elif i < j:
                dens_mean[i, j] = np.mean(np.array(dist_min[j, i]) * np.array(card1[j, i]))
                dist_mm[i, j] = np.mean(dist_min[j, i])
            dens_mean[i, j] /= (n_points_in_cl[i] + n_points_in_cl[j])
    inter_dens = np.sum(np.max(dens_mean, axis=0)) / (2 * n_clusters * stdev)
    dist_m = np.sum(np.min(dist_mm, axis=0)) / n_clusters
    sep = dist_m / (1 + inter_dens)
    return sep


def CDbw(X, labels, metric="euclidean", s=3, eps=None):
    """
    Calculate CDbw-index for cluster validation, as defined in [1]

    CDbw = compactness*cohesion*separation

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.  (-1 - for noise)
    metric : str,
        The distance metric, can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
        ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’,
        ‘yule’.
    s : int,
        Number of art representative points. (>2)
    Returns
    -------
    cdbw : float,
        The resulting CDbw validity index.

    References:
    -----------
    .. [1] M. Halkidi and M. Vazirgiannis, “A density-based cluster validity approach using multi-representatives”
        Pattern Recognition Letters 29 (2008) 773–786.

    """
    if len(set(labels)) < 2 or len(set(labels)) > len(X) - 1:
        raise ValueError("No. of unique labels must be > 1 and < n_samples")
    if s < 2:
        raise ValueError("Parameter s must be > 2")
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    distvec = gen_dist_func(metric)
    n_clusters, stdev, dimension = prep(X, labels)
    if eps is not None:
        stdev = eps
    rep_dic, mean_arr, n_rep, n_points_in_cl = rep(X, labels, n_clusters, dimension)
    middle_point, dist_min, n_cl_rep = closest_rep(X, n_clusters, rep_dic, n_rep, metric, distvec)
    try:
        a_rep_shell = art_rep(X, n_clusters, rep_dic, n_rep, mean_arr, s)
    except ValueError:
        return 0
    compact, cohesion = compactness(X, labels, n_clusters, stdev, a_rep_shell, n_rep, n_points_in_cl, distvec, s)
    if np.isinf(compact) or np.isnan(compact):
        return 0
    sep = separation(X, labels, n_clusters, stdev, middle_point, dist_min, n_cl_rep, n_points_in_cl, distvec)
    cdbw = compact * cohesion * sep
    return cdbw
