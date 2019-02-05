# -*- coding: utf-8 -*-
"""Created by lashkov on 31.10.18"""

import math

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import check_X_y


def calc_nearest_points(X: np.ndarray, labels: np.ndarray, k: int, centroids: np.ndarray, metric: str) -> np.ndarray:
    """
    Calculation of coordinates of clusters points closest to their geometric centers

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.
    k : int,
        No. of clusters (k > 0)
    centroids : array-like, shape (n_clusters, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single centroid.
    metric : str,
        The distance metric, can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
        ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’,
        ‘yule’. Default is ‘euclidean’.

    Returns
    -------
    centroids : array-like, shape (n_clusters, n_features)
        List of n_features-dimensional points. Each row corresponds
        to a single centroid.
    """
    centroids_p = []
    for i in range(k):
        dist = cdist(X[labels == i], np.array(centroids[i], ndmin=2), metric=metric)
        centroids_p.append(X[labels == i][np.argmin(dist)])
    return np.array(centroids_p)


def calc_centroids(X: np.ndarray, k: int, labels: np.ndarray, centr: str) -> np.ndarray:
    """
    Calculation of coordinates of centroids

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    k : int,
        No. of clusters (k > 0)
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample.
    centr : str,
        Cluster center calculation method (mean (default) or median)

    Returns
    -------
    centroids : array-like, shape (n_samples', n_features)
         List of n_features-dimensional data points. Each row corresponds
         to a single data point - center if a cluster.
    """
    centers = []
    for i in range(k):
        if centr == "mean":
            centers.append(np.mean(X[labels == i], axis=0))
        elif centr == "median":
            centers.append(np.median(X[labels == i], axis=0))
    return np.array(centers)


def density(X: np.ndarray, centroids: np.ndarray, labels: np.ndarray,
            stdev: float, clusters_list: list, method: str, density_list: list = None, lambd: float = 0.7) -> float:
    """
    Compute the density of one or two cluster(depend on cluster_list)

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    centroids : array-like, shape (n_samples', n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single centroid.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample (-1 - for noise).
    stdev : float,
        Average standard deviation of clusters (for Halkidi method [1])
    clusters_list : list,
        List contains one or two No. of cluster
    method : str,
        S_Dbw calc method:
        'Halkidi' - original paper [1]
        'Kim' - see [2]
        'Tong' - see [3]
    density_list : list,
        List contains density of each cluster for calculate muij [3]
    lambd : float,
        Lambda coefficient is a positive constant between 0 and 1 (default: 0.7, see [3]

    Returns
    -------
    score : float
        Density of one or two cluster

    References:
    -----------
    .. [1] M. Halkidi and M. Vazirgiannis, “Clustering validity assessment: Finding the optimal partitioning
        of a data set,” in ICDM, Washington, DC, USA, 2001, pp. 187–194.
    .. [2] Youngok Kim and Soowon Lee. A clustering validity assessment Index. PAKDD’2003, Seoul, Korea, April 30–May 2,
        2003, LNAI 2637, 602–608
    .. [3] Tong, J. & Tan, H. J. Electron.(China) (2009) 26: 258. https://doi.org/10.1007/s11767-007-0151-8
    """

    density = 0
    center_p1 = centroids[clusters_list[0]]
    if len(clusters_list) == 2:
        center_p2 = centroids[clusters_list[1]]
        if method == 'Kim' or method == 'Tong':
            sigmai = np.std(X[labels == clusters_list[0]], axis=0)
            sigmaj = np.std(X[labels == clusters_list[1]], axis=0)
            sigmaij = (sigmai + sigmaj) / 2
            ni = X[labels == clusters_list[0]].shape[0]
            nj = X[labels == clusters_list[1]].shape[0]
            nij = ni + nj
            if method == 'Tong':
                center_v = lambd * (center_p1 * nj + center_p2 * ni) / nij + \
                           (1 - lambd) * ((center_p1 * density_list[clusters_list[0]] +
                                           center_p2 * density_list[clusters_list[1]]) /
                                          (density_list[clusters_list[0]] + density_list[clusters_list[1]]))
            else:
                center_v = (center_p1 + center_p2) / 2
        else:
            center_v = (center_p1 + center_p2) / 2
    else:
        center_v = center_p1
        if method == 'Kim' or method == 'Tong':
            sigmaij = np.std(X[labels == clusters_list[0]], axis=0)
            nij = X[labels == clusters_list[0]].shape[0]
    if method == 'Halkidi':
        for i in clusters_list:
            temp = X[labels == i]
            for j in temp:
                if np.linalg.norm(j - center_v) <= stdev:
                    density += 1
    elif method == 'Kim' or method == 'Tong':
        CI = 1.96 * sigmaij / math.sqrt(nij)
        for i in clusters_list:
            temp = X[labels == i]
            for j in temp:
                if np.all(np.abs(j - center_v) <= CI):
                    density += 1
    return density


def Dens_bw(X: np.ndarray, centroids: np.ndarray, labels: np.ndarray, k: int, method: str = 'Halkidi') -> float:
    """
    Compute Inter-cluster Density (ID) - It evaluates the average density in the region among clusters in relation
    with the density of the clusters. The goal is the density among clusters to be significant low in comparison with
    the density in the considered clusters [1].

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    centroids : array-like, shape (n_samples', n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single centroid.
    k : int,
        No. of clusters
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample (-1 - for noise).
    method : str,
        S_Dbw calc method:
        'Halkidi' - original paper [1]
        'Kim' - see [2]
        'Tong' - see [3]

    Returns
    -------
    score : float
        Inter-cluster Density

    References:
    -----------
    .. [1] M. Halkidi and M. Vazirgiannis, “Clustering validity assessment: Finding the optimal partitioning
        of a data set,” in ICDM, Washington, DC, USA, 2001, pp. 187–194.
    .. [2] Youngok Kim and Soowon Lee. A clustering validity assessment Index. PAKDD’2003, Seoul, Korea, April 30–May 2,
        2003, LNAI 2637, 602–608
    .. [3] Tong, J. & Tan, H. J. Electron.(China) (2009) 26: 258. https://doi.org/10.1007/s11767-007-0151-8
    """
    density_list = []
    result = 0
    stdev = 0
    if method == 'Halkidi':
        for i in range(k):
            std_matrix_i = np.std(X[labels == i], axis=0)
            stdev += math.sqrt(np.dot(std_matrix_i.T, std_matrix_i))
        stdev = math.sqrt(stdev) / k
    for i in range(k):
        density_list.append(density(X, centroids, labels, stdev, [i], method))
    if density_list.count(0) > 1:
        raise ValueError('The density for two or more clusters to equal zero.')
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            result += density(X, centroids, labels, stdev, [i, j], method, density_list) / max(density_list[i],
                                                                                               density_list[j])
    return result / (k * (k - 1))


def Scat(X: np.ndarray, k: int, labels: np.ndarray, method: str) -> float:
    """
    Calculate intra-cluster variance (Average scattering for clusters).
    Lower value -> better clustering.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.
    k : int,
        No. of clusters
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample (-1 - for noise).
    method : str,
        S_Dbw calc method:
        'Halkidi' - original paper [1]
        'Kim' - see [2]
        'Tong' - see [3]

    Returns
    -------
    score : float
        Average scattering for clusters

    References:
    -----------
    .. [1] M. Halkidi and M. Vazirgiannis, “Clustering validity assessment: Finding the optimal partitioning
        of a data set,” in ICDM, Washington, DC, USA, 2001, pp. 187–194.
    .. [2] Youngok Kim and Soowon Lee. A clustering validity assessment Index. PAKDD’2003, Seoul, Korea, April 30–May 2,
     2003, LNAI 2637, 602–608
    .. [3] Tong, J. & Tan, H. J. Electron.(China) (2009) 26: 258. https://doi.org/10.1007/s11767-007-0151-8
    """
    theta_s = np.std(X, axis=0)
    theta_s_2norm = math.sqrt(np.dot(theta_s.T, theta_s))
    sum_theta_2norm = 0
    if method == 'Halkidi':
        for i in range(k):
            theta_i = np.std(X[labels == i], axis=0)
            sum_theta_2norm += math.sqrt(np.dot(theta_i.T, theta_i))
        result = sum_theta_2norm / (theta_s_2norm * k)
    else:
        n = len(labels)
        for i in range(k):
            ni = X[labels == i].shape[0]
            theta_i = np.std(X[labels == i], axis=0)
            sum_theta_2norm += ((n - ni) / n) * math.sqrt(np.dot(theta_i.T, theta_i))
        result = sum_theta_2norm / (theta_s_2norm * k)
        if method == 'Tong':
            result = sum_theta_2norm / (theta_s_2norm * (k - 1))
    return result


def S_Dbw(X: np.ndarray, labels: np.ndarray,
          method: str = 'Tong', centr: str = 'mean', metric: str = 'euclidean') -> float:
    """
    Compute the S_Dbw validity index
    S_Dbw validity index is defined by equation:
    S_Dbw = scatt + dens
    where scatt - means average scattering for clusters and dens - inter-cluster density.
    Lower value -> better clustering.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (n_samples,)
        Predicted labels for each sample (-1 - for noise).
    method : str,
        S_Dbw calc method:
        'Halkidi' - original paper [1]
        'Kim' - see [2]
        'Tong' - see [3]
    centr : str,
        cluster center calculation method (mean (default) or median)
    metric : str,
        The distance metric, can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
        ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘wminkowski’,
        ‘yule’. Default is ‘euclidean’.

    Returns
    -------
    score : float
        The resulting S_DBw score.

    References:
    -----------
    .. [1] M. Halkidi and M. Vazirgiannis, “Clustering validity assessment: Finding the optimal partitioning
        of a data set,” in ICDM, Washington, DC, USA, 2001, pp. 187–194.
    .. [2] Youngok Kim and Soowon Lee. A clustering validity assessment Index. PAKDD’2003, Seoul, Korea, April 30–May 2,
        2003, LNAI 2637, 602–608
    .. [3] Tong, J. & Tan, H. J. Electron.(China) (2009) 26: 258. https://doi.org/10.1007/s11767-007-0151-8
    """
    if len(set(labels)) < 2 or len(set(labels)) > len(X) - 1:
        raise ValueError("No. of unique labels must be > 1 and < n_samples")
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    k = len(le.classes_)
    centroids = calc_centroids(X, k, labels, centr)
    centroids = calc_nearest_points(X, labels, k, centroids, metric)
    if k < 2:
        raise ValueError('Only one cluster!')
    sdbw = Dens_bw(X, centroids, labels, k, method) + Scat(X, k, labels, method)
    return sdbw
