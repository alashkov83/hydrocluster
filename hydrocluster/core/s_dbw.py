#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 31.10.18"""

import numpy as np
import math
from scipy.spatial.distance import cdist

from sklearn.preprocessing import LabelEncoder
from sklearn.utils import check_X_y


def center_of_mass(str_nparray):
    mass_sum = str_nparray.shape[0]
    dim = str_nparray.shape[1]
    center_m = str_nparray.sum(axis=0) / mass_sum
    center_m.shape = (-1, dim)
    return center_m


def get_center_id(data, labels):
    center_id = []
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    for i in range(n_clusters):
        c_mass = center_of_mass(data[labels == i])
        dist = cdist(data[labels == i], c_mass)
        center_id.append(np.where(np.all(data == data[labels == i][np.argmin(dist), :], axis=1))[0][0])
    return np.array(center_id)


def density(data, centers_id, labels, stdev, density_list):
    """
    compute the density of one or two cluster(depend on density_list)
    """
    density = 0
    centers_id1 = centers_id[density_list[0]]
    if len(density_list) == 2:
        centers_id2 = centers_id[density_list[1]]
        center_v = (data[centers_id1] + data[centers_id2]) / 2
    else:
        center_v = data[centers_id1]
    for i in density_list:
        temp = data[labels == i]
        for j in temp:
            if np.linalg.norm(j - center_v) <= stdev:
                density += 1
    return density


def Dens_bw(data, centers_id, labels, stdev, k):
    density_list = []
    result = 0
    for i in range(k):
        density_list.append(density(data, centers_id, labels, stdev, density_list=[i]))
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            result += density(data, centers_id, labels, stdev, [i, j]) / max(density_list[i], density_list[j])
    return result / (k * (k - 1))


def Scat(data, k, labels):
    theta_s = np.std(data, axis=0)
    theta_s_2norm = math.sqrt(np.dot(theta_s.T, theta_s))
    sum_theta_2norm = 0
    for i in range(k):
        matrix_data_i = data[labels == i]
        theta_i = np.std(matrix_data_i, axis=0)
        sum_theta_2norm += math.sqrt(np.dot(theta_i.T, theta_i))
    return sum_theta_2norm / (theta_s_2norm * k)


def prepare_data(X, labels):
    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    if len(set(labels)) < 2 or len(set(labels)) > len(X) - 1:
        raise ValueError("No. of unique labels must be > 1 and < n_samples-1")
    centers_id = get_center_id(X, labels)
    k = len(centers_id)
    stdev = 0
    for i in range(k):
        std_matrix_i = np.std(X[labels == i], axis=0)
        stdev += math.sqrt(np.dot(std_matrix_i.T, std_matrix_i))
    stdev = math.sqrt(stdev) / k
    return X, centers_id, labels, stdev, k


def S_Dbw(X, labels):
    """
    X --> raw data
    data_cluster --> The category that represents each piece of data(the number of category should begin 0)
    centers_id --> the center_id of each cluster's center
    """
    X, centers_id, labels, stdev, k = prepare_data(X, labels)
    sdbw = Dens_bw(X, centers_id, labels, stdev, k) + Scat(X, k, labels)
    return sdbw
