#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 01.06.18"""
import os
import os.path
import sys
from queue import Queue
from threading import Lock, Thread
from .testlist import opendb, download, clusterThread


def main(namespace):
    """

    :param namespace:
    """
    output = namespace.output
    inp = namespace.input
    if not output:
        print("Output name is empty!")
        sys.exit()
    else:
        outputDir = '{:s}_data'.format(output)
    if not inp:
        print("Input name is empty!")
        sys.exit()
    min_eps = namespace.emin
    max_eps = namespace.emax
    step_eps = namespace.estep
    min_min_samples = namespace.smin
    max_min_samples = namespace.smax
    try:
        os.makedirs(outputDir, exist_ok=True)
    except OSError:
        print('Unable to create folder', outputDir)
        sys.exit()
    try:
        with open(inp) as inpf:
            filelist = list(map(lambda x: x.strip().lower(), inpf.readlines()))
            if not filelist:
                raise ValueError
    except (UnicodeDecodeError, ValueError):
        print("File {:s} does not contain text!".format(inp))
        sys.exit()
    except (OSError, FileNotFoundError):
        print("File {:s} is unavailable!".format(inp))
        sys.exit()
    q = Queue()
    lock = Lock()
    cursor, conn = opendb(output)
    dTask = Thread(target=download, args=(filelist, q, lock, cursor, conn, outputDir))
    dTask.start()
    n = 1
    while True:
        item = q.get()
        if item is None:
            break
        clusterThread(item, outputDir, cursor, conn, lock, min_eps, max_eps,
                      step_eps, min_min_samples, max_min_samples, n_jobs=0)
        print("All tasks completed for {:s} ({:d}/{:d})".format(item, n, len(filelist)))
        n += 1
    dTask.join()
    conn.close()
