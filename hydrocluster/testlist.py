#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 01.06.18"""
import os
import os.path
import sqlite3
import sys
import warnings
from multiprocessing import Queue, Process, Lock

import psutil
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from matplotlib.backends.backend_agg import FigureCanvasAgg

try:
    from .pdbcluster import ClusterPdb
except ImportError:
    print('Error! Scikit-learn not installed!')
    sys.exit()

warnings.filterwarnings("ignore")
pH = 7.0
htables = ['hydropathy', 'nanodroplet', 'menv', 'fuzzyoildrop', 'positive', 'negative']
metrics = ['si_score', 'calinski']


def opendb(fiename: str):
    """

    :param fiename:
    :return:
    """
    isDbExists = os.path.exists("{:s}.db".format(fiename))
    # Создаем соединение с нашей базой данных
    conn = sqlite3.connect("{:s}.db".format(fiename), check_same_thread=False)
    # Создаем курсор - это специальный объект который делает запросы и получает их результаты
    cursor = conn.cursor()
    # cursor.execute("PRAGMA synchronous = OFF")
    # cursor.execute("PRAGMA journal_mode = MEMORY")
    if not isDbExists:
        try:
            cursor.execute("""CREATE TABLE "Structures"  (IDPDB TEXT PRIMARY KEY UNIQUE,
NAME TEXT,
HEAD TEXT,
METHOD TEXT,
RESOLUTION REAL DEFAULT 0,
NCOMP INTEGER DEFAULT 1,
NCHAIN INTEGER DEFAULT 1,
NRES INTEGER,
MMASS INTEGER,
EC TEXT)""")
            cursor.execute("""CREATE TABLE "Results" (ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
IDPDB TEXT NOT NULL,
HTABLE TEXT NOT NULL,
NTRES INTEGER NOT NULL,
MIND REAL NOT NULL,
MAXD REAL NOT NULL,
MEAND REAL NOT NULL,
SCORING_FUNCTION TEXT NOT NULL,
SCORE_MAX REAL NOT NULL,
EPS REAL NOT NULL,
MIN_SAMPLES INTEGER NOT NULL,
N_CLUSTERS INTEGER NOT NULL,
P_NOISE REAL NOT NULL,
FOREIGN KEY(IDPDB) REFERENCES Structures(IDPDB))""")
        except sqlite3.DatabaseError as err:
            print("Error: ", err)
            conn.close()
            sys.exit()
        else:
            conn.commit()
    return cursor, conn


def download(filelist, q, lock, cursor, conn, dir):
    """

    :param filelist:
    :param q:
    :param lock:
    :param cursor:
    :param conn:
    :param dir:
    """
    for file in filelist:
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(file, pdir=os.path.join(dir, file), file_format='pdb')
        if not os.path.exists(os.path.join(dir, file, 'pdb{:s}.ent'.format(file))):
            print("File with ID PDB: {:s} not found!".format(file))
            continue
        parser = PDBParser()
        structure = parser.get_structure('{:s}', os.path.join(dir, file, 'pdb{:s}.ent'.format(file)))
        name = parser.header.get('name', '')
        head = parser.header.get('head', '')
        method = parser.header.get('structure_method', '')
        res = parser.header.get('resolution', '')
        ncomp = 0
        nchain = 0
        eclist = []
        for values in parser.header['compound'].values():
            ncomp += 1
            nchain += len(values['chain'].split(','))
            eclist.append(values.get('ec', '') or values.get('ec_number', ''))
        ec = ", ".join(eclist)
        nres = 0
        mmass = 0
        ppb = PPBuilder()
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()
            nres += len(seq)
            seqan = ProteinAnalysis(str(seq))
            mmass += int(seqan.molecular_weight())
        lock.acquire()
        try:
            cursor.execute("""INSERT INTO Structures (IDPDB, NAME, HEAD, METHOD, RESOLUTION, NCOMP, NCHAIN,
NRES, MMASS, EC) VALUES ("{:s}", "{:s}", "{:s}", "{:s}", {:.2f}, {:d}, {:d},{:d}, {:d}, "{:s}")""".format(
                file, name, head, method, res, ncomp, nchain, nres, mmass, ec))
        except sqlite3.DatabaseError as err:
            print("Error: ", err)
            continue
        else:
            print("Download Done for ID PDB: {:s}".format(file))
            conn.commit()
            q.put(file)
        finally:
            lock.release()
    q.put(None)


def graph(cls, dir, basefile):
    """

    :param dir:
    :param cls:
    :param basefile:
    :return:
    """
    grid, legend = True, True
    sa = os.path.join(dir, '{:s}'.format(basefile + '.png'))
    try:
        fig, ax = cls.graph(grid, legend)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError:
        print('Error! Failed to plot!\n')


def colormap(cls, newdir: str, basefile: str):
    """

    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    sa = os.path.join(newdir, '{:s}'.format(basefile + '.cm.png'))
    try:
        fig = cls.colormap(grid_state=True)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError:
        print('Error! Failed to plot!\n')


def save_state(cls, newdir: str, basefile: str):
    """

    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    st = os.path.join(newdir, '{:s}'.format(basefile + '.dat'))
    try:
        cls.savestate(st)
    except FileNotFoundError:
        return


def save_pymol(cls, newdir: str, basefile: str):
    """

    :param cls:
    :param newdir:
    :param basefile:
    :return:
    """
    pymol = os.path.join(newdir, '{:s}'.format(basefile + '.py'))
    try:
        cls.save_pymol_script(pymol)
    except FileNotFoundError:
        return


def db_save(con, curr, lock, file, htable, ntres, mind, maxd,
            meand, metric, score, eps, min_samples, n_clusters, p_noise):
    """

    :param p_noise:
    :param con:
    :param curr:
    :param lock:
    :param file:
    :param htable:
    :param ntres:
    :param mind:
    :param maxd:
    :param meand:
    :param metric:
    :param score:
    :param eps:
    :param min_samples:
    :param n_clusters:
    :return:
    """
    lock.acquire()
    try:
        curr.execute("""INSERT INTO Results (IDPDB, HTABLE, NTRES, MIND, MAXD, MEAND, SCORING_FUNCTION, SCORE_MAX, EPS,
MIN_SAMPLES, N_CLUSTERS, P_NOISE) VALUES ("{:s}", "{:s}", {:d}, {:.2f}, {:.2f}, {:.2f}, "{:s}", {:.3f}, {:.2f}, {:d}, {:d}, {:.2f})""".format(
            file, htable, ntres, mind, maxd, meand, metric, score, eps, min_samples, n_clusters, p_noise))
    except sqlite3.DatabaseError as err:
        print("Error: ", err)
        return
    else:
        con.commit()
    finally:
        lock.release()


def clusterThread(file, dir, cursor, conn, lock, min_eps, max_eps, step_eps,
                  min_min_samples, max_min_samples, n_jobs=1):
    """

    :param n_jobs:
    :param file:
    :param dir:
    :param cursor:
    :param conn:
    :param lock:
    :param min_eps:
    :param max_eps:
    :param step_eps:
    :param min_min_samples:
    :param max_min_samples:
    :return:
    """
    cls = ClusterPdb()
    cls.open_pdb(os.path.join(dir, file, 'pdb{:s}.ent'.format(file)))
    for htable in htables:
        try:
            ntres, mind, maxd, meand = cls.parser(htable=htable, pH=pH)
        except ValueError:
            print('Error! Invalid file format\nor file does not contain {:s} residues\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'nanodroplet')
                else 'negative' if htable == 'negative' else 'positive'))
            return
        dir_ptable = os.path.join(dir, file, htable)
        try:
            os.makedirs(dir_ptable, exist_ok=True)
        except OSError:
            print('Unable to create folder ' + dir_ptable)
            continue
        for metric in metrics:
            try:
                cls.init_cycles(min_eps, max_eps, step_eps, min_min_samples, max_min_samples,
                                n_jobs=n_jobs, metric=metric)
                for _ in cls.auto_yield():
                    pass
                eps, min_samples = cls.auto()
            except ValueError:
                print('Error! Could not parse file or clustering failed\n')
                continue
            else:
                print("Job completed for ID_PDB: {:s}, ptable: {:s}, metric: {:s}".format(file, htable, metric))
                db_save(conn, cursor, lock, file, htable, ntres, mind, maxd, meand, metric,
                        cls.score, eps, min_samples, cls.n_clusters, cls.noise_percent())
                dir_metric = os.path.join(dir_ptable, metric)
                try:
                    os.makedirs(dir_metric, exist_ok=True)
                except OSError:
                    print('Unable to create folder ' + dir_metric)
                    continue
                save_pymol(cls, dir_metric, file)
                graph(cls, dir_metric, file)
                try:
                    colormap(cls, dir_metric, file)
                    save_state(cls, dir_metric, file)
                except ValueError as err:
                    print(err)


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
    dTask = Process(target=download, args=(filelist, q, lock, cursor, conn, outputDir))
    dTask.start()
    TaskDone = False
    clusterTasks = []
    n = 1
    tmp_dict = {}
    while not TaskDone:
        while len(clusterTasks) < psutil.cpu_count() - 1:
            item = q.get()
            if item is None:
                TaskDone = True
                break
            p = Process(target=clusterThread, args=(item, outputDir, cursor, conn, lock, min_eps, max_eps, step_eps,
                                                    min_min_samples, max_min_samples))
            p.start()
            tmp_dict[p.pid] = item
            clusterTasks.append(p)
        for process in clusterTasks:
            if not process.is_alive():
                clusterTasks.remove(process)
                print("All tasks completed for {:s} ({:d}/{:d})".format(tmp_dict[process.pid], n, len(filelist)))
                del tmp_dict[process.pid]
                n += 1
    for p in clusterTasks:
        p.join()
    dTask.join()
    conn.close()
