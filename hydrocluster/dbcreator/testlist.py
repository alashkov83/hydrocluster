#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by lashkov on 01.06.18"""
import os
import os.path
import sqlite3
import sys
import warnings
from queue import Queue
from threading import Lock, Thread

import hjson
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from matplotlib.backends.backend_agg import FigureCanvasAgg

try:
    from ..core.pdbcluster import ClusterPdb, n_usecpus
except ImportError:
    print('Error! Scikit-learn not installed!')
    sys.exit()

warnings.filterwarnings("ignore")


def read_config(input_fn: str) -> dict:
    """
    :param input_fn:
    :return:
    """
    valid_keys = ('Input file name', 'Project name', 'Property tables list', 'pH', 'Minimum eps value (A)',
                  'Maximum eps value (A)', 'Step of eps value (A)', 'Minimum min_samples', 'Maximum min_samples',
                  'Scoring coefficients list', 'Save states')
    ptables = ('hydropathy', 'nanodroplet', 'menv', 'rekkergroup', 'fuzzyoildrop', 'aliphatic_core', 'hydrophilic',
               'positive', 'negative', 'ngroup', 'pgroup')
    scores = ('si_score', 'calinski', 's_dbw', 'cdbw')
    if not input_fn:
        raise ValueError("Input file name is not defined!")
    if not os.path.exists(input_fn):
        raise ValueError("Input file was not found!")
    with open(input_fn) as f:
        try:
            options_dict = hjson.load(f)
        except hjson.HjsonDecodeError:
            raise ValueError('Invalid hjson file!')
    if set(options_dict.keys()).issubset(set(valid_keys)):
        print("Parameters in config file:")
        for key in options_dict:
            print("Parameter: {:s}, value: {:s}".format(key, str(options_dict[key])))
    else:
        raise ValueError("Error! Parameter(s): {:s} is not include in supported parameters: {:s}!\n".format(
            ', '.join(set(options_dict.keys()).difference(set(valid_keys))),
            ','.join(valid_keys)))
    if 'Input file name' not in options_dict:
        raise ValueError('Filename of ID PDBs list is not defined!')
    if not os.path.exists(options_dict['Input file name']):
        raise ValueError('File {:s} contains ID PDBs list is unavailable!'.format(options_dict['Input file name']))
    if 'Project name' not in options_dict or not options_dict['Project name']:
        pname = os.path.splitext(os.path.basename(options_dict['Input file name']))[0]
        print('Project name is not defined! Set to {:s}'.format(pname))
        options_dict['Project name'] = pname
    if 'Property tables list' not in options_dict or not options_dict['Property tables list']:
        print('Property tables are not defined! Set to all possible tables')
        options_dict['Property tables list'] = ptables
    if set(options_dict['Property tables list']).issubset(set(ptables)):
        print("Selected ptables: {:s}\n".format(', '.join(ptables)))
    else:
        raise ValueError("Error! Table(s): {:s} is not include in supported tables: {:s}!\n".format(
            ', '.join(set(options_dict['Property tables list']).difference(set(ptables))),
            ','.join(ptables)))
    if 'pH' not in options_dict:
        print('pH is not defined! Set to 7.0')
        options_dict['pH'] = 7.0
    if (type(options_dict['pH']) != int and type(options_dict['pH']) != float) or not (0 < options_dict['pH'] < 14):
        raise ValueError('Invalid pH value')
    else:
        options_dict['pH'] = float(options_dict['pH'])
    if 'Minimum eps value (A)' not in options_dict:
        print('Minimum eps value is not defined! Set to 3.0 \u212B')
        options_dict['Minimum eps value (A)'] = 3.0
    if (type(options_dict['Minimum eps value (A)']) != int and type(options_dict['Minimum eps value (A)']) != float) \
            or options_dict['Minimum eps value (A)'] <= 0:
        raise ValueError('Invalid minimum eps value')
    else:
        options_dict['Minimum eps value (A)'] = float(options_dict['Minimum eps value (A)'])
    if 'Maximum eps value (A)' not in options_dict:
        print('Maximum eps value is not defined! Set to 15.0 \u212B')
        options_dict['Maximum eps value (A)'] = 15.0
    if (type(options_dict['Maximum eps value (A)']) != int and type(options_dict['Maximum eps value (A)']) != float) \
            or options_dict['Maximum eps value (A)'] <= 0:
        raise ValueError('Invalid maximum eps value')
    else:
        options_dict['Maximum eps value (A)'] = float(options_dict['Maximum eps value (A)'])
    if 'Step of eps value (A)' not in options_dict:
        print('Step of eps value is not defined! Set to 0.1 \u212B')
        options_dict['Step of eps value (A)'] = 0.1
    if (type(options_dict['Step of eps value (A)']) != int and type(options_dict['Step of eps value (A)']) != float) \
            or options_dict['Step of eps value (A)'] <= 0:
        raise ValueError('Invalid step of eps value')
    else:
        options_dict['Step of eps value (A)'] = float(options_dict['Step of eps value (A)'])
    if 'Minimum min_samples' not in options_dict:
        print('Minimum min_samples is not defined! Set to 3')
        options_dict['Minimum min_samples'] = 3
    if type(options_dict['Minimum min_samples']) != int or options_dict['Minimum min_samples'] < 2:
        raise ValueError('Invalid minimum min_samples eps value')
    if 'Maximum min_samples' not in options_dict:
        print('Maximum min_samples is not defined! Set to 50')
        options_dict['Maximum min_samples'] = 50
    if type(options_dict['Maximum min_samples']) != int or options_dict['Maximum min_samples'] < 2:
        raise ValueError('Invalid maximum min_samples eps value')
    if 'Scoring coefficients list' not in options_dict or not options_dict['Scoring coefficients list']:
        print('Scoring coefficients are not defined! Set to all scoring coefficients')
        options_dict['Scoring coefficients list'] = scores
    if set(options_dict['Scoring coefficients list']).issubset(set(scores)):
        print("Selected scoring coefficients: {:s}\n".format(', '.join(scores)))
    else:
        raise ValueError("Error! Scoring(s): {:s} is not include in supported scorings: {:s}\n".format(
            ', '.join(set(options_dict['Scoring coefficients list']).difference(set(scores))),
            ','.join(scores)))
    if 'Save states' not in options_dict:
        print('Save states option is not defined! Set to False')
        options_dict['Save states'] = False
    return options_dict


def opendb(fiename: str) -> tuple:
    """

    :param fiename:
    :return:
    """
    isDbExists = os.path.exists("{:s}.db".format(fiename))
    conn = sqlite3.connect("{:s}.db".format(fiename), check_same_thread=False)
    cursor = conn.cursor()
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
CL REAL NOT NULL,
CR REAL NOT NULL,
R2L REAL NOT NULL,
R2R REAL NOT NULL, 
FOREIGN KEY(IDPDB) REFERENCES Structures(IDPDB))""")
        except sqlite3.DatabaseError as err:
            print("Error: ", err)
            conn.close()
            sys.exit()
        else:
            conn.commit()
    return cursor, conn


def download(filelist: list, q: Queue, lock: Lock, cursor: sqlite3.Cursor, conn: sqlite3.Connection, dir_name: str):
    """
    :param filelist:
    :param q:
    :param lock:
    :param cursor:
    :param conn:
    :param dir_name:
    """
    with open('status_tmp.txt', 'w') as f:
        f.write('')
    for file in filelist:
        if file in open('status_tmp.txt').readlines():
            continue
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(file, pdir=os.path.join(dir_name, file), file_format='pdb')
        if not os.path.exists(os.path.join(dir_name, file, 'pdb{:s}.ent'.format(file))):
            print("File with ID PDB: {:s} not found!".format(file))
            continue
        parser = PDBParser()
        structure = parser.get_structure('{:s}', os.path.join(dir_name, file, 'pdb{:s}.ent'.format(file)))
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
            with open('status_tmp.txt', 'at') as f:
                f.write((file + '\n'))
    os.remove('status_tmp.txt')
    q.put(None)


def graph(cls: ClusterPdb, dir_name: str, basefile):
    """

    :param dir_name:
    :param cls:
    :param basefile:
    :return:
    """
    grid, legend = True, True
    sa = os.path.join(dir_name, '{:s}'.format(basefile + '.png'))
    try:
        fig, ax = cls.graph(grid, legend)
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(sa)
    except AttributeError:
        print('Error! Failed to plot!\n')


def colormap(cls: ClusterPdb, newdir: str, basefile: str):
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


def save_state(cls: ClusterPdb, newdir: str, basefile: str):
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


def save_pymol(cls: ClusterPdb, newdir: str, basefile: str):
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


def db_save(con: sqlite3.Connection, curr: sqlite3.Cursor, lock: Lock, file: str, htable: str, ntres: int,
            mind: float, maxd: float, meand: float, metric: str, score: float, eps: float, min_samples: int,
            n_clusters: int, p_noise: float, cl: float, cr: float, r2l: float, r2r: float):
    """

    :param r2r:
    :param cr:
    :param r2l:
    :param cl:
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
MIN_SAMPLES, N_CLUSTERS, P_NOISE, CL, CR, R2L, R2R) VALUES ("{:s}", "{:s}", {:d}, {:.2f}, {:.2f}, {:.2f}, "{:s}", {:.3f}, {:.2f}, {:d}, {:d}, {:.2f}, {:.4f},{:.4f},{:.2f},{:.2f})""".format(
            file, htable, ntres, mind, maxd, meand,
            metric, score, eps, min_samples, n_clusters, p_noise,
            cl, cr, r2l, r2r))
    except sqlite3.DatabaseError as err:
        print("Error: ", err)
        return
    else:
        con.commit()
    finally:
        lock.release()


def clusterThread(file: str, dir_name: str, cursor: sqlite3.Cursor, conn: sqlite3.Connection, lock: Lock,
                  min_eps: float, max_eps: float, step_eps: float, min_min_samples: int, max_min_samples: int,
                  htables: list, metrics: list, pH: float, ss: bool):
    """

    :param ss:
    :param pH:
    :param metrics:
    :param htables:
    :param file:
    :param dir_name:
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
    cls.open_pdb(os.path.join(dir_name, file, 'pdb{:s}.ent'.format(file)))
    for htable in htables:
        try:
            ntres, mind, maxd, meand = cls.parser(htable=htable, pH=pH)
        except ValueError:
            print('\nError! Invalid file format\nor file does not {:s} contain {:s}\n'.format(
                'hydrophobic' if htable in ('hydropathy', 'menv', 'nanodroplet', 'fuzzyoildrop', 'rekkergroup')
                else 'negative' if htable in ('ngroup', 'negative')
                else 'positive' if htable in ('pgroup', 'positive')
                else 'aliphatic', 'groups' if htable in ('rekkergroup', 'pgroup', 'ngroup') else 'residues'))
            return
        dir_ptable = os.path.join(dir_name, file, htable)
        try:
            os.makedirs(dir_ptable, exist_ok=True)
        except OSError:
            print('Unable to create folder ' + dir_ptable)
            continue
        for metric in metrics:
            try:
                cls.init_cycles_old(min_eps, max_eps, step_eps, min_min_samples, max_min_samples, metric=metric)
                for _ in cls.auto_yield():
                    pass
                eps, min_samples = cls.auto()
            except ValueError as e:
                print('Error! Could not parse file or clustering failed for ID_PDB: {:s}, ptable: {:s}, metric: {:s}'
                      '\nError: {:s}'
                      .format(file, htable, metric, str(e)))
                continue
            else:
                cl, cr, r2l, r2r = cls.get_conc
                print("Subjob completed for ID_PDB: {:s}, ptable: {:s}, metric: {:s}".format(file, htable, metric))
                db_save(conn, cursor, lock, file, htable, ntres, mind, maxd, meand, metric, cls.score, eps, min_samples,
                        cls.n_clusters, cls.noise_percent(), cl, cr, r2l, r2r)
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
                    if ss:
                        save_state(cls, dir_metric, file)
                except ValueError as err:
                    print(err)


def db_main(namespace):
    """

    :param namespace:
    """
    inp_fn = namespace.input
    try:
        options_dict = read_config(inp_fn)
    except ValueError as e:
        print(str(e))
        sys.exit(-1)
    inp = options_dict['Input file name']
    output = options_dict['Project name']
    outputDir = '{:s}_data'.format(output)
    htables = options_dict['Property tables list']
    metrics = options_dict['Scoring coefficients list']
    ss = options_dict['Save states']
    pH = options_dict['pH']
    min_eps = options_dict['Minimum eps value (A)']
    max_eps = options_dict['Maximum eps value (A)']
    step_eps = options_dict['Step of eps value (A)']
    min_min_samples = options_dict['Minimum min_samples']
    max_min_samples = options_dict['Maximum min_samples']
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
        clusterThread(item, outputDir, cursor, conn, lock, min_eps, max_eps, step_eps, min_min_samples, max_min_samples,
                      htables, metrics, pH, ss)
        print("All tasks completed for {:s} ({:d}/{:d})".format(item, n, len(filelist)))
        n += 1
    dTask.join()
    conn.close()
