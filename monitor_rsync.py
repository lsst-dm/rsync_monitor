#!/usr/bin/env python

# special non-LSST stack reqs:
# configargparse(PyPi, conda)

# Assumptions:
# ** Multiple places assume that only one instance of this code is running on same
#    camera, src_site, and dayobs at same time
# ** Within a particulr camera, src_site and dayobs, filenames are unique without
#    including paths
# ** really simple dump of database tables to html tables

import configargparse 

from datetime import datetime, timedelta
import time
import pprint
#import sqlalchemy
import sqlite3
import sys
import os
import glob
import re
import statistics
import traceback


def get_config():
    """Parse command line args, config file, environment variables to get configuration.
    Returns
    -------
    config: `dict`
        Configuration values for program.
    """
    parser = configargparse.ArgParser(config_file_parser_class=configargparse.YAMLConfigFileParser,
                                      auto_env_var_prefix="RSYMON_")
    parser.add('-c', '--config', action='store', type=str, required=False,
               is_config_file=True,
               help='Config file with paths, etc')
    parser.add('--camera', action='store', type=str, required=True,
               help='camera name')
    parser.add('--src_site', action='store', type=str, required=True,
               help='source site (e.g., Tucson, SLAC)')
    parser.add('--src_log', action='store', type=str, required=True,
               help='path at ncsa to log files containing src site info')
    parser.add('--stor_path_split', action='store', type=str, required=True,
               help='string over which to split full path to get relative path')
    parser.add('--ncsa_storage_prefix', action='store', type=str, required=True,
               help='path to storage directory at NCSA')
    parser.add('--ncsa_gen2_registry', action='store', type=str, required=True,
               help='connection string for NCSA gen2 registry')
    parser.add('--ncsa_gen2_prefix', action='store', type=str, required=True,
               help='path to images in gen2 repo at NCSA')
    parser.add('--ncsa_rsync_log', action='store', type=str, required=True,
               help='glob pattern for rsync log files inside log prefix')
    parser.add('--monitordb', action='store', type=str, 
               help='connection string for database with monitor tables')
    parser.add('--storage_filename_pattern', action='store', type=str,
                help='filename pattern for creating storage filenames')
    parser.add('--gen2_filename_pattern', action='store', type=str,
                help='filename pattern for creating gen2 repo filename')
    parser.add('--write_html', action='store_true', 
               help='Whether to output html for specified dates')
    parser.add('--html_prefix', action='store', type=str, required=True,
               help='path to top level directory for html report')
    parser.add('--init_db', action='store_true', 
               help='initialize DB schema')
    parser.add('--dry_run', action='store_true', 
               help='gather information, but do not save to monitor database or write html')
    parser.add('--debug', action='store_true', 
               help='output very verbose details')
    parser.add('--db_root', action='store', type=str, required=False, 
               help='Base path for sqlite3 files')
    parser.add('--html_root', action='store', type=str, required=False, 
               help='Base path for html files')
    parser.add('dayobs', action='store', type=str, nargs='+', 
               help='YYYY-MM-DD [YYYY-MM-DD] (treated as inclusive date range if 2 dates given)')
    config = vars(parser.parse_args())
    if len(config['dayobs']) > 2:
        parser.error('dayobs must be 1 date or date range (2 dates)')
    return config


def init_db(dbconn):
    """Initialize rsync monitor database tables
    Assumes tables do not already exist
    Parameters
    ----------
    dbconn: `str`
        Database connection string
    """
    try:
        # since sqlite is a file, make directory
        os.makedirs(os.path.dirname(dbconn), exist_ok=True)
        conn = sqlite3.connect(dbconn)
    except sqlite3.OperationalError:
        print("Error connecting to '%s'" % dbconn)
        raise
    c = conn.cursor()

    c.execute('''create table monitor_unexpected_file (
                camera text, 
                src_site text, 
                dayobs text, 
                loctype text, 
                path text, 
                filename text)''')

    c.execute('''create table monitor_note (
                camera text, 
                src_site text, 
                start_dayobs text, 
                end_dayobs text, 
                last_update text, 
                notes text)''')

    c.execute('''create table monitor_transfer (
                transfer_id integer primary key autoincrement, 
                camera text, 
                src_site text, 
                dayobs text, 
                log_filename text, 
                batch_size integer, 
                total_bytes_sent integer, 
                total_file_bytes integer, 
                bytes_per_sec real, 
                start_epoch integer, 
                end_epoch integer)''')

    c.execute('''create table monitor_rsync_file (
                dayobs text, 
                camera text, 
                src_site text, 
                filename text, 
                relpath text, 
                time_obs text,
                time_src text, 
                time_ncsa text, 
                id_gen2 integer, 
                gen2_ingest_error text, 
                time_gen2 text)''')

    c.execute('''create table monitor_rsync_summary (
                camera text, 
                src_site text, 
                dayobs text, 
                cnt_obs integer, 
                cnt_src integer, 
                cnt_ncsa integer, 
                cnt_exp_gen2 integer,
                cnt_gen2 integer, 
                cnt_gen2_error integer, 
                time_xfer_min integer, 
                time_xfer_max integer, 
                time_xfer_mean real, 
                time_xfer_stdev real, 
                time_gen2_min integer, 
                time_gen2_max integer, 
                time_gen2_mean real, 
                time_gen2_stdev real, 
                last_update text default CURRENT_TIMESTAMP)''')
    conn.commit()
    conn.close()


def unexpected_file(step, config, filename, finfo):
    """Gather information about unexpected file for later saving to database.

    Parameters
    ----------
    step: `str`
        Name of the step in file lifecycle to help in debugging.
    config: `dict`
        Program configuration values.
    filename: `str`
        Name of the unexpected file.
    path: `str`
        Path to the unexpected file.

    Returns
    -------
    unexp: `dict`
        Dictionary with information about unexpected file.
    """
    unexp = {'camera': config['camera'],
             'src_site': config['src_site'],
             'dayobs': config['curr_dayobs'],
             'path': finfo['path'],
             'filename': filename}

    return unexp


def get_obs_info(config):
    """Get observing information about what images should exist
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    obsInfo: `dict`
        Source of truth about what images should exist.
    """
    # Postponed until have way to query EFD or ETraveler (SLAC) or equiv
    obsInfo = {}
    return obsInfo


def get_src_info(config):
    """Get information about images already written at source site.
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    src_info: `dict`
        Information gathered about images at source site. (key: relpath/filename)
    """
    stime = time.time()
    src_info = {}
    dayobs_year = config['curr_dayobs'][0:4]
    # src_info_file format: filename;size in bytes;modification time in epoch
    src_info_file = config['src_log'].format(dayobs=config['curr_dayobs'],
                                             dayobs_nohyphens=config['curr_dayobs_nohyphens'],
                                             dayobs_year=dayobs_year)
    if os.path.exists(src_info_file):
        with open(src_info_file) as srcfh:
            for line in srcfh:
                lparts = line.strip().split(';')
                fname = os.path.basename(lparts[0])
                fpath = os.path.dirname(lparts[0])
                try:
                    endpath = fpath.split(config['splitstr'], 1)[1]
                except IndexError:
                    endpath = ""
                relpath = os.path.join(config['splitstr'], endpath)
                relfname = os.path.join(relpath, fname)
 
                src_info[relfname] = {'filename': fname,
                                      'path': fpath,
                                      'relpath': relpath,
                                      'size': lparts[1],
                                      'mtime': lparts[2]
                                     }
    etime = time.time()
    print("TIMING: get_src_storage_info %f secs" % (etime-stime))
    return src_info


def get_ncsa_storage_info(config, ncsa_xfer_info):
    """Get information about images already copied to ncsa storage
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    ncsa_storage: `dict`
        Information about images already copied to ncsa storage. (key: relpath/filename)
    """
    stime = time.time()
    ncsa_storage = {}
    storprefix = config['ncsa_storage_prefix'].format(dayobs=config['curr_dayobs'], 
                                                   dayobs_nohyphens=config['curr_dayobs_nohyphens'])

    print(storprefix)
    for dirroot, _, files in os.walk(storprefix, topdown=False):
        print(dirroot)
        for fname in files:
            print(fname)
            try:
                endpath = dirroot.split(config['splitstr'], 1)[1]
            except IndexError:
                endpath = ""
            relpath = os.path.join(config['splitstr'], endpath)
            relfname = os.path.join(relpath, fname)
            finfo = {'filename': fname, 
                    'path': relpath,
                    'relpath': relpath}
            #sstat = os.stat(os.path.join(dirroot, fname))
            #finfo['mtime'] = int(sstat.st_mtime)
            # size could also be stored here 

            ncsa_storage[relfname] = finfo
        
            if relfname in ncsa_xfer_info:
                finfo['mtime'] = ncsa_xfer_info[relfname]['logtime']
            else:
                print("Missing from ncsa_xfer_info: %s" % relfname)
                sstat = os.stat(os.path.join(dirroot, fname))
                finfo['mtime'] = int(sstat.st_mtime)

    etime = time.time()
    print("TIMING: get_ncsa_storage_info %f secs" % (etime-stime))
    return ncsa_storage


def get_ncsa_transfer_info(config):
    """Get information about the image transfer to ncsa.
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    ncsa_xfer_chunk_info: `dict`
        Information (esp. timing) about transfer of files to ncsa.
            (key: log filename)
    ncsa_xfer_file_info: `dict`
        Information (esp. timing) about transfer of files to ncsa.
            (key: image filename)
    """
    stime = time.time()
    ncsa_xfer_chunk_info = {}
    ncsa_xfer_file_info = {}
    rlogs = glob.glob(config['ncsa_rsync_log'].format(dayobs=config['curr_dayobs'], 
                                                      dayobs_nohyphens=config['curr_dayobs_nohyphens']))
    for rlogname in rlogs:
        logstat = os.stat(os.path.join(rlogname))
        cinfo = {'camera': config['camera'],
                 'src_site': config['src_site'],
                 'dayobs': config['curr_dayobs'],
                 'log_filename': os.path.basename(rlogname),
                 'total_bytes_sent': 0,
                 'bytes_per_sec': 0,
                 'batch_size': 0}
        
        with open(rlogname, 'r') as lfh:
            for line in lfh:
                if line.startswith('send  <f'):
                    #send  <f+++++++++ <relpath>/<filename> <bytes> <md5sum>
                    parts = line.split()[2:]
                    path = os.path.dirname(parts[0])
                    try:
                        endpath = path.split(config['splitstr'], 1)[1]
                    except IndexError:
                        endpath = ""
                    relpath = os.path.join(config['splitstr'], endpath)
                    fname = os.path.basename(parts[0])
                    relfname = os.path.join(relpath, fname)

                    finfo = {'log_filename': cinfo['log_filename'],
                             'filename': fname,
                             'path': relpath,
                             'filesize': int(parts[1].replace(',','')),
                             'chksum': parts[2],
                             'logtime': int(logstat.st_mtime)}
                    ncsa_xfer_file_info[relfname] = finfo
                    cinfo['batch_size'] += 1
                elif line.startswith('sent'):
                    #sent 151,196,813 bytes  received 58 bytes  27,490,340.18 bytes/sec
                    m = re.match('sent ([0-9][0-9,.]+) bytes .* ([0-9][0-9,.]+) bytes/sec', line) 
                    if m is not None:
                        cinfo['total_bytes_sent'] = int(m.group(1).replace(',',''))
                        cinfo['bytes_per_sec'] = float(m.group(2).replace(',',''))
                elif line.startswith('total size'):
                    # total size is 755,798,400  speedup is 1.00
                    cinfo['total_file_bytes'] = int(line.split()[3].replace(',',''))
        ncsa_xfer_chunk_info[cinfo['log_filename']] = cinfo

    if config['debug']:
        pp = pprint.PrettyPrinter(indent=4)
        print('ncsa xfer  ==================================================')
        pp.pprint(ncsa_xfer_chunk_info)
        pp.pprint(ncsa_xfer_file_info)
                
    etime = time.time()
    print("TIMING: get_ncsa_transfer_info %f secs" % (etime-stime))
    return ncsa_xfer_chunk_info, ncsa_xfer_file_info


def save_transfer_info(config, xfer_chunk, xfer_file):
    """Save transfer information into monitor database getting
       database ids for later use
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    xfer_chunk: `dict`
        Information (esp. timing) about transfer at transfer batch level
            (key: log filename)
    xfer_file: `dict`
        Information (esp. timing) about transfer at file level
            (key: image filename)
    Returns
    -------
    xfer_file: `dict`
        Updated dictionary with transfer_id from database
    """
    if not config['dry_run']:
        if len(xfer_chunk) > 0:
            mconn = sqlite3.connect(config['monitordb'])
            curs = mconn.cursor()
            curs.execute('delete from monitor_transfer where dayobs="%s" and camera="%s" and src_site="%s"' % 
                            (config['curr_dayobs'], config['camera'], config['src_site']))
            print('\tNumber of rows deleted from monitor_transfer = %d' % curs.rowcount)
            xfer = next(iter(xfer_chunk.values()))
            keys = xfer.keys()
            sql = 'insert into monitor_transfer (%s) values (%s)' % (','.join(keys), 
                                                                 ','.join([':%s' % x for x in keys]))
            for xfer in xfer_chunk.values():
                curs.execute(sql, xfer)
                xfer['transfer_id'] = curs.lastrowid
            mconn.commit()
            mconn.close()
    else:
        # fake database ids (using negative values)
        cnt = -1
        for xfer in xfer_chunk.values():
            xfer['transfer_id'] = cnt
            cnt -= 1

    # update file-level info with transfer chunk ids
    for xfile in xfer_file.values():
        if xfile['log_filename'] in xfer_chunk:
            xfile['transfer_id'] = xfer_chunk[xfile['log_filename']]['transfer_id']
        else:
            raise KeyError('Inconsitency in transfer information.  Cannot find %s entry in xfer_chunk' % xfile['log_filename'])


def get_ncsa_gen2_registry_info(config):
    """Get registry info for images already ingested into ncsa gen2 repo
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    ncsa_gen2: `dict`
        Dictionary containing information about images from gen2 repo
    srmapping: `dict`
        Mapping of storage filename to gen2 repo id
    """
    stime = time.time()
    repoprefix = os.path.join(config['ncsa_gen2_prefix'])
    nohyphens = config['curr_dayobs_nohyphens']

    registry = config['ncsa_gen2_registry']
    dbh = sqlite3.connect(registry)
    curs = dbh.cursor()
    gstime = time.time()
    sql = 'select * from raw where dayobs="%s"' % config['curr_dayobs']
    curs.execute(sql)
    print("\tTiming gen2repo select * from raw %f" % (time.time()-gstime))
    colnames = [d[0].lower() for d in curs.description]
    ncsa_gen2 = {}
    for row in curs:
        finfo = dict(zip(colnames, row))
        finfo['orig_date'] = datetime.strptime(finfo['date'], '%Y-%m-%dT%H:%M:%S.%f')
        # create repo fullname
        finfo['path'] = repoprefix
        finfo['filename'] = config['gen2_filename_pattern'].format(**finfo, dayobs_nohyphens=nohyphens)
        linkfull = os.path.join(finfo['path'], finfo['filename'])
        gs2time = time.time() 
        try:
            finfo['orig_filename'] = os.path.realpath(linkfull)
        except OSError:
            finfo['orig_filename'] = config['storage_filename_pattern'].format(**finfo, dayobs_nohyphens=nohyphens)
        #print("\tTiming gen2repo realpath %f" % (time.time()-gs2time))
        try:
            endpart = finfo['orig_filename'].split(config['splitstr'], 1)[1]
        except IndexError:
            endpart = ""
        finfo['orig_relfname'] = os.path.join(config['splitstr'], endpart)
        gs3time = time.time()
        try:
            stat_results = os.lstat(os.path.join(finfo['path'], finfo['filename']))
            #print(os.path.join(finfo['path'], finfo['filename']))
            #print(os.path.realpath(os.path.join(finfo['path'], finfo['filename'])))
            finfo['mtime'] = int(stat_results.st_mtime)
        except OSError:
            print("Couldn't find repo file: %s" % os.path.join(finfo['path'], finfo['filename']))
            pass
        #print("\tTiming gen2repo lstat %f" % (time.time()-gs3time))
        ncsa_gen2[finfo['id']] = finfo
    
    # create storage filename to gen2 id mapping (assumes unique filename in dayobs) 
    srmapping = {}
    for finfo in ncsa_gen2.values():
        #srmapping[finfo['orig_relfname']] = finfo['id']
        srmapping[finfo['id']] = finfo['orig_relfname']

    if config['debug']:
        pp = pprint.PrettyPrinter(indent=4)
        print('ncsa_gen2  ==================================================')
        pp.pprint(ncsa_gen2)
        pp.pprint(srmapping)

    etime = time.time()
    print("TIMING: get_ncsa_gen2_registry_info %f secs" % (etime-stime))
    return ncsa_gen2, srmapping



def prepare_unexpected_file(location, config, fname, finfo):
    """Create dictionary containing info about unexpected file
    Parameters
    ----------
    location: `str`
        Label for where unexpected file was found (storage, repo, etc)
    config: `dict`
        Program configuration values.
    fname: `str`
        Name of unexpected file
    finfo: `dict`
        Dictionary with information about unexpected file

    Returns
    -------
    unexp: `dict`
        Dictionary containing information about unexpected file
    """
    unexp = {'camera': config['camera'],
             'src_site': config['src_site'],
             'dayobs': config['curr_dayobs'],
             'loctype': location,
             'path': finfo['path'],
             'filename': fname}
    return unexp


def prepare_file_info(config, obs_info, src_info, ncsa_storage, ncsa_gen2, srmapping):
    """Prepare file-level information for saving into database later.
    Parameters
    ----------
    config: `dict`
        Program configuration values.

    Returns
    -------
    monitor_img_info: `dict`
        Dictionary containing file-level information 
    """
    monitor_img_info = {}
    unexp_files = []
    
    # process information from obs/ccs
    #   Placeholder for later source of truth
    
    # process information from src site
    for relfname, finfo in src_info.items():
        if relfname not in monitor_img_info:
            # change when have obsInfo, but for now src site is source of truth
            mdict = {'dayobs': config['curr_dayobs'],
                     'camera': config['camera'],
                     'src_site': config['src_site'],
                     'filename': finfo['filename'],
                     'relpath': finfo['relpath'],
                     'time_obs': None,
                     'time_src': int(finfo['mtime']),
                     'time_ncsa': None,
                     'id_gen2': None,
                     'gen2_ingest_error': None,
                     'time_gen2': None
                    }
            monitor_img_info[os.path.join(relfname)] = mdict
        
    # process information from NCSA storage
    for relfname, finfo in ncsa_storage.items():
        if relfname not in monitor_img_info:
            unexp_files.append(prepare_unexpected_file('ncsa_storage', config, relfname, finfo))
        else:
            monitor_img_info[relfname]['time_ncsa'] = finfo['mtime']
    
    # process information from NCSA gen2 repo
    for id_gen2, finfo in ncsa_gen2.items():
        if finfo['orig_relfname'] not in monitor_img_info:
            print("Not in monitor_img_info: %s" % finfo['orig_relfname'])
            if id_gen2 not in srmapping:
                unexp_files.append(prepare_unexpected_file('ncsa_gen2', config, finfo['orig_filename'], finfo))
            else:
                print("found it by id_gen2 %d" % id_gen2)
        else:
            monitor_img_info[finfo['orig_relfname']]['id_gen2'] = finfo['id']
            monitor_img_info[finfo['orig_relfname']]['time_gen2'] = finfo['mtime']
    
    if config['debug']:
        pp = pprint.PrettyPrinter(indent=4)
        print('monitor_img_info ==================================================')
        pp.pprint(monitor_img_info)
        pp.pprint(unexp_info)
        
    return monitor_img_info, unexp_files

    
def prepare_summary_info(config, monitor_info, unexp_info):
    """Prepare the summary info for ingestion to database.
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    Returns
    -------
    suminfo: `dict`
        Summary of current lifecycle info for images taken on dayobs
    """
    # create summary info
    suminfo = {'camera': config['camera'],
               'src_site': config['src_site'],
               'dayobs': config['curr_dayobs'],
               'cnt_obs': 0,
               'cnt_src': 0,
               'cnt_ncsa': 0,
               'cnt_exp_gen2': 0,
               'cnt_gen2': 0,
               'cnt_gen2_error': 0,
               'time_xfer_min': None,
               'time_xfer_max': None,
               'time_xfer_mean': None,
               'time_xfer_stdev': None,
               'time_gen2_min': None,
               'time_gen2_max': None,
               'time_gen2_mean': None,
               'time_gen2_stdev': None,
              }
    time_xfer = []
    time_gen2 = []
    for fname, finfo in monitor_info.items():
        start_time = None
        if fname.endswith('fits'):
            suminfo['cnt_exp_gen2'] += 1
        if finfo['time_obs'] is not None:
            suminfo['cnt_obs'] += 1 
            start_time = finfo['time_obs']
        if finfo['time_src'] is not None:
            suminfo['cnt_src'] += 1 
            start_time = finfo['time_src']
        if finfo['time_ncsa'] is not None:
            suminfo['cnt_ncsa'] += 1 
            if start_time is not None:
                time_xfer.append(finfo['time_ncsa'] - start_time)
        if finfo['time_gen2'] is not None:
            suminfo['cnt_gen2'] += 1 
            if start_time is not None:
                time_gen2.append(finfo['time_gen2'] - start_time)
    
    if len(time_xfer) > 0:
        suminfo['time_xfer_min'] = min(time_xfer)
        suminfo['time_xfer_max'] = max(time_xfer)
        suminfo['time_xfer_mean'] = round(statistics.mean(time_xfer), 2)
        if len(time_xfer) > 1:
            suminfo['time_xfer_stdev'] = round(statistics.stdev(time_xfer), 2)
        print("why %0.2f" % suminfo['time_xfer_min'])
    
    if len(time_gen2) > 0:
        suminfo['time_gen2_min'] = min(time_gen2)
        suminfo['time_gen2_max'] = max(time_gen2)
        suminfo['time_gen2_mean'] = round(statistics.mean(time_gen2), 2)
        if len(time_gen2) > 1:
            suminfo['time_gen2_stdev'] = round(statistics.stdev(time_gen2), 2)
    
    if config['debug']:
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(time_xfer)
        pp.pprint(time_gen2)
        pp.pprint(suminfo)

    return suminfo


def save_file_info(config, curs, file_level_info):
    """Save file-level information to database.
    Parameters
    ----------
    config: `dict`
        Program configuration values.
    curs: 
        Database cursor
    file_level_info: `dict` of `dict`
        Monitor info for each file
    """
    
    # get info from monitor_rsync table
    curs.execute('select * from monitor_rsync_file where dayobs="%s" and camera="%s" and src_site="%s"' % 
                 (config['curr_dayobs'], config['camera'], config['src_site']))
    colnames = [d[0].lower() for d in curs.description]
    rsyncDBinfo = {}
    for row in curs:
        finfo = dict(zip(colnames, row))
        rsyncDBinfo[finfo['filename']] = finfo
    
    #### TODO add check differences between prev and curr data
    
    # delete info in monitor_rsync_file table
    curs.execute('delete from monitor_rsync_file where dayobs="%s" and camera="%s" and src_site="%s"' % 
              (config['curr_dayobs'], config['camera'], config['src_site']))
    print('\tNumber of rows deleted from monitor_rsync_file = %d' % curs.rowcount)
    
    # insert/update info in monitor_rsync_file table
    for finfo in file_level_info.values():
        keys = finfo.keys()
        sql = 'insert into monitor_rsync_file (%s) values (%s)' % (','.join(keys), ','.join([':%s' % x for x in keys]))
        curs.execute(sql, finfo)


def save_unexpected_files(config, curs, unexpInfo):
    """Save new unexpected file info to database.
    Parameters
    ----------
    curs:
        Database cursor.
    unexpInfo: `list` of `dict`
        Information about unexpected files discovered.
    """
    # delete info in monitor_unexpected_file table
    curs.execute('delete from monitor_unexpected_file where dayobs="%s" and camera="%s" and src_site="%s"' % 
              (config['curr_dayobs'], config['camera'], config['src_site']))
    print('\tNumber of rows deleted from monitor_unexpected_file = %d' % curs.rowcount)
    
    for finfo in unexpInfo:
        keys = finfo.keys()
        sql = 'insert into monitor_unexpected_file (%s) values (%s)' % \
              (','.join(keys), ','.join([':%s' % x for x in keys]))
        curs.execute(sql, finfo)


def save_summary_info(config, curs, sumInfo):
    """Save summary info for day to database
    Parameters
    ----------
    curs:
        Database cursor.
    sumInfo: `dict`
        Summary of current lifecycle info for images taken on dayobs.
    """
    # update info in monitor_rsync_summary table
    #(camera text, src_site text, dayobs text, cnt_obs integer, cnt_src integer, cnt_ncsa integer, cnt_gen2 integer, cnt_gen2_error integer, time_xfer_min integer, time_xfer_max integer, time_xfer_mean integer, time_gen2_min integer, time_gen2_max integer, time_gen2_mean integer, last_update text)

    # delete info in monitor_rsync_summary table
    curs.execute('delete from monitor_rsync_summary where dayobs="%s" and camera="%s" and src_site="%s"' % 
              (config['curr_dayobs'], config['camera'], config['src_site']))
    print('\tNumber of rows deleted from monitor_rsync_summary = %d' % curs.rowcount)

    keys = sumInfo.keys()
    sql = 'insert into monitor_rsync_summary (%s) values (%s)' % \
          (','.join(keys), ','.join([':%s' % x for x in keys]))
    curs.execute(sql, sumInfo)
#


def create_dayobs(orig_dayobs):
    """Create list of dayobs strings from command line arg

    Parameters
    ----------
    orig_dayobs: `list` of strings
        1 or 2 YYY-mm-dd indicating endpoints of dayobs range.
    
    Returns
    -------
    all_dayobs: `list` of strings
        List of dayobs strings (YYYY-mm-dd) for report process to be run over.
    """
    
    if len(orig_dayobs) != 1 and orig_dayobs[0] != orig_dayobs[1]:
        startdate = datetime.strptime(orig_dayobs[0], '%Y-%m-%d')
        enddate = datetime.strptime(orig_dayobs[1], '%Y-%m-%d')
        currdate = startdate
        all_dayobs = []
        while currdate <= enddate:
            all_dayobs.append(currdate.strftime('%Y-%m-%d'))
            currdate += timedelta(days=1)
    else:
        all_dayobs = orig_dayobs

    return all_dayobs
    

def read_table(curs, tname, order_by=None):
    sql = 'select * from %s' % tname
    if order_by is not None:
        sql += ' order by %s' % order_by
    curs.execute(sql)
    desc = [tuple[0] for tuple in curs.description]
    rows = []
    for row in curs:
        d = dict(list(zip(desc, row)))
        rows.append(d)
    return (desc, rows)


def read_data(dbstr):
    data = {}
    mconn = sqlite3.connect(dbstr)
    curs = mconn.cursor()
    
    # monitor_rsync_summary
    data['summary'] = read_table(curs, 'monitor_rsync_summary', 'dayobs DESC')

    # monitor_rsync_file
    data['files'] = read_table(curs, 'monitor_rsync_file')

    # monitor_transfer
    data['xfer'] = read_table(curs, 'monitor_transfer')

    # monitor_unexpected_file
    data['unexp'] = read_table(curs, 'monitor_unexpected_file')

    # monitor_note
    data['notes'] = read_table(curs, 'monitor_note')

    mconn.close()
    return data
    

def write_html_table(hfh, tname, tinfo):
    try:
        # monitor_rsync_file
        hfh.write('<p>\n')
        hfh.write('<table border="1" cellpadding="5">\n')
        hfh.write('<tr><th align="center" colspan=100>%s</th></tr>\n' % tname.upper())
        hfh.write('<tr>')
        for c in tinfo[0]:
            hfh.write('<th>%s</th>\n' % c.replace('_', ' '))
        hfh.write('</tr>\n')

        for r in tinfo[1]:
            hfh.write('<tr>')
            for c in tinfo[0]:
                hfh.write('<td>%s</td>' % r[c])
            hfh.write('</tr>\n')
        hfh.write('</table>\n')
        hfh.write('</p>\n')
    except Exception:
        traceback.print_exc()


def write_dayobs_html(prefix, dayobs, data):
    os.makedirs(prefix, exist_ok=True)
    
    htmlfile  = os.path.join(prefix, "%s.html" % dayobs)
    with open(htmlfile, 'w') as hfh:
        hfh.write('<html>\n')
        hfh.write('<head>\n')
        hfh.write('<title>LSST Teststand Rsync Monitor - %s</title>\n' % dayobs)
        hfh.write('<meta http-equiv="refresh" content="300">\n')
        hfh.write('</head>\n')
        hfh.write('<body>\n')
        hfh.write('<h1>LSST Teststand Rsync Monitor - %s</h1>\n' % dayobs)

        # monitor_rsync_summary
        print("summary: %s" % data['summary'])
        cols = data['summary'].keys()
        rows = [data['summary']]
        write_html_table(hfh, 'Summary', (cols, rows))

        # monitor_rsync_file
        rows = data['files'].values()
        cols = []
        if len(rows) > 0:
            cols = next(iter(rows)).keys()
        write_html_table(hfh, 'Files', (cols, rows))

        # monitor_transfer
        rows = data['xfer'].values()
        cols = []
        if len(rows) > 0:
            cols = next(iter(rows)).keys()
        write_html_table(hfh, 'Transfer', (cols, rows))

        # monitor_unexpected_file
        rows = data['unexp']
        cols = []
        if len(rows) > 0:
            cols = rows[0].keys()
        write_html_table(hfh, 'Unexpected Files', (cols, rows))

        # monitor_note
        # todo grab notes from db that pertain to this dayobs

        hfh.write('</body>\n')
        hfh.write('</html>\n')


def write_summary_html(htmlfile, data):
    os.makedirs(os.path.dirname(htmlfile), exist_ok=True)
    with open(htmlfile, 'w') as hfh:
        hfh.write('<html>\n')
        hfh.write('<head>\n')
        hfh.write('<title>LSST Teststand Rsync Monitor</title>\n')
        hfh.write('<meta http-equiv="refresh" content="300">\n')
        hfh.write('</head>\n')
        hfh.write('<body>\n')
        hfh.write('<h1>LSST Teststand Rsync Monitor</h1>\n')

        # monitor_rsync_summary
        write_html_table(hfh, 'Summary', data['summary'])

        # monitor_note
        write_html_table(hfh, 'Notes', data['notes'])

        hfh.write('</body>\n')
        hfh.write('</html>\n')

def main():
    """Program entry point.  Control process that gathers information for 
       each step in the image lifecycle.
    """
    config = get_config()
    all_dayobs = create_dayobs(config['dayobs']) 
    #curr_time = datetime.now()
    curr_time = time.time()

    config['monitordb'] = config['monitordb'].format(**config)
    config['html_prefix'] = config['html_prefix'].format(**config)

    if config['init_db']:
        init_db(config['monitordb'])

    for dayobs in all_dayobs:
        print("Working on dayobs=%s" % dayobs)
        config['curr_dayobs'] = dayobs
        config['curr_dayobs_nohyphens'] = config['curr_dayobs'].replace('-', '')
        config['splitstr'] = config['stor_path_split'].format(dayobs=config['curr_dayobs'],
                                                              dayobs_nohyphens=config['curr_dayobs_nohyphens'])
        # later joins need the relative path not to start with /
        if not config['splitstr'].endswith('/'):
            config['splitstr'] += '/'

        # get and save transfer info (need database ids for use later)
        ncsa_xfer_chunk, ncsa_xfer_file = get_ncsa_transfer_info(config)
        save_transfer_info(config, ncsa_xfer_chunk, ncsa_xfer_file)

        # gather all info
        obs_info = get_obs_info(config)
        src_info = get_src_info(config)
        ncsa_stor_info = get_ncsa_storage_info(config, ncsa_xfer_file)
        ncsa_gen2_registry, srmapping = get_ncsa_gen2_registry_info(config)
        file_level_info, unexp_info = prepare_file_info(config, obs_info, src_info, 
                                                        ncsa_stor_info, ncsa_gen2_registry,
                                                        srmapping)
        summary_info = prepare_summary_info(config, file_level_info, unexp_info)

        if not config['dry_run']:
            # save remaining info to database as single transaction
            mconn = sqlite3.connect(config['monitordb'])
            curs = mconn.cursor()
            save_file_info(config, curs, file_level_info)
            save_unexpected_files(config, curs, unexp_info)  
            save_summary_info(config, curs, summary_info)
            mconn.commit()
            mconn.close()
            # write dayobs html
            if config['write_html']:
                write_dayobs_html(config['html_prefix'], dayobs, {'summary': summary_info,
                                                   'files': file_level_info,
                                                   'xfer': ncsa_xfer_chunk,
                                                   'unexp': unexp_info})

    if not config['dry_run'] and config['write_html']:
        data = read_data(config['monitordb'])
        htmlfile = os.path.join(config['html_prefix'], 'monitor_rsync.html')
        write_summary_html(htmlfile, data)


if __name__ == '__main__':
    main()
