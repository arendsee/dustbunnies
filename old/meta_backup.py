#! /usr/bin/python3

import re
import sys
from collections import defaultdict

import orphanlib.sqlite_interface as misc
import orphanlib.entrez_interface as entrez
import orphanlib.initialize as initialize
from   orphanlib.lineage import Lineage, MRCA

# ==================
# EXPORTED FUNCTIONS
# ==================

def update_dbinfo(cur, deep=False, destroy=False, verbose=False):
    if(not misc.table_exists('blastdatabase', cur) or destroy):
        initialize.init_blastdatabase(cur, verbose)
        deep = True
    if(deep):
        print("Retrieving database names from blastreport (this may take awhile)")
        db = misc.get_fields('Blastoutput_db', 'blastreport', cur, is_distinct=True)
        for d in db:
            misc.insert({'database':d}, 'blastdatabase', cur)


    dbfiles = misc.get_fields('database', 'blastdatabase', cur, is_distinct=True)

    for f in dbfiles:
        ori = misc.get_fields(('species','taxid'), 'blastdatabase',
                               cur, ident='database', value=f)[0]

        if(ori[1] is not None): continue

        if(ori[0] is None):
            name = re.sub("\..*", "", f)
            name = re.sub("_", " ", name)
        else:
            name = ori[0]

        taxid = entrez.sciname2taxid(name)

        if(taxid is None):
            print("No taxid found for '{}' (taxon parsed as '{}')".format(f, name))
            name, taxid = _void_taxid(name)

        misc.update({'species': name, 'taxid': taxid}, 'blastdatabase',
                    ('database', f), cur)

def update_mrca(cur, sync=True, taxids=None, verbose=False):
    if(not misc.table_exists('mrca', cur)):
        initialize.init_mrca(cur, verbose)
    if(not misc.table_exists('Taxid2Name', cur)):
        initialize.init_taxid2name(cur, verbose)

    taxid_in = set()
    if(taxids):
        taxid_in.update(taxids)
    if(sync):
        db_taxids = misc.get_fields('taxid', 'blastdatabase', cur, is_distinct=True)
        taxid_in.update(db_taxids)

    # Retrieve taxonomy xml file from entrez
    lin = entrez.taxid2lineage(taxid_in)

    # Find mrca
    mrca = {} # A dict that holds mrca for pairs of taxids
    for k1 in lin:
        for k2 in lin:
            mrca[(k1, k2)] = MRCA(k1, k2)

    # Update MRCA SQL database |taxid1|taxid2|mrca|phylostratum|
    for key in mrca:
        dic = {'taxid_1':key[0].taxid,
               'taxid_2':key[1].taxid,
               'mrca':mrca[key].taxid,
               'phylostratum':mrca[key].level}
        misc.insert(dic, 'MRCA', cur, replace=True)

    # Update Taxid2Name SQL database |taxid|name|
    def t2n_insert(t, s):
        dic = {'sciname': s, 'taxid': t}
        misc.insert(dic, 'Taxid2Name', cur, replace=True)

    for l in lin:
        t2n_insert(l.taxid, l.sciname)
        for pair in l.lineage:
            t2n_insert(*pair)

def update_besthits(cur, verbose=False):
    if(not misc.table_exists('besthits', cur)):
        initialize.init_besthits(cur, verbose)
    col = (
        ('blastoutput_db', 'database'),
        ('query_gb', 'qgb'),
        ('hit_num', 'hit'),
        ('hsp_num', 'hsp'),
        ('hsp_query_from', 'qfrom'),
        ('hsp_query_to', 'qto'),
        ('hsp_hit_from', 'hfrom'),
        ('hsp_hit_to', 'hto'),
        ('hsp_gaps', 'gaps'),
        ('hsp_positive', 'pos'),
        ('hsp_identity', 'ident'),
        ('hsp_bit_score', 'score')
    )
    selection = ', '.join([x[0] for x in col])
    cond = "where hit_num = 1 and query_gb in ('NP_001005221.2', 'NP_001005277.1', 'NP_689699.2')"
    hspdat_cmd = "select {} from blastreport {}".format(selection, cond)
    result = misc.fetch(hspdat_cmd, cur)

    dat = defaultdict(lambda: defaultdict(list))
    for r in result:
        dat[(r[0], r[1])][r[2]].append({
            'hsp'  :r[3],
            'qfrom':r[4],
            'qto'  :r[5],
            'hfrom':r[6],
            'hto'  :r[7],
            'gaps' :r[8],
            'pos'  :r[9],
            'ident':r[10],
            'score':r[11]
        })

    def print_dict_diagnostic(d):
        for pair, hits in d.items():
            for hit, hsps in hits.items():
                print("{} hit:{}".format(pair, hit))
                for hsp in hsps:
                    print("\t{}".format(hsp))

    # Debugging function that neatly print results
    def diagnostics(path, total):
        for key in path:
            print("{} nhsps: {}".format(key, path[key]['nhsp']))
            for field in ['gaps', 'pos', 'ident', 'score', 'alen']:
                path_val = path[key][field]
                total_val = total[key][field]
                print("\t{}: {} {}".format(field, path_val, total_val))

    best = _calculate_besthits(dat, cutoff=20)

    dic = defaultdict(list)
    for pair in best['path']:
        dic['database'].append(pair[0])
        dic['qgb'].append(pair[1])
        for pair, hits in best['path'].items():
            for field, val in hits.items():
                dic['p' + field].append(val)
        for pair, hits in best['total'].items():
            for field, val in hits.items():
                dic['t' + field].append(val)

    col = sorted(dic.keys())
    rows = tuple(tuple(dic[k][i] for k in col) for i in range(len(dic['qgb'])))

    misc.insertmany(col, rows, 'besthits', cur, replace=True)


# =================
# UTILITY FUNCTIONS
# =================

def _set_taxid(name=None):
    qstr = "Please enter taxid (e.g. 3702): "
    taxid = input(qstr)
    try:
        n = entrez.taxid2sciname(taxid)
        if(taxid is None):
            print("Entrez does not recognize this taxon id")
            _set_taxid(name)
        qstr = "You mean this database is comprised solely of '{}' (y/n)? "
        r = input(qstr.format(n)).lower()
        if('n' in r):
            qstr = "Would you like to leave the taxid null (y/n)? "
            r = input(qstr).lower()
            if('y' in r):
                return(name, None)
            else:
                return(_set_taxid())
        qstr = "Would you like to set the taxon name to '{}': "
        r = input(qstr.format(n)).lower()
        if('y' in r):
            name = n
    except:
        print("Error...")
        return(_set_taxid())
    return(name, taxid)

def _void_taxid(name):
    taxid = None
    print("Entrez does not recognize taxon {}".format(name))
    r = input("Do you want to reset the name and try again (y/n)? ").lower()
    if('y' in r):
        name = input("Then give me a new taxon name: ")
        taxid = entrez.sciname2taxid(name)
        if(taxid is None):
            name, taxid = _void_taxid(name)
    else:
        r = input("Would you like to manually assign a taxon id (y/n)? ").lower()
        if('y' in r):
            name, taxid = _set_taxid()
        else:
            print("The taxid and name will be left NULL in the SQL database")
            name, taxid = (None, None)
    return((name, taxid))

def _calculate_besthits(d, cutoff=20, tosum=['gaps', 'pos', 'ident', 'score']):
    """ \
    Arguments:
        d      - data strucutre from update_besthits method
        cutoff - hsp count above which the approximate score algorithm
            will be used
    """
    # Functions used for calculating optiml paths
    def _avg_val(hsps):
        # An approximated score is produced by multiplying the total
        # alignment length by the average score per aligned character
        # Other values (e.g. positive, identical, gaps) are dealt with
        # similarly

        # This sorting step is essential to the algorithm
        # that calculates the non-overlapping alignment length
        hsps = sorted(hsps, key=lambda x: x['qfrom'])

        # Number of query characters that are aligned against at least
        # one HSPs
        alen = 0

        # Beginning of query alignment
        qfrom = hsps[0]['qfrom']

        # End of query alignment
        qto = hsps[0]['qto']

        # The true alignment length is found be combining all overlapping
        # segments and then summing their lengths
        sum_fields = defaultdict(int)
        for i in range(1, len(hsps)):
            for field in tosum:
                sum_fields[field] += hsps[i][field]
            if(hsps[i]['qfrom'] <= qto):
                qto = max(hsps[i]['qto'], qto)
            else:
                alen += qto - qfrom + 1
                qfrom, qto = hsps[i]['qfrom'], hsps[i]['qto']
        alen += qto - qfrom + 1

        # Summed length of all HSPs
        tlen = sum([x['qto'] - x['qfrom'] + 1 for x in hsps])

        avg_fields = {k:((v / tlen) * alen) for k,v in sum_fields.items()}
        avg_fields['alen'] = alen

        return(avg_fields)

    # Otherwise initialize an algorithm to search all possible paths
    # through the directed acyclic graph of non-overlapping subsequences
    # (HSPs) with weights equal to the bitscores. I will use a very crude
    # bruteforce algorithm. The time could be improved with a bit of branch
    # pruning. It also might be worthwhile to implement this is C and
    # optimize the hell out of it.
    # More importantly, this algorithm ensures a path with no overlaps in
    # query sequence, however it still allows overlaps in hits. This should
    # be remedied.
    def _bestpath(v, end=-1, s=0, p=[]):
        scores = []
        eol = True
        for i in range(len(v)):
            nv = v[0:i] + v[(i+1):]
            if(v[i]['qfrom'] > end):
                eol = False
                scores.append(_bestpath(
                    nv,
                    v[i]['qto'],
                    s + v[i]['score'],
                    p + [v[i]['hsp']]))
        if(eol):
            return((s,p))
        else:
            return(max(scores))

    def _global_val(hsps):
        # If there is only one HSP, simply return the single HSPs relevant
        # values. In biological cases, the number of HSPs very often is 0 or 1.

        path = _bestpath(hsps)
        out = {}
        for field in tosum:
            out[field] = sum(hsps[i-1][field] for i in path[1])
            out['alen'] = sum([hsps[i-1]['qto'] - hsps[i-1]['qfrom'] + 1 for i in path[1]])
        out['nhsp'] = len(hsps)
        return(out)

    # The _global_val function is very slow for large numbers of hsps, so the
    # approximating function, _avg_val, will be called when size is above
    # the given cutoff
    def _path(hsps):
        if(len(hsps) <= cutoff):
            out = _global_val(hsps)
        else:
            out = _avg_val(hsps)
        return(out)

    # Straight sum of all hsps summable values
    def _total(hsps):
        out = {k:sum([x[k] for x in hsps]) for k in tosum}
        out['nhsp'] = len(hsps)
        out['alen'] = sum([x['qto'] - x['qfrom'] + 1 for x in hsps])
        return(out)

    # Loop through the data calling the input function
    def _loop(func):
        out = {}
        for pair, hits in d.items():
            highscore = -1
            for hit, hsps in hits.items():
                if(hsps[0]['score'] == 0):
                    out[pair] = {x:0 for x in tosum + ['alen', 'hit']}
                    break
                if(len(hsps) == 1):
                    result = {x:hsps[0][x] for x in tosum}
                    result['alen'] = hsps[0]['qto'] - hsps[0]['qfrom'] + 1
                else:
                    result = func(hsps)
                if(result['score'] > highscore):
                    out[pair] = result
                    out[pair]['hit'] = hit
                    highscore = result['score']
        return(out)

    return({'total':(_loop(_total)),
            'path':(_loop(_path))})
