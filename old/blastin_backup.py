#! /usr/bin/python3

import argparse
import os
import re
import sqlite3 as sql
import sys
import xml.etree.cElementTree as et

import orphanlib.initialize as initialize
import orphanlib.sqlite_interface as misc
import orphanlib.meta as meta


# ==================
# EXPORTED FUNCTIONS
# ==================

def parse(parent, *args, **kwargs):
    parser = parent.add_parser(
        'blast',
        help="Convert BLAST report XML to SQL db",
        parents=args,
        description="stub")
    parser.add_argument(
        '-c', '--collection',
        help="blast collection")
    parser.add_argument(
        '-m', '--db_desc',
        help="BLAST database description")
    parser.set_defaults(func=parse_blast_xml)

def parse_blast_xml(args, cur):
    if(not misc.table_exists('blastreport', cur)):
        initialize.init_blastreport(cur, verbose=False)
    if(not misc.table_exists('blastdatabase', cur)):
        initialize.init_blastdatabase(cur, verbose=False)
    root = Root(cur, args)
    root.write()
    meta.update_dbinfo(cur)
    meta.update_mrca(cur)


# =============================
# UTILITY FUNCTIONS AND CLASSES
# =============================

def _parse_fasta_header(row, header):
    try:
        for match in re.finditer('([^|]+)\|([^|]+)', header):
            for tag in ('locus', 'gi', 'taxon', 'gb'):
                if(match.group(1) == tag and match.group(2) != None):
                    row['Query_' + tag] = match.group(2)
    except:
        print("Cannot parse header {}".format(header), file=sys.stderr)

def _addvalue(dat, event, elem):
    if('end' == event and elem.text is not None and not elem.text.isspace()):
        dat[elem.tag] = elem.text

class Root:
    def __init__(self, cur, args):
        self.cur = cur
        self.dat = {}
        self.con = et.iterparse(args.input, events=('start', 'end'))
        if(args.db_desc):     self.dat['db_desc']    = args.db_desc
        if(args.collection):  self.dat['collection'] = args.collection

    def write(self):
        for event, elem in self.con:
            if('start' == event and 'BlastOutput_iterations' == elem.tag):
                for sub_event, sub_elem in self.con:
                    if('start' == sub_event and 'Iteration'):
                        Iteration(self.dat).write(self.con, self.cur)
                    sub_elem.clear()
            else:
                if('reference' in elem.tag or 'query' in elem.tag):
                    pass
                elif('BlastOutput_db' in elem.tag):
                    base = os.path.basename(elem.text)
                    self.dat[elem.tag] = base
                    if(not misc.entry_exists('blastdatabase', 'database', base, self.cur)):
                        misc.insert({'database': base}, 'blastdatabase', self.cur)
                else:
                    _addvalue(self.dat, event, elem)

class Iteration:
    def __init__(self, dat):
        self.dat = dat
        self.children = []

    def write(self, context, cur):
        for event, elem in context:
            if('end' == event and 'Iteration_query-def' == elem.tag):
                _parse_fasta_header(self.dat, elem.text)
                self.dat[elem.tag] = elem.text
            elif('start' == event and 'Iteration_hits' == elem.tag):
                for subevent, subelem in context:
                    if('start' == subevent and 'Hit' == subelem.tag):
                        self.children.append(Hit().parse(context))
                    elif('end' == subevent and 'Iteration_hits' == subelem.tag):
                        break
            elif('end' == event and 'Iteration' == elem.tag):
                break
            else:
                _addvalue(self.dat, event, elem)
        alldat = {'prior': self.dat}
        if(len(self.children) == 0):
            self._insert(alldat, cur)
        for hit in self.children:
            alldat['hit'] = hit.dat
            for hsp in hit.children:
                alldat['hsp'] = hsp.dat
                self._insert(alldat, cur)

    def _insert(self, dat, cur):
        table = {}
        for key in dat:
            for subkey in dat[key]:
                cleankey = re.sub('-', '_', subkey)
                table[cleankey] = dat[key][subkey]
        misc.insert(table, 'BlastReport', cur, True)

class Hit:
    def __init__(self):
        self.dat = {}
        self.children = []

    def parse(self, context):
        for event, elem in context:
            if('start' == event and 'Hit_hsps' == elem.tag):
                self.children.append(Hsp().parse(context))
            elif('end' == event and elem.tag == 'Hit'):
                break
            else:
                _addvalue(self.dat, event, elem)
        return(self)

class Hsp:
    def __init__(self):
        self.dat = {}

    def parse(self, context):
        for event, elem in context:
            if('end' == event and elem.tag == 'Hsp'):
                break
            else:
                _addvalue(self.dat, event, elem)
        return(self)
