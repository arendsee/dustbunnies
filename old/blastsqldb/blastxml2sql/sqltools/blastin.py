#! /usr/bin/python3

import argparse
import sqlite3 as sql
import xml.etree.ElementTree as et
from time import gmtime, strftime
import sqltools.misctools as misc
import sys
import os
import re


######## PUBLIC FUNCTION ###################################

def parse(parser):
	parser.add_argument(
		'-c', '--collection', 
		help="parent blast collection")
	parser.add_argument(
		'-q', '--query_taxid', 
		help="set query taxon id") 
	parser.add_argument(
		'-m', '--db_desc', 
		help="BLAST database description")
	parser.set_defaults(func=__blastin__)



def __blastin__(args):
    tree = et.parse(args.input)
    filename = misc.set_db_path(args.sqldb)
    con = sql.connect(filename)
    with con:
        cursor = con.cursor()
        __parse_blast_xml__(tree.getroot(), cursor, args)



def __clean_tag__(tag):
    tag = re.sub('^.*_', '', tag)
    tag = re.sub('-', '_', tag)
    return(tag)



def __last_insert_rowid__(cur):
    cur.execute("select last_insert_rowid()")
    return(cur.fetchall())



def __initialize_collection__(cur, args):
    if(not args.collection): return
    sqlcmd = ("INSERT OR IGNORE INTO BlastCollection (name) " + 
             "VALUES('" + args.collection + "');")
    try:
        cur.execute(sqlcmd)
    except Exception as e:
        print(e)
        sys.exit(1)



def __parse_root__(root, table, cur, args):
    dat = {}
    dat['date_added'] = strftime("%y-%b-%d %H:%M:%S", gmtime())
    if(args.collection):  dat['parent']      = args.collection
    if(args.db_desc):     dat['db_desc']     = args.db_desc
    if(args.query_taxid): dat['query_taxid'] = args.query_taxid
    for child in root:
        if('reference' in child.tag): continue
        if('query' in child.tag): continue
        if('iterations' in child.tag): continue
        if('param' in child.tag):
            for par in child.findall('./Parameters/'):
                dat[__clean_tag__(par.tag)] = par.text
            continue
        dat[__clean_tag__(child.tag)] = child.text
    misc.insert(dat, table, cur) 


 
def __parse_iteration__(iteration, table, cur, args):
    dat = {}
    dat['parent'] = __last_insert_rowid__(cur)[0][0]
    for child in iteration:
        if('Iteration_hits'    == child.tag): continue
        if('Iteration_message' == child.tag): continue
        if('Iteration_stat'    == child.tag):
            for stat in child.findall('./Statistics/'):
                dat[__clean_tag__(stat.tag)] = stat.text
            continue
        dat[__clean_tag__(child.tag)] = child.text
    misc.insert(dat, table, cur)



def __parse_hit__(hit, table, cur, args):
    dat = {}
    dat['parent'] = __last_insert_rowid__(cur)[0][0]
    for child in hit:
        if("Hit_hsps" == child.tag): continue
        dat[__clean_tag__(child.tag)] = child.text
    misc.insert(dat, table, cur)



def __parse_hsp__(hsp, table, cur, args):
    dat = {}
    dat['parent'] = __last_insert_rowid__(cur)[0][0]
    for child in hsp:
        dat[__clean_tag__(child.tag)] = child.text
    misc.insert(dat, table, cur)



# Descend through BLAST XML report structure
def __parse_blast_xml__(root, cur, args):
    __initialize_collection__(cur, args)
    __parse_root__(root, 'BlastOutput', cur, args)
    for iteration in root.findall("./BlastOutput_iterations/Iteration"):
        __parse_iteration__(iteration, 'Iteration', cur, args)
        for hit in iteration.findall("Iteration_hits/Hit"):
            __parse_hit__(hit, 'Hit', cur, args)
            for hsp in hit.findall("Hit_hsps/Hsp"):
                __parse_hsp__(hsp, 'Hsp', cur, args)
