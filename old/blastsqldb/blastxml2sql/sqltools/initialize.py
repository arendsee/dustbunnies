#! /usr/bin/python3

#######################################################3
#
# Initializes the blast sql database, creating or clearing
# as necessary
#
########################################################

import argparse
import sqlite3 as sql
import sqltools.misctools as misc
import sys

def parse(parser):
	parser.add_argument(
		'-b', '--blast', 
		help=("Initialize (or clear) BLAST output tables"),
                default=False, action='store_true')
	parser.add_argument(
		'-d', '--database',
		help='Initialize (or clear) BLAST database tables',
		default=False, action='store_true')
	parser.add_argument(
		'-m', '--mrca',
		help='Initialize (or clear) MRCA tables',
		default=False, action='store_true')
	parser.set_defaults(func=__initialize__)

def __initialize__(args):
    if(not args.blast and not args.database and not args.mrca):
        print("Please select '-b', '-d' or 'm' options")
        sys.exit(0)

    filename = misc.set_db_path(args.sqldb)
    con = sql.connect(filename)
    with con:
        cur = con.cursor()
        cur.execute("pragma foreign_keys = ON")
        if(args.blast):    __init_blast__(cur)
        if(args.database): __init_dbinfo__(cur)
        if(args.mrca):     __init_mrca__(cur)



def __init_blast__(cur):
    COLLECTION_VAL = ','.join((
        "name varchar primary key",
        "desc varchar",
    ))

    BLAST_VAL = ','.join((
            # keys
        "parent varchar",
        "sid   integer primary key autoincrement",
            # Added fields
        "db_desc     varchar",
        "query_taxid int check(query_taxid >= 0)",
            # Automatic fields
        "date_added varchar not null",
            # Main tags
        "program varchar not null",
        "version varchar not null",
        "db      varchar not null",
            # Parameter tags
        "matrix     varchar not null",
        "expect     float   not null check(expect >= 0)",
        "gap_open   tinyint not null check(gap_open >= 0)",
        "gap_extend tinyint not null check(gap_extend >= 0)",
        "filter     varchar not null",
            # Constraints
        """foreign key(parent) references BlastCollection(name) 
            on delete set null 
            on update cascade"""
    ))

    ITERATION_VAL = ','.join((
            # keys
        "parent int",
        "sid    integer primary key autoincrement",
            # Main tags
        "iter_num  int     not null check(iter_num >= 0)",
        "query_ID  varchar not null",
        "query_def varchar not null",
        "query_len int     not null check(query_len >= 0)",
            # Statistic tags
        "db_num    int   not null check(db_num >= 0)",
        "db_len    int   not null check(db_len >= 0)",
        "hsp_len   int   not null check(hsp_len >= 0)",
        "eff_space float not null check(eff_space >= 0)",
        "kappa     float not null check(kappa >= 0)",
        "lambda    float not null check(lambda >= 0)",
        "entropy   float not null check(entropy >= 0)",
            # Constraints
        """foreign key(parent) references BlastOutput(sid) 
            on delete cascade
            on update cascade"""
    ))

    HIT_VAL = ','.join((
            # keys
        "parent int",
        "sid   integer primary key autoincrement",
            # Main tags
        "num       int     not null check(num >= 0)",
        "id        varchar not null",
        "def       varchar not null",
        "accession varchar not null",
        "len       int     not null check(num >= 0)",
            # Constraints
        """foreign key(parent) references Iteration(sid) 
            on delete cascade
            on update cascade"""
    ))

    HSP_VAL = ','.join((
            # keys
        "parent int",
            # Main tags
        "num         int   not null check(num >= 0)", 
        "bit_score   float not null check(bit_score >= 0)",
        "score       float not null check(score >= 0)",
        "evalue      float not null check(evalue >= 0)",
        "query_from  int   not null check(query_from >= 0)",
        "query_to    int   not null check(query_to >= 0)",
        "hit_from    int   not null check(hit_from >= 0)",
        "hit_to      int   not null check(hit_to >= 0)",
        "query_frame int   not null check(query_frame >= 0)",
        "hit_frame   int   not null check(hit_frame >= 0)",
        "identity    int   not null check(identity >= 0)",
        "positive    int   not null check(positive >= 0)",
        "align_len   int   not null check(align_len >= 0)",
        "gaps        int            check(gaps >= 0)",
        "qseq        varchar not null",
        "hseq        varchar not null",
        "midline     varchar not null",
            # Constraints
        """foreign key(parent) references Hit(sid) 
            on delete cascade
            on update cascade"""
    ))

    # The drop commands must be performed in this order to avoid
    # to avoid foreign key errors
    cur.execute("DROP TABLE IF EXISTS Hsp")
    cur.execute("DROP TABLE IF EXISTS Hit")
    cur.execute("DROP TABLE IF EXISTS Iteration")
    cur.execute("DROP TABLE IF EXISTS BlastOutput")
    cur.execute("DROP TABLE IF EXISTS BlastCollection")

    # Again these commands must be performed in the correct order
    cur.execute("CREATE TABLE BlastCollection(" + COLLECTION_VAL + ")")
    cur.execute("CREATE TABLE BlastOutput(" + BLAST_VAL + ")")
    cur.execute("CREATE TABLE Iteration(" + ITERATION_VAL + ")")
    cur.execute("CREATE TABLE Hit(" + HIT_VAL + ")")
    cur.execute("CREATE TABLE Hsp(" + HSP_VAL + ")")



def __init_dbinfo__(cur):
    DATABASE_VAL = ','.join((
        "database varchar primary key",
        "taxid int not null check(taxid >= 0)",
        "species varchar not null"
    ))

    cur.execute("DROP TABLE IF EXISTS BlastDatabase")
    cur.execute("CREATE TABLE BlastDatabase(" + DATABASE_VAL + ")")



def __init_mrca__(cur):
    MRCA_VAL = ','.join((
        "focal_taxid  int not null check(focal_taxid >= 0)",
        "outer_taxid  int not null check(outer_taxid >= 0)",
        "mrca_taxid   int not null check(mrca_taxid >= 0)",
        "phylostratum int not null check(phylostratum >= 0)",
        "mrca_name    varchar not null",
        "primary key(focal_taxid, outer_taxid)"
    ))

    cur.execute("DROP TABLE IF EXISTS MRCA")
    cur.execute("CREATE TABLE MRCA(" + MRCA_VAL + ")")


