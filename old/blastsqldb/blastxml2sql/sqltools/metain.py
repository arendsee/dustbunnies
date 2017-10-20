#! /usr/bin/python3

import argparse
import csv
import sqlite3 as sql
import sys
import sqltools.misctools as misc

def parse(parser):
	parser.add_argument(
		'-d', '--delimiter', 
		help="csv delimiter (default=';')",
                default=';')
	parser.set_defaults(func=__metain__)

def __metain__(args):
    with args.input as f:
        reader = csv.DictReader(f, delimiter=args.delimiter)
        filename = misc.set_db_path(args.sqldb)
        with sql.connect(filename) as con:
            cursor = con.cursor()
            for line in reader:
                misc.insert(line, 'BlastDatabase', cursor)
