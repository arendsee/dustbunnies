#! /usr/bin/python3

import argparse
import sys
import sqltools.misctools as misc

def parse(parser):

    parser.add_argument(
        '-d', '--delimiter',
        help="CSV output delimiter (default=',')",
        type=str,
        default=',')

    sub = parser.add_subparsers(
        help="add_subparsers STUB",
        title="Retrieval commands",
        dest="retrieval_function")

    raw = sub.add_parser(
        'raw',
        help="Submit raw SQLite command")
    raw.add_argument('sqlcmd')

    mat = sub.add_parser(
        'mat',
        help="Fetch a matrix")
    mat.add_argument(
        '-f', '--filling',
        help=
            """
            Fill matrix with this hsp value. The chosen value must
            correspond to a one of the Hsp table fields.
            """,
        default="bit_score")
    mat.add_argument(
        '-c', '--collection',
        help="Blast collection")

    parser.set_defaults(func=__dispatch__)

def __dispatch__(args):
    call = {'raw': __fetch_and_print__,
            'mat': __score_matrix__}
    call[args.retrieval_function](args)


def __fetch_and_print__(args):
    rows = misc.fetch(args.sqlcmd, args.sqldb)
    for row in rows:
        print((args.delimiter).join(map(str, row)))

def __score_matrix__(args):
    print("stub")
#    sqlcmd = 
#        """
#        
#        """



