#! /usr/bin/python3

import sys
import argparse

import sqltools.blastin    as blastin
import sqltools.initialize as initialize
import sqltools.metain     as metain
import sqltools.retrieval  as retrieval

# Top parser
parser = argparse.ArgumentParser(
    description="Manages io for blast result database",
    epilog='For help on any sub-command, enter: %(prog)s [sub-command] -h')

# Input parent parser
_input_parser = argparse.ArgumentParser(add_help=False)
_input_parser.add_argument(
    '-i', '--input',
    help="Input file (from stdin)",
    default=sys.stdin)

sub = parser.add_subparsers(
    help='sub-command help')

# sqltools.blast parser
blast_parser = sub.add_parser(
    'blast', 
    help="Convert BLAST report XML to SQL db",
    parents=[_input_parser],
    description="stub")
blastin.parse(blast_parser)

# sqltools.metadata parser
meta_parser = sub.add_parser(
    'meta', 
    help="Add meta data to SQL db (BUG: can take input only from pipe)",
    parents=[_input_parser],
    description="Description stub")
metain.parse(meta_parser)

# sqltools.initialize parser
init_parser = sub.add_parser(
    'init', 
    help="Initialize SQL database",
    description="Description stub")
initialize.parse(init_parser)

# sqltools.retrieve parser
retr_parser = sub.add_parser(
    'retrieval',
    help="Retrieve data from database",
    aliases=["ret", "retr", "retrieve"])
retrieval.parse(retr_parser)

# This position ensures the sql database will always be the last positional
# argument (don't change this again for Christ's sake)
parser.add_argument(
    'sqldb',
    help="SQL database name")

# Parse arguments
args = parser.parse_args()

# Pass arguments to proper sub-command
args.func(args)
