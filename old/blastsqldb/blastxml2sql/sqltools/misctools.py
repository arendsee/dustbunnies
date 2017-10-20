#! /usr/bin/python3

import sqlite3 as sql
import os, os.path, sys, argparse

# Builds an INSERT SQL command from a dict input (columns as keys)
def insert(dic, table, cur):
    row = list(map(str, dic.values()))
    col = list(map(str, dic.keys()))
    sqlcmd = ''.join(("INSERT INTO ", table,
                      " (", ','.join(col), ") ",
                      " VALUES('", "','".join(row), "');"))    
    try:
        cur.execute(sqlcmd)
    except Exception as e:
        print(e)
        sys.exit(1)

# Given a SQL database filename, builds a path to the correct location
# in the package homefolder
def set_db_path(filename):
    head, tail = os.path.split(filename)
    if(head and os.path.isfile(filename)): 
        return(filename)
    else:
        path = os.path.join(os.path.abspath('.'), 'data')
        if(not os.path.exists(path)):
            os.mkdir(path)
        dbfile = os.path.join(path, tail)
        return(dbfile)

# Submit a sql cmd
def fetch(cmd, dbname):
    db = set_db_path(dbname)
    con = sql.connect(db)
    with con:
        cur = con.cursor()
        try:
            cur.execute(cmd)
            result = cur.fetchall()
        except Exception as e:
            print(e)
            sys.exit(1)
    return(result)

