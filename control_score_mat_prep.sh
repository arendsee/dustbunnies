#! /bin/bash

db="$HOME/db/sqlitedb/control.db"
wsp="$HOME/ohome/workspace/control"
orp="$HOME/ohome/lib/sqldbm/orphandbm.py"
col=(
    '1000' 
    '1000.tantan' 
    '1000.segmasker'
    '1000.reverse'
    '1000.reverse.tantan'
    '1000.reverse.segmasker'
    )
filling='hsp_bit_score'
focal_taxid=3702

echo "Updating database"
$orp update -q $db;

echo "Writing dbinfo file"
$orp dump --mrca $focal_taxid -q "$db" > "$wsp/control_dbinfo.csv";

echo "Writing matrices"
for c in ${col[*]}; do
    echo -e "\t$c";
    $orp query mat -l $filling -c $c -q $db -o "$wsp/$c.$filling.csv"; 
done
