select
    max(qalign)
from 
    (
    select
        qlocus, min(((malen - mgaps) / qlen)) as qalign
    from 
        besthits
    group by qlocus
    ) as q
;



-- find universally conserved proteins
-- score > 200
-- percent aligned > 0.5
select
    qgb
from
    (
    select
        qgb, malen, mgaps, qlen, min(score) as score
    from
        (
        select 
            qgb, phylostratum as ps, malen, mgaps, qlen, max(mscore) as score
        from 
            mrca 
        inner join 
            (
            select 
                *
            from 
                blastdatabase
            natural join
                besthits
            )
            as bhbd 
        on 
            mrca.taxid_1 == bhbd.taxid
            and
            mrca.taxid_2 == bhbd.qtaxon
        group by
            qgb, ps
        )
    group by 
        qgb
    )
where
    (malen - mgaps) / qlen > 0.8
;
