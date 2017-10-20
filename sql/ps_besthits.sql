-- Extracts all useful info from besthits along with mrca data for each match
.separator ','
.header on
.output 'unmasked_at_50g_2013-10-15.csv'
select 
    qgb,
    qlocus,
    species as hspecies, 
    phylostratum,
    mrca as mrca,
    sciname as mrca_sciname,
    qlen,
    palen,
    malen,
    salen,
    mevalue,
    pscore,
    mscore,
    sscore 
from 
    taxid2name
inner join
    (
    --1----------------------------------------------------------
    select
        *
    from
    mrca 
    inner join 
        (
        --2----------------------------------------------------------
        select 
            *
        from 
            blastdatabase
        inner join
            (
            --3----------------------------------------------------------
            select
                *
            from
                besthits
            where
                qgb in
                (
                --4----------------------------------------------------------
                select
                    qgb
                from
                    (
                    select 
                        qgb, max(pscore)
                    from
                        besthits
                    group by 
                        besthits.qgb
                    )
                --4----------------------------------------------------------
                )
            --3----------------------------------------------------------
            ) as bh
        on
            blastdatabase.database = bh.database
        --2----------------------------------------------------------
        )
        as bhbd 
    on 
        mrca.taxid_1 == bhbd.taxid
        and
        mrca.taxid_2 == bhbd.qtaxon
    --1----------------------------------------------------------
    ) 
    as mbhbd
on
    mbhbd.mrca == taxid2name.taxid
group by
    qgb, database
;

