-- Get scores and phylostratum for each loci
select 
    qgb, min(phylostratum) as phylostratum, pscore, mscore, sscore 
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
where
    pscore > 100
group by
    qgb
;

-- get counts of each phylostratum
select phylostratum, count(phylostratum) as pscount 
from
    (
    select 
        qgb, min(phylostratum) as phylostratum, pscore, mscore, sscore 
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
    where
        sscore > 100
    group by
        qgb
    )
group by
    phylostratum
;

-- Count ghosts
select count(qgb) from (
    select
        qgb, qlocus, pscore, mscore, sscore
    from
        (
        select 
            qgb, qlocus, max(mscore) as mscore, pscore, sscore
        from
            besthits
        where
            database like 'Arabidopsis_thaliana%'
        and
            mscore < 100
        group by 
            qgb
        )
)
;

-- Find proteins that do not significantly self-match
-- OUPUT: ghost gis
.header on
.separator ,
.output at_ghosts_50g.csv
select
    qgb
from
    (
    select 
        qgb, qlocus, max(mscore) as mscore, pscore, sscore
    from
        besthits
    where
        database like 'Arabidopsis_thaliana%'
    and
        mscore <= 100
    group by 
        qgb
    )
;

-- Writes a csv file containing phylostratum data along with scores
-- Reads from masked database
-- OUTPUT:
-- qlocus|species|phylostratum|mrca|sciname|qlen|palen|malen|salen|pscore|mscore|sscore 
-- LIMITED TO: ghosts
attach 'masked_at_ps.db' as masked;
.separator ','
.header on
.output masked_at.csv
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
                        -- masked.besthits
                    where
                        besthits.database like 'Arabidopsis_thaliana%'
                        -- masked.besthits.database like 'Arabidopsis_thaliana%'
                    and
                        besthits.pscore < 100
                        -- masked.besthits.pscore < 100
                    group by 
                        besthits.qgb
                        -- masked.besthits.qgb
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


-- Get everything from a database
.header on
.separator ,
.output masked_at_g50.csv
select 
    qgb,
    qlocus,
    species as hspecies, 
    phylostratum,
    mrca as mrca,
    sciname as mrca_sciname,
    mevalue,
    qlen,
    palen,
    malen,
    salen,
    pscore,
    mscore,
    sscore 
from 
    taxid2name
inner join
    (
    select
        *
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
    ) 
    as mbhbd
on
    mbhbd.mrca == taxid2name.taxid
group by
    qgb, database
;


-- Get corporeal orphans
.header on
.separator ,
.output strata.csv

select
    count(qgb)
from 
    (
    select 
        qgb, min(phylostratum) as phylostratum 
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
    where
        mscore > 75
    group by
        qgb
    )
where
    phylostratum == 19
;

-- Get corporeal orphans by e-value
.header on
.separator ,
.output masked_orphans_at_g50.csv

select
    qgb
from 
    (
    select 
        qgb, min(phylostratum) as phylostratum 
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
    where
        mevalue < 1E-5
    group by
        qgb
    )
where
    phylostratum == 19
;

-- get ghosts
.header on
.separator ,
.output ghosts_at_g50.csv
select 
    qgb
from
    besthits
where
    database like 'Arabidopsis_thaliana%'
and
    mevalue > 1E-5
;
