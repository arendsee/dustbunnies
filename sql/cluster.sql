-- For a given gb, find all putative paralogs
select
    gb,
    count(Hit_num) as clusters
from
(
    select
        Hit_num,
        gb,
        qlen,
        similarity,
        coverage,
        ave_score,
        score
    from
    (
        select
            Hit_num,
            Query_gb as gb,
            Iteration_query_len as qlen,
            round(((1.0 * Hsp_identity) / (Hsp_align_len - Hsp_gaps)), 3) as similarity,
            round(((1.0 * (Hsp_align_len - Hsp_gaps)) / Iteration_query_len), 3) as coverage,
            round(((1.0 * Hsp_bit_score) / (Hsp_align_len - Hsp_gaps)), 3) as ave_score,
            -- max(Hsp_bit_score) as score
            Hsp_bit_score as score
        from
            blastreport
        where
            Query_gb in (
                'NP_199142.1',
                'NP_200952.1',
                'NP_568494.1'
            )
            and
            BlastOutput_db like 'Arabidopsis_thali%'
            and
            Hsp_evalue < 1e-5
        group by
            Hit_num, gb
    )
    where
        similarity > 0
        and
        coverage > 0.3
        and
        ave_score > 0.4
)
group by
    gb
order by 
    clusters
;

                'NP_176378.4',
                'NP_179885.1',
                'NP_199142.1',
                'NP_200952.1',
                'NP_568494.1',
                'NP_199407.1',
                'NP_199186.1',
                'NP_179990.1',
                'NP_200840.1',
                'NP_192745.2',
                'NP_001117619.1',
                'NP_196694.1',
                'NP_563970.1',
                'NP_563640.1',
                'NP_176642.1',
                'NP_565288.1',
                'NP_181813.1',
                'NP_849955.1',
                'NP_001190595.1',
                'NP_568071.1',
                'NP_176517.1',
                'NP_568139.1',
                'NP_973453.2',
                'NP_565456.2',
                'NP_198901.1',
                'NP_187235.1',
                'NP_851211.1',
                'NP_179281.3',
                'NP_196620.1',
                'NP_974675.1',
                'NP_567633.2',
                'NP_177021.1',
                'NP_001154765.1',
                'NP_001031387.1',
                'NP_200211.1',
                'NP_568924.1',
                'NP_001185259.1',
                'NP_568979.1',
                'NP_001031239.1',
                'NP_187437.1',
                'NP_198529.3',
                'NP_181807.1',
                'NP_566321.1',
                'NP_200946.1',
                'NP_179943.1',
                'NP_001190794.1',
                'NP_567062.1',
                'NP_001117605.1',
                'NP_175894.4',
                'NP_974076.1',
                'NP_197822.1',
                'NP_001077894.1',
                'NP_182244.2',
                'NP_173961.3',
                'NP_565098.1',
                'NP_176088.2',
                'NP_176925.1',
                'NP_568532.1',
                'NP_567724.1',
                'NP_850717.1',
                'NP_178532.2',
                'NP_850465.1',
                'NP_191255.1',
                'NP_199358.1',
                'NP_564611.1',
                'NP_001190079.1',
                'NP_564798.1',
                'NP_187156.1',
                'NP_195585.1',
                'NP_568213.2',
                'NP_177142.1',
                'NP_564070.1',
                'NP_200650.1',
                'NP_180201.1',
                'NP_850089.1',
                'NP_178080.2',
                'NP_565557.1',
                'NP_194326.2',
                'NP_181082.1',
                'NP_188503.2',
                'NP_001077451.1',
                'NP_192894.1',
                'NP_001190685.1',
                'NP_181358.1',
                'NP_567394.1',
                'NP_001154582.1',
                'NP_180641.1',
                'NP_850840.1',
                'NP_188596.1',
                'NP_178291.1',
                'NP_197963.1',
                'NP_191404.2',
                'NP_199599.1',
                'NP_564107.1',
                'NP_567088.1',
                'NP_177334.1',
                'NP_001117417.1',
                'NP_191325.1',
                'NP_973993.1',
                'NP_001030838.1'
