-- Get some value from maximum scoring rows
-- select blastoutput_db, query_locus, max(hsp_bit_score), hsp_positive, min(mrca.phylostratum) 
-- from (
--         (blastreport inner join blastdatabase on blastreport.blastoutput_db=blastdatabase.database)
--         inner join mrca on blastreport.query_taxon=mrca.taxid_1 and blastdatabase.taxid=mrca.taxid_2
--      )
-- where query_locus in ('At1G10960') 
-- group by blastoutput_db, query_locus;


-- select blastoutput_db, query_locus, max(hsp_bit_score),
--        hsp_bit_score, hsp_positive, iteration_query_len, hit_len, hsp_align_len 
-- from blastreport
-- where query_locus in ('NP_000005') 
-- group by blastoutput_db, query_locus;


-- Get some value from maximum scoring rows
-- select blastdatabase.species, query_locus, max(hsp_bit_score)
-- from (
--         (blastreport inner join blastdatabase on blastreport.blastoutput_db=blastdatabase.database)
--         inner join mrca on blastreport.query_taxon=mrca.taxid_1 and blastdatabase.taxid=mrca.taxid_2
--      )
-- --where query_locus in ('NP_000005', 'NP_000006')
-- group by blastoutput_db, query_locus
-- order by query_locus, mrca.phylostratum, mrca.mrca
-- limit 200;


-- Get phylostrata for list of loci at given bitscore cutoff
-- DON"T USE YET, CORRECTNESS NOT TEST
-- select query_locus, min(phylostratum)
-- from
-- (
--     select blastreport.blastoutput_db,
--            blastreport.query_locus,
--            blastreport.hsp_positive,
--            mrca.phylostratum,
--            max(hsp_bit_score)
--     from (
--             (blastreport inner join blastdatabase on blastreport.blastoutput_db=blastdatabase.database)
--             inner join mrca on blastreport.query_taxon=mrca.taxid_1 and blastdatabase.taxid=mrca.taxid_2
--          )
--     where hsp_bit_score > 100 and query_locus in ('AT1G10960', 'AT4G19940', 'AT5G39860') 
--     group by blastoutput_db, query_locus
-- )
-- group by query_locus;
