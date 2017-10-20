library(ggplot2)
library(reshape2)

filename <- '~/ohome/workspace/at/chosen/all.csv'

d <- read.csv(filename)
d$align_len <- abs(d$hsp_query_from - d$hsp_query_to)
d$qcov <-  d$align_len / d$iteration_query_len

d.melt <- melt(d, id.vars=c('query', 'database'))

# A 3d array with x = query locus, y = database species name, z = field
# (e.g. hsp_bit_score)
d.array <- acast(d.melt, query ~ database ~ variable, value.var='value')

Mat2Bin <- function(da, column, cutoff){ #, f=function(x){x > cutoff}){
    mat <- da[ , , column]
    bin <- t(apply(mat, 1, function(x) as.numeric(x) > cutoff))
    return(bin)
}

Stratify <- function(bin, ps){
    apply(bin, 1, function(x) ifelse(sum(x) == 0, '-1', min(as.numeric(ps$mrca.phylostratum[which(x)]))))
}

ps <- as.data.frame(d.array[1, , c('mrca.phylostratum', 'mrca.mrca', 'taxid2name.sciname')],
                    stringsAsFactors=FALSE)

sciname <- c(unique(ps)$taxid2name.sciname, 'Unclassifiable')
rows <- c(unique(ps)$mrca.phylostratum, '-1')
ps2name <- data.frame(sciname, row.names=rows)

smat <- Mat2Bin(d.array, 'hsp_bit_score', 100)
cmat <- Mat2Bin(d.array, 'qcov', 0.5)

bin <- smat & cmat

strata <- Stratify(bin, ps)
strata.sciname <- as.factor(ps2name[strata,])
summary(strata.sciname)

# Assert that the number of classified objects is identical to the total number
# of queryies
stopifnot(length(unique(d$query)) == sum(as.numeric(summary(as.factor(ps2name[strata,])))))


# orphans 30 for Hs, 19 for At
orphans <- strata[which(as.numeric(strata) == 19)]
orphan.loci <- attributes(orphans)$names

unclassified <- strata[which(as.numeric(strata) == -1)]
unclassified.loci <- attributes(unclassified)$names
