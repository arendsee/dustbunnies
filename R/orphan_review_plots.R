home <- '~/ohome/lib/R/'
source(paste0(home, 'visual_functions.R'), chdir=TRUE)
source(paste0(home, 'io.R'), chdir=TRUE)
source(paste0(home, 'analyses.R'), chdir=TRUE)
source(paste0(home, 'phylogenetic_functions.R'), chdir=TRUE)

###### A basic phylostratigraphy plot
Phylostratigraph_review.2013.10.23 <- function(m, cutoffs=c(1e-5), fill='umeval'){
    stratdat <- ToStratCountEvalue(m, cutoffs=cutoffs, fill=fill)
    g <- ggplot(stratdat, aes(x=stratum, y=count, group=cutoff)) +
        theme_bw() +
        geom_line() +
        theme(
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm"), 
            legend.position='none',
            axis.text.x = element_text(angle=330, hjust=0),
            axis.title = element_blank()
             )
    return(g)
}



###### A conservative/liberal ribbon plot
Ribbon_review.2013.10.23 <- function(m, lowcut=1e-4, highcut=1e-7){
    lowstrat <- ToStratCountEvalue(m, cutoffs=lowcut, fill='umeval')
    highstrat <- ToStratCountEvalue(m, cutoffs=highcut, fill='mmeval')
    ribbondat <- data.frame(ps=1:nrow(lowstrat),
                            stratum=lowstrat$stratum,
                            min=pmin(lowstrat$count, highstrat$count),
                            max=pmax(lowstrat$count, highstrat$count))
    g <- ggplot() +
        geom_blank(data=ribbondat, aes(x=stratum)) +
        geom_ribbon(data=ribbondat, aes(ymin=min, ymax=max, x=ps), alpha=0.30) +
        theme_bw() +
        theme(
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm"), 
            legend.position='none',
            axis.text.x = element_text(angle=330, hjust=0),
            axis.title = element_blank()
             )

    return(g)
}



###### Orphan count: changes with parameters
# Assumes none, and segmasked
# e-value based
# OrphanCounts_review.2013.10.23 <- function(m, cutoffs=10^(-(1:10))){
#     ud <- ToStratCountEvalue(m, cutoffs=cutoffs, fill='umeval')
#     ud <- ud[which(ud$stratum == 'Arabidopsis thaliana'), ]
#     md <- ToStratCountEvalue(m, cutoffs=cutoffs, fill='mmeval')
#     md <- md[which(md$stratum == 'Arabidopsis thaliana'), ]
#     strat <- cbind(ud$count, md$count)
#     colnames(strat) <- c('none', 'seg')
#     rownames(strat) <- md$cutoff
#     return(strat)
# }



# obo plot of orphans (conventional, by score)
Plotobo_review.2013.10.23  <- function(orphans, title_='unmasked_at_orphans_e-5.pdf'){
    plotobo(orphans, title_=title_)
}



###### Regulon data
# LoadRegulonData_review.2013.10.23 <- function(filename='~/ohome/projects/orphans/regulon/regulon.csv'){
#     rd <- read.csv(filename)
#     rd$LocusID <- toupper(rd$LocusID)
#     return(rd)
# }



# A weird plot, probably won't use
PlotOrphanCounts_review.2013.10.23 <- function(m, ...){
    counts <- as.data.frame(OrphanCounts_review.2013.10.23(m, ...))
    counts.melt <- data.frame(evalue=as.factor(as.numeric(rep(rownames(counts), 2))),
                         mask=c(rep('none', nrow(counts)), rep('seg', nrow(counts))),
                         count=c(counts$none, counts$seg))
    g <- ggplot(counts.melt, aes(x=evalue, y=count, group=mask)) +
        geom_point(aes(color=mask)) +
        geom_path(alpha=I(0.2))
    g
    return(g) 
}



# WriteLingCsv_review.2013.10.23 <- function(gbs, full, file){
#     out <- full[which(full$qgb %in% gbs), c('qgb', 'qlocus', 'qlen')]
#     out <- unique(out)
#     colnames(out) <- c('id', 'locus', 'length')
#     write.table(out, file=file, quote=FALSE, sep=',', row.names=FALSE)
# }


# Binomial test for significance of deviation from regulon expectation
# One test for each of 72 bins
# Input: d    - all data set
#        m    - number of permutation trials to run
#        last - last regulon to consider as independent bin, all later ones will
#               be gouped together
RegulonAnlysis_review.2013.10.23 <- function(d., loci, last=71){
    # Setup regulon dataset
    reg <- d.[-which(is.na(d.$regulon.number)), c('locus', 'regulon.number', 'regulon.function')]
    colnames(reg) <- c('LocusID', 'Regulon', 'Function')
    lowset <- paste0(last + 1, ':', max(d.$regulon.number, na.rm=TRUE))
    reg$Regulon <- ifelse(reg$Regulon <= last, reg$Regulon, lowset)
    # Number of elements in each regulon
    reg.n <- as.numeric(summary(factor(reg$Regulon), maxsum=1e6))
    reg.l <- levels(factor(reg$Regulon))
    k <- length(reg.l)

    reg.loci <- intersect(loci, reg$LocusID)
    n <- length(reg.loci)
    probs <- reg.n / sum(reg.n)

    # Function that counts the elements in the input factor
    pile <- function(a){
        as.numeric(summary(factor(a, levels=reg.l), maxsum=1e6))
    }
    observed <- pile(reg[which(reg$LocusID %in% loci), 'Regulon'])
    bin <- data.frame(x=observed,
                      n=rep(n, length(reg.n)),
                      prob=probs)    
    pvals <- apply(bin, 1, function(x) binom.test(x[1], x[2], x[3], 
                                                  alternative='two.sided')$p.val)
    pvals <- data.frame(Regulon=as.character(reg.l),
                        Size=reg.n,
                        Observed=observed,
                        Expected=round(n * probs),
                        P.val=pvals,
                        Sign=ifelse(observed > n * probs, '+', '-'))
    pvals <- merge(pvals, unique(reg[, 2:3]))
    pvals <- pvals[order(pvals$P.val), ]
    return(pvals)
}

OrphanRegulonAnlysis <- function(d., ...){
    loci <- as.character(d.[which(d.$orphan), 'locus'])
    pval <- RegulonAnlysis_review.2013.10.23(d., loci, ...)
    return(pval)
}



CompOne_review.2013.10.23 <- function(d., field., title.='Unnamed', log.=FALSE){
    dat <- d.[, c('gb', 'orphan', field.)]
    dat$orphan <- ifelse(dat$orphan, 'Orphan', 'Not Orphan')
    colnames(dat) <- c('gb', 'variable', 'value')
    g <- ggplot(dat, aes(x=variable, y=value)) +
        geom_violin(fill='gray') +
        theme_bw() +
        ggtitle(title.) +
        theme(
            axis.title=element_blank()
             )
    if(log.){
        g <- g + scale_y_continuous(trans='log2',
                               breaks=trans_breaks('log2', function(x) round(2^x)))
    }
    return(g)
}



BarplotOne_review.2013.10.23 <- function(d., field., title.='Unnamed'){
    dat <- d.[, c('gb', 'orphan', field.)]
    dat$orphan <- ifelse(dat$orphan, 'Orphan', 'Not Orphan')
    colnames(dat) <- c('gb', 'variable', 'value')
    g <- ggplot(dat, aes(x=value)) +
        geom_bar() +
        ggtitle(title.) +
        theme(
            axis.title=element_blank(),
            axis.text.x = element_text(angle=330, hjust=0)
             ) +
        facet_grid(variable ~ ., scale='free')
    return(g)
}


###### Plots for review

masked <- '~/ohome/projects/orphans/blast_csvs/masked_at_g50.csv'
unmasked <- '~/ohome/projects/orphans/blast_csvs/unmasked_at_g50.csv'
d <- PrepareDataframe(masked=masked, unmasked=unmasked)
m <- melt(d, id.vars=c('qgb', 'hspecies', 'mrca', 'ps'))

# Ribbon plot, shows uncertainty in phylostratigraph
g.ribbon <- Ribbon_review.2013.10.23(m, lowcut=1e-4, highcut=1e-7)

# Basic phylostratigraph
g.phylo <- Phylostratigraph_review.2013.10.23(m, cutoffs=c(1e-5), fill='umeval')

ggsave('ribbon_phylostratigraph.pdf', g.ribbon)
ggsave('simple_phylostratigraph_e5_unmasked.pdf', g.phylo)

big <- read.csv('~/ohome/data/genomes/at/supp_data/at_alldata_2013-20-27.csv')

# Violin: for length, mw, and pI
# TODO codon optimization, aa composition, complexity
g.length <- CompOne_review.2013.10.23(big, field.='length',     log.=TRUE, title.='Length')
g.MW     <- CompOne_review.2013.10.23(big, field.='MW',         log.=TRUE, title.='Molecular Weight')
g.pI     <- CompOne_review.2013.10.23(big, field.='pI',                    title.='pI')
g.go     <- CompOne_review.2013.10.23(big, field.='GO.terms',   log.=TRUE, title.='Number of GO-terms')
g.po     <- CompOne_review.2013.10.23(big, field.='PO.terms',   log.=TRUE, title.='Number of PO-terms')
g.exon   <- CompOne_review.2013.10.23(big, field.='exon.count', log.=TRUE, title.='Number of exons')

# Barplot: localization, structure
big$location <- ifelse(big$location == 'undetermined', 'undefined', as.character(big$location))
g.loc <- BarplotOne_review.2013.10.23(big, 'location', title.='Location')
g.struct <- BarplotOne_review.2013.10.23(big, 'structural.class', title.='Protein Structure')

pdf('plots_2013-10-28.pdf')
g.length
g.MW
g.pI
g.go
g.po
g.exon
g.loc
g.struct
dev.off()

# Table: regulon
pvals <- OrphanRegulonAnlysis(big) 
write.csv(pvals, 'uo_reg.csv')


###### Diagnostic plots

# Usefull for seeing difference between masked and unmasked searches 
# mu.dif <- d[which(d$qgb %in% setdiff(mo$qgb, uo$qgb)), ]
# mu.dif <- melt(mu.dif, id.vars=c('qgb', 'hspecies', 'mrca', 'ps'))
# Plotobo_review.2013.10.23(mu.dif, title_='mask_only_orphans_1e-5.pdf')


# BLink data
# Add is_my_orphan line to Ling's stupid BLink file
# b <- read.csv(file='At_cut100-org2-set0_zo_2013-10-24.csv')
# bo.gb <- unique(b[which(b$nspecies <= 1), 'gb'])
# b$zeb_orphan_e5_unmasked <- b$gb %in% uo.gb
# b$blink_orphan <- b$nspecies <= 1
# write.csv(b, 'At_cut100-org2-set0_zo_2013-10-24.csv', quote=FALSE)
# 
# bor <- rd[which(rd$LocusID %in% b[which(b$gb %in% bo.gb), 'locus']), ]
# 
# bo.pval <- RegulonPermutationAnlysis_review.2013.10.23(bor$LocusID, rd, m=1e5)
# write.csv(bo.pval, 'bo_e5_pval.csv', quote=TRUE)
# 
# b.only <- bo.gb[-which(bo.gb %in% uo.gb)]
# u.only <- uo.gb[-which(uo.gb %in% bo.gb)]
# WriteLingCsv_review.2013.10.23(b.only, full, 'blink_only_2013-10-24.csv')
# WriteLingCsv_review.2013.10.23(u.only, full, 'mine_only_2013-10-24.csv')
# WriteLingCsv_review.2013.10.23(uo.gb, full, 'my_orphans_1e-5_unmasked_2013-10-24.csv')
# WriteLingCsv_review.2013.10.23(bo.gb, full, 'blink_orphans_2013-10-24.csv')


###### Time based phylostratigraph
# INPUT: 
# qgb|qlocus|hspecies|phylostratum|mrca|mrca_sciname|qlen|palen|malen|salen|pscore|mscore|sscore
# Of these, qgb, mrca_sciname, and evalue are needed
blastcsv.filename = '~/ohome/data/blast_csvs/unmasked_at_50g_2013-11-15.csv'
d <- OrderedFullDataframe(blastcsv.filename)
colnames(d)[4] = 'ps'
dm <- melt(d, id.vars=c('qlocus', 'qgb', 'hspecies', 'ps', 'mrca', 'mrca_sciname'))
strat <- GetLSG.by.evalue(dm, cutoff=1e-5, sciname=TRUE)
counts <- as.numeric(summary(strat$mrca)[unique(d$mrca_sciname)])
plot(counts, type="l")

# Divergence times, median times from timetree.org
# 19, 0
# 18, 5.3
# 17, NA
# 16, 16.4
# 15, NA
# 14, 95, 108.5, 90, (Cc, Eg, Gr) I'll settle with 100
# 13, 108, 114 (Cs, Vv) Settle with 110
# 12, 125
# 10, 147.8
# 9, 329
# 7, 406, 593.4 (median, expert)
# 6, 593.4 (expert)
# 3, 936
# 2, 1628
# 1, 2520
# 0, 3480 (Nora Noffke, 2013)
 
# lump( 17, 16), (15,14), (7,6)
d2 <- d
join.strata <- function(d., l1, l2){
    sn <- unique(d.[which(d$ps %in% c(l1, l2)),c('ps','mrca','mrca_sciname')])
    sn$mrca_sciname <- as.vector(sn$mrca_sciname)
    d.$ps <- ifelse(d.$ps == sn[1,1], sn[2,1], d.$ps)
    d.$mrca <- ifelse(d.$mrca == sn[1,2], sn[2,2], d.$mrca)
    d.$mrca_sciname <- ifelse(as.vector(d.$mrca_sciname) == sn[1,3], sn[2,3], as.vector(d.$mrca_sciname))
    return(d.)
}
d2 <- join.strata(d2, 17, 16)
d2 <- join.strata(d2, 15, 14)
d2 <- join.strata(d2, 7, 6)
dm2 <- melt(d2, id.vars=c('qlocus', 'qgb', 'hspecies', 'ps', 'mrca', 'mrca_sciname'))
strat2 <- GetLSG.by.evalue(dm2, cutoff=1e-5, sciname=TRUE)
counts2 <- as.numeric(summary(strat2$mrca)[unique(d2$mrca_sciname)])
plot(counts2, type="l")

d.time <- rbind(
    c(19, 0),
    c(18, 5.3),
    c(16, 16.4),
    c(14, 100),
    c(13, 110),
    c(12, 125),
    c(10, 147.8),
    c(9, 329),
    c(6, 593.4),
    c(3, 936),
    c(2, 1628),
    c(1, 2520),
    c(0, 3480))
d.time <- data.frame(
            ps=d.time[1:(nrow(d.time)-1), 1],
            time=d.time[1:(nrow(d.time)-1), 2],
            dif=d.time[2:nrow(d.time), 2] - d.time[1:(nrow(d.time)-1), 2])

d2 <- merge(d2, d.time, sort=FALSE)

logtime <- log2(sort(unique(d2$time)))
logtime <- ifelse(logtime < 0, 0, logtime)
plot(x=logtime, y=counts2, type="l")
plot(x=d.time$time, y=counts2, type="l")
perMA <- counts2/d.time$dif
plot(x=logtime, y=perMA, type="l")
plot(x=d.time$time, y=perMA, type="l")
plot(x=1:length(perMA), y=perMA, type="l")




