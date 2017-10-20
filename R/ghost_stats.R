#! /usr/bin/Rscript

require(reshape2)
require(ggplot2)

at.stat.file <- 'at_100-bit-ghosts_stat.csv'
hs.stat.file <- 'hs_100-bit-ghosts_stat.csv'

aas <- c('L','S','E','V','G','K','A','D','R','I',
         'T','P','N','F','Q','Y','M','H','C','W')
ids <- c('gi', 'locus', 'species')

myread <- function(f, species){
    d <- read.csv(f, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
    for(aa in aas){
        eval(parse(text=paste0('d$', aa, ' <- d$', aa, ' / d$length')))
    }
    d$species <- rep(species, nrow(d))
    return(d)
}

PropMasked <- function(d){ return(d$X / d$length) }

at.stat <- myread(at.stat.file, 'Arabidopsis')
hs.stat <- myread(hs.stat.file, 'Human')

d <- rbind(at.stat, hs.stat)
d <- melt(d, id.vars=ids)

d.aa <- d[-which(d$variable == 'length'), ]
d.aa$variable <- factor(d.aa$variable, levels=aas)
d.len <- d[which(d$variable == 'length'), ]

g <- ggplot(d.aa, aes(x=variable, y=value)) +
    geom_boxplot() +
    scale_y_log10() + 
    facet_grid(species ~ .)
g

g <- ggplot(d.len, aes(x=species, y=value)) +
    geom_boxplot(notch=TRUE) +
    ylim(0,200)
g

g <- ggplot(d.len, aes(x=species, y=value)) +
    geom_violin() +
    scale_y_continuous(trans="log2") +
    xlab("Species") +
    ylab("Protein length") +
    ggtitle("100-bit-ghosts") +
    theme(legend.position="none") +
    geom_point(aes(alpha=0.05), position='jitter')
g


