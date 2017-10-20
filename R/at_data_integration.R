# Merge data from many Arabidopsis thaliana datasets. None of this has been
# tested (and some I know won't work) on data from other species.

home <- '~/ohome/lib/R/'
source(paste0(home, 'visual_functions.R'), chdir=TRUE)
source(paste0(home, 'io.R'), chdir=TRUE)
source(paste0(home, 'analyses.R'), chdir=TRUE)
source(paste0(home, 'phylogenetic_functions.R'), chdir=TRUE)

# Input filenames
UNMASKED.FILE <- '~/ohome/projects/orphans/blast_csvs/unmasked_at_g50.csv'
MASKED.FILE   <- '~/ohome/projects/orphans/blast_csvs/masked_at_g50.csv'
REGULON.FILE  <- '~/ohome/projects/orphans/regulon/regulon.csv'
BLINK.FILE    <- '~/ohome/projects/blink/At_cut100-org2-set0_2013-09-11.csv'
CRE.FILE      <- '~/ohome/data/genomes/at/supp_data/confidence_rankings/cr_exon_TAIR10.csv'
CRG.FILE      <- '~/ohome/data/genomes/at/supp_data/confidence_rankings/cr_gene_TAIR10.csv'
PROTEIN.FILE  <- '~/ohome/data/genomes/at/supp_data/protein-data_2013-10-25.csv'
GO.FILE       <- '~/ohome/data/genomes/at/supp_data/GO/GO.csv'
PO.FILE       <- '~/ohome/data/genomes/at/supp_data/PO/PO.csv'
DESC.FILE     <- '~/ohome/data/genomes/at/supp_data/desc/desc.csv' 
OUTPUT.FILE   <- '~/ohome/data/genomes/at/supp_data/at_alldata_2013-xx-xx.csv'

# Required dataframe columns
BLINK.COL <- c('locus', 'model', 'length', 'nhits', 'nproteins', 'nspecies', 'gi', 'gb',
              'Archae', 'Bacteria', 'Metazoa', 'Fungi', 'Plants', 'Viruses', 'Other_eukaryotes')
CRE.COL <- c('Model', 'Start', 'Stop', 'Overall_confidence', 'AGI', 'Proteomics', 'X.species', 'VISTA')
CRG.COL <- c('Model', 'Exon_class', 'Confidence', 'Overall', 'AGI', 'Proteomics', 'X.species', 'VISTA')
PROTEIN.COL <- c('Model', 'MW', 'pI', 'Location', 'TM.Domains', 'Structural.Class')
REGULON.COL <- c('AffyID', 'LocusID', 'Regulon', 'Function')
GO.COL <- c('Locus', 'GO.Term', 'GO.ID', 'cat', 'GO.Slim', 'Code', 'Reference', 'Date')
PO.COL <- c('Locus', 'PO.Term', 'PO.ID', 'cat', 'Code', 'Reference', 'Date')
DESC.COL <- c('model', 'desc')



LoadFullBesthits <- function(file.=UNMASKED.FILE){
    out. <- read.csv(file.)
    return(out.)
}

LoadPrepData <- function(masked.=MASKED.FILE,
                         unmasked.=UNMASKED.FILE,
                         melt.=FALSE){
    out. <- PrepareDataframe(masked=masked., unmasked=unmasked.)
    if(melt.){
        out. <- MeltPrepData(out.)
    }
    return(out.)
}

LoadRegulonData <- function(file.=REGULON.FILE){
    rd. <- read.csv(file.)
    stopifnot(setequal(colnames(rd.), REGULON.COL))
    rd.$LocusID <- toupper(rd.$LocusID)
    return(rd.)
}

LoadBlinkData <- function(file.=BLINK.FILE){
    b. <- read.csv(file.) 
    stopifnot(setequal(colnames(b.), BLINK.COL))
    b.$Orphan <- ifelse(b.$nspecies <= 1, TRUE, FALSE)
    return(b.)
}

LoadProteinData <- function(file.=PROTEIN.FILE){
    p. <- read.csv(file., stringsAsFactors=FALSE)
    stopifnot(setequal(colnames(p.), PROTEIN.COL))
    s. <- as.character(p.$Structural.Class)
    p.$Structural.Class <- ifelse(s. == "", 'undefined', s.)
    return(p.)
}

# By locus
# ';' delimited csv
LoadGO <- function(file.=GO.FILE){
    go. <- read.csv2(file., stringsAsFactors=TRUE) 
    stopifnot(setequal(colnames(go.), GO.COL))
    return(go.)
}

# By locus
LoadPO <- function(file.=PO.FILE){
    po. <- read.csv(file., stringsAsFactors=TRUE) 
    stopifnot(setequal(colnames(po.), PO.COL))
    return(po.)
}

LoadGeneConfidence <- function(file.=CRG.FILE){
    crg. <- read.csv(file.)
    stopifnot(setequal(colnames(crg.), CRG.COL))
    return(crg.)
}

# Used in calculating exon numbers
LoadExonConfidence <- function(file.=CRE.FILE){
    cre. <- read.csv(file.)
    stopifnot(setequal(colnames(cre.), CRE.COL))
    return(cre.)
}

LoadDescData <- function(file.=DESC.FILE){
    desc <- read.csv(file.)
    stopifnot(setequal(colnames(desc), DESC.COL))
    return(desc)
}

MeltPrepData <- function(d.){
    m. <- melt(d., id.vars=c('qgb', 'hspecies', 'mrca', 'ps'))
    return(m.)
}

# This works in Arabidopsis' specific naming scheme ONLY, where the third
# character in the ID is the chromosome number or organelle letter.
GetBlinkOrganelleOrphans <- function(d.){
    orphan <- d.[which(d.$Orphan), ]
    orphan$chr <- substr(orphan$locus, 3, 3) 
    org.orphan <- orphan[which(orphan$chr %in% c('M', 'C')), 'gb']
    return(as.character(org.orphan))
}

# Union of my orphans and Blink organelle orphans
GetOrphanset <- function(m, cutoff=1e-5, fill='umeval'){
    # Unmasked e-val=1e-5 orphanset
    uo <- GetLSG.by.evalue(m, cutoff=cutoff, fill=fill)
    uo <- uo[which(uo$mrca == 'Arabidopsis thaliana'), 'names']
    # Add blink organelle orphans to my orphanset
    blink <- LoadBlinkData()
    boo <- GetBlinkOrganelleOrphans(blink)
    buo <- union(uo, boo) # Blink Unmasked Orphanset
    return(buo)
}

GetOntologyCount <- function(o.){
    count <- data.frame(locus=levels(o.$Locus),
                        count=as.numeric(summary(o.$Locus, maxsum=Inf)))
    return(count)
}

GetExonCounts <- function(cre.){
    ec <- data.frame(model=levels(cre.$Model),
                     count=as.numeric(summary(cre.$Model, maxsum=Inf)))
    return(ec)
}

# winnow - remove any loci that have null orphan values (in my specific dataset
# there are 434 loci in the regulon data that do not appear in TAIR10, most
# seem to be transposable elements and pseudogenes.
LoadAll <- function(winnow=TRUE, ...){
    # Base dataframe
    d <- LoadBlinkData()[, c('locus', 'model', 'gi', 'gb', 'length')]

    # Add orphan column
    dm       <- LoadPrepData(melt=TRUE)
    oset     <- GetOrphanset(dm, ...)
    d$orphan <- ifelse(d$gb %in% oset, TRUE, FALSE)

    # Add exon counts
    ce <- GetExonCounts(LoadExonConfidence())
    d <- merge(d, ce, by.x='model', by.y='model', all=TRUE)

    # Add protein data
    protein  <- LoadProteinData()
    d <- merge(d, protein, by.x='model', by.y='Model', all=TRUE)

    # BY LOCUS
    # Add PO and GO counts
    go <- GetOntologyCount(LoadGO())
    po <- GetOntologyCount(LoadPO())
    colnames(go) <- c('locus', 'GO.terms')
    colnames(po) <- c('locus', 'PO.terms')
    d <- merge(d, go, all=TRUE)
    d <- merge(d, po, all=TRUE)

    # BY LOCUS
    # Add regulon data
    regulon  <- LoadRegulonData()[, c('LocusID', 'Regulon', 'Function')]
    d <- merge(d, regulon, by.x='locus', by.y='LocusID', all=TRUE)

    # Add gene desctiption strings
    desc <- LoadDescData()
    d <- merge(d, desc, all=TRUE)

    if(winnow){
        d <- d[-which(is.na(d$orphan)), ]
    }

    # Final adjustment of names
    # KLUDGY way of doing things...
    colnames(d) <- c("model","locus","gi","gb","length","orphan","exon.count","MW","pI","location",
                     "transmembrane.domains","structural.class","GO.terms","PO.terms","regulon.number",
                     "regulon.function","desc")

    return(d)
}

WriteAll <- function(d., file=OUTPUT.FILE){
    write.csv(d., file=OUTPUT.FILE)
}

d <- LoadAll()
WriteAll(d, OUTPUT.FILE)
