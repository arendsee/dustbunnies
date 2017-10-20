library(ggplot2)
library(reshape2)

# Read in database information file
# Tests for correct format
LoadInfo <- function(filename){
  d <- read.csv(file=filename, row.names=1, header=TRUE, 
                 stringsAsFactors=TRUE)
  
  stopifnot(nrow(d) > 0,
            all(c("MRCA_Name", "Stratum") %in% colnames(d)),
            all(d$Stratum %% 1 == 0) && all(d$Stratum > 0))
  
  d <- d[order(-d$Stratum, d$Taxid, rownames(d)), ] # Order dbinfo by stratum
  return(d)
}

# Read in BLAST results csv file (matrix of scores, evalues, or identity)
# Tests for correct format and correlation between database file and score
# matrix
LoadMat <- function(filename, d, low.bound=0, high.bound=Inf, na2=NA){
  m <- read.csv2(file=filename, row.names=1, header=TRUE, 
                 stringsAsFactors=FALSE)
  
  # make all entries numeric
  m <- apply(m, c(1,2), as.numeric)
  
  # Set NA to desired value (for searches in which no hit was found evalue 
  # might be reported as NA, but could almost equivalently be given a high
  # numeric value like 100)
  m[which(is.na(m))] = na2
  
  # Set data matrix columns to match dbinfo rows
  m <- m[,rownames(d)]
  
  # Stops if database names do not match or if any values are outside the
  # specified bounds
  stopifnot(identical(colnames(m), rownames(d)),
            max(m) <= high.bound,
            min(m) <= low.bound)
  
  # Stop if query headers do not contain locus information
  stopifnot(all(grepl("locus\\|\\w+", rownames(m), perl=TRUE)))
  
  return(m)
}

# Get ordered list of phylostrata
GetStrata <- function(strata){
  strata <- as.factor(strata)
  lev <- levels(strata)
  firsts <- rep(0, length(lev))
  for(i in 1:length(lev)){
    firsts[i] <- min(which(strata == lev[i]))
  }
  return(firsts)
}

# Assumes the following format for fasta headers: ...|gb|<<gb>>|...
ExtractGb <- function(header){
  header <- gsub(".*gb\\|", "", header, perl=TRUE)
  header <- gsub("\\|.*", "", header, perl=TRUE)
  
  stopifnot(length(header) > 0)
  
  return(header)
}

# Assumes the following format for fasta headers: ...|locus|<<locus>>|...
ExtractLocus <- function(header){
  header <- gsub(".*locus\\|", "", header, perl=TRUE)
  header <- gsub("\\|.*", "", header, perl=TRUE)
  
  stopifnot(length(header) > 0)
  
  return(header)
}

# The purpose of this function is to conflate alternative splicing products to 
# the highest scoring forms for each subject. This is done by creating an 
# inclusion matrix.
# BuildInclusionMat <- function(m, max.is.significant=TRUE){
#   loci <- as.character(sapply(rownames(m), ExtractLocus))
#   
#   highest <- if (max.is.significant) {which.max} else {which.min}
#   
#   column.inclusion <- function(x){
#     column <- rep(FALSE, length(x))
#     unique.loci <- unique(loci)
#     for(i in 1:length(unique.loci)){
#       ul <- which(loci == unique.loci[i]) # Find indices of locus i (1 or more)
#       mul <- highest(x[ul]) # Determine which of these forms has the best score
#       muli <- ul[mul] # Get the index of the best score
#       column[muli] = TRUE # Set this index to 1 (others will be 0) in inclusion matrix
#     }
#     return(column)
#   }
#   
#   inc.mat <- apply(m, 2, function(x) column.inclusion(x))
#   
#   # Die if all columns do not sum to the number of unique loci
#   stopifnot(all(apply(inc.mat, 2, sum) == length(unique(loci))))
#   
#   return(inc.mat)
# }
# 
# Creates a pdf with 5X5 genes ploted per page.
# PlotFew <- function(d, sm, em, im){
#   nw <- 5
#   nh <- 5
#   
#   # Draw phylostratum divisions
#   strata <- GetStrata(d$Stratum)
#   
#   # Color selection function
#   choose.colors <- function(x){
#     colors <- ifelse(x <  1e-10, "violet",   "black")
#     colors <- ifelse(x >= 1e-10, "blue",     colors)
#     colors <- ifelse(x >  1e-8,  "green",    colors)
#     colors <- ifelse(x >  1e-5,  "yellow",   colors)
#     colors <- ifelse(x >  1e-3,  "red",      colors)
#     colors <- ifelse(x >  1e-1,  "dark red", colors)
#     return(colors)
#   }
#   
#   for(i in 0:((nrow(sm) - 1) %/% (nw * nh))){
#     par(mfrow=c(nw,nh), mar=c(0,0,1,0))
#     for(j in 1:(nw*nh)){
#       index <- nw * nh * i + j
#       if(index > nrow(sm)) next
#       
#       # Convert percent identity to rounded ASCII values
#       # e.g. 0.56 --> 6 --> 54 (ASCII code for "6")
#       # 10 is converted to A (hexadecimal for 10)
#       chars <- round(im[index, ],1) * 10 + 48
#       chars <- ifelse(chars == 58, 65, chars)
#       
#       gb <- ExtractGb(rownames(sm)[index])
#       
#       plot(1:nrow(d), sm[index, ], col=choose.colors(em[index, ]), tck=0, 
#            xaxt="n", yaxt="n", pch=chars, main=paste0(gb))
#       
#       abline(v=(strata[2:(length(strata))] - 0.5), col="gray")
#     }
#   }
# }

PlotMany <- function(d, sm, em){
  rd <- d[colnames(sm), ]
  nw <- 5
  nh <- 5
  
  # Draw phylostratum divisions
  strata <- GetStrata(rd$Stratum)
  
  # Color selection function
  choose.colors <- function(x){
    colors <- ifelse(x <  1e-10, "violet",   "black")
    colors <- ifelse(x >= 1e-10, "blue",     colors)
    colors <- ifelse(x >  1e-8,  "green",    colors)
    colors <- ifelse(x >  1e-5,  "yellow",   colors)
    colors <- ifelse(x >  1e-3,  "red",      colors)
    colors <- ifelse(x >  1e-1,  "dark red", colors)
    return(colors)
  }
  
  for(i in 0:((nrow(em) - 1) %/% (nw * nh))){
    par(mfrow=c(nw,nh), mar=c(0,0,1,0))
    for(j in 1:(nw*nh)){
      index <- nw * nh * i + j
      if(index > nrow(em)) next
      
      gb <- ExtractGb(rownames(em)[index])
      
      plot(1:ncol(sm), sm[index, ], col=choose.colors(em[index, ]), tck=0, 
           xaxt="n", yaxt="n", main=paste0(gb), pch=".")
      
      abline(v=(strata - 0.5), 
             col=rgb(0, 100, 0, 50, maxColorValue=255))
    }
  }
}

# Counts the constituents of each strata at each cutoff
Phylostratify <- function(m, d, cut, high.is.significant=TRUE, 
                          contains.self=TRUE, no.first=FALSE){
  strata.names <- unique(d$MRCA_Name)
  s.num <- length(strata.names)
  c.num <- length(cut)
  strat.mat <- matrix(rep(0, s.num*c.num), nrow=c.num)
  rownames(strat.mat) <- as.character(cut)
  colnames(strat.mat) <- strata.names

  make.binary <- if(high.is.significant) 
                 { function(mat, cutoff) ifelse(mat > cutoff, 1, 0) } else
                 { function(mat, cutoff) ifelse(mat < cutoff, 1, 0) }
  
  for(i in 1:c.num){
    bin.mat <- make.binary(m, cut[i])
    if(contains.self) bin.mat[,1] <- 1
    most.distant <- apply(bin.mat, 1, function(x) max(which(x == 1)))
    strata <- d$MRCA_Name[most.distant]
    strat.mat[i, ] <- summary(strata)[strata.names]
  }
  # If do not default seqs found everywhere to the species level
  if(no.first) strat.mat <- strat.mat[, -1]
  return(strat.mat)
}

# Draws a phylostratigraph with optionally multiple cutoffs
PlotPhylostratigraph <- function(m, d, cut, title="Phylostratigraph", ...){
  rd <- d[colnames(m), ]
  strat.mat <- Phylostratify(m, rd, cut, ...)
  strata.names <- unique(rd$MRCA_Name)
  melted.strat <- melt(strat.mat)
  colnames(melted.strat) <- c("Cutoff", "Phylostratum", "Count")
  melted.strat$Phylostratum <- factor(melted.strat$Phylostratum, 
                                      levels=strata.names, ordered=TRUE)
  melted.strat$Cutoff <- as.factor(melted.strat$Cutoff)
  
  q <- ggplot(melted.strat, aes(Phylostratum, Count, group=Cutoff, 
                                colour=Cutoff)) + 
    geom_line() +
    labs(title=title) +
    theme(axis.text.x = element_text(angle=330, hjust=0))
  q
}

# Numerically summarizes results
SummarizePhylostratigraph <- function(m, d, cut, ...){
  strat.mat <- Phylostratify(m, d, cut, ...)
  strata.names <- unique(d$MRCA_Name)
  
}

# Takes a score matrix and reduces it in accordance with the inclusion matrix I
# Trim <- function(m, I) {
#   trimmed <- matrix(rep(0, ncol(m)*sum(I[,1])), nrow=sum(I[,1]))
#   rownames(trimmed) <- unique(as.character(sapply(rownames(m), ExtractLocus)))
#   colnames(trimmed) <- colnames(m)
#   
#   for(i in 1:ncol(m)){ trimmed[,i] <- m[which(I[,i]),i] }
#   
#   return(trimmed)
# }
# 


##################################################

# Read And Order a score matrix
rao <- function(f, d) {
    m <- read.csv(file=f, row.names=1)
    m <- m[, intersect(rownames(d), colnames(m))] 
    m <- m[sort(rownames(m)),]
    return(m)
}


##################################################
# Particular

setwd("~/ohome/workspace/hs/")

d <- LoadInfo("dbinfo.csv")
sm <- rao("score.csv", d)
em <- rao("evalue.csv", d)

# Extreme cutoffs
xecut <- (10 ^ (2 * c(-0.5, -1:-20)))
xscut <- (10 * (1:30))

# moderate cutoffs
mecut <- (10 ^ (-12:-3))
mscut <- (5 * (8:20))


pdf("all.pdf")
PlotPhylostratigraph(em, d, xecut, "Phylostratigrah (e-value cutoff)",
                     high.is.significant=FALSE, contains.self=TRUE)
PlotPhylostratigraph(sm, d, xscut, "Phylostratigrah (bit-score cutoff)",
                     high.is.significant=TRUE, contains.self=TRUE)
PlotPhylostratigraph(em, d, mecut, "Phylostratigrah (e-value cutoff)",
                     high.is.significant=FALSE, contains.self=TRUE)
PlotPhylostratigraph(sm, d, mscut, "Phylostratigrah (bit-score cutoff)",
                     high.is.significant=TRUE, contains.self=TRUE)
dev.off()


pdf('all_loci.pdf')
PlotMany(d, sm, em)
dev.off()



###### Scoring Schemes ######
setwd('~/ohome/workspace/score_schemes/')
d <- read.csv('at_besthits.csv', stringsAsFactors=FALSE)
dbinfo <- LoadInfo('at.dbinfo.csv')


# TODO specify the id.vars
d.melt <- melt(d)
adat <- acast(d.melt, qgb ~ database ~ variable)
adat <- adat[ , rownames(dbinfo), ]

mmat <- adat[,,'mscore']
smat <- adat[,,'sscore']
pmat <- adat[,,'pscore']

mscut <- c(30, 35, 40, 50, 75, 100, 200, 500, 1000)
PlotPhylostratigraph(pmat, dbinfo, mscut, "Best Path (bit-score cutoff)",
                     high.is.significant=TRUE, contains.self=TRUE)
ggsave("at_bestpath_extreme.jpeg")
# PlotPhylostratigraph(smat, dbinfo, mscut, "Summed Scores (bit-score cutoff)",
#                      high.is.significant=TRUE, contains.self=TRUE)
# PlotPhylostratigraph(mmat, dbinfo, mscut, "Max Hsp (bit-score cutoff)",
#                      high.is.significant=TRUE, contains.self=TRUE)

CompareScoreMethods <- function(a, dbinfo, cutoff=c(50, 75, 100),
                                ggtitle='Method comparison',
                                filename='methods.tiff'){
    mmat <- a[,,'mscore']
    smat <- a[,,'sscore']
    pmat <- a[,,'pscore']

    stratmelt <- function(m, d, method){
        strat <- Phylostratify(m, d, cutoff)
        strat.melt <- melt(strat)
        colnames(strat.melt) <- c('Cutoff', 'Phylostratum', 'Count')
        strat.melt$Method <- rep(method, length(strat.melt))
        return(strat.melt)
    }

    mstrat <- stratmelt(mmat, dbinfo, 'Maximum') 
    pstrat <- stratmelt(pmat, dbinfo, 'Bestpath')
    sstrat <- stratmelt(smat, dbinfo, 'Summed')

    stratdat.melt <- rbind(mstrat, pstrat, sstrat)

    stratdat.melt$Phylostratum <- factor(stratdat.melt$Phylostratum, levels=unique(dbinfo$MRCA_Name), ordered=TRUE)

    number_ticks <- function(n) {function(limits) pretty(limits, n)}

    ggplot(stratdat.melt, aes(x=Phylostratum, y=Count, group=Method, colour=Method)) + 
        geom_point(alpha=0.7) +
        geom_line(alpha=0.2) +
        theme(axis.text.x = element_text(angle=330, hjust=0)) +
        scale_y_continuous(breaks=number_ticks(10), limits=c(0,1000)) +
        scale_x_discrete(limits=unique(dbinfo$MRCA_Name)[1:15]) +
        ggtitle(ggtitle) + 
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        facet_grid(. ~ Cutoff)
    ggsave(filename)
}

CompareScoreMethods(adat, dbinfo, ggtitle='Scoring method comparison (Hs)', filename='hs_methods_lim.jpeg')

###### Controls #######

setwd("~/ohome/workspace/control")

d <- LoadInfo('control_dbinfo.csv')
srtk = rao('1000.reverse.tantan.score.csv', d)
srsk = rao('1000.reverse.segmasker.score.csv', d)
srk  = rao('1000.reverse.score.csv', d)
sk   = rao('1000.score.csv', d)
stk  = rao('1000.tantan.score.csv', d)
ssk  = rao('1000.segmasker.score.csv', d)

jpeg('masking_differences.jpeg', width=960, height=960)
par(mfrow=c(2,1))
boxplot(sk - stk, names=c(1:ncol(sk)), main="Tantan")
boxplot(sk - ssk, names=c(1:ncol(sk)), main="Seq")
dev.off()

threshold <- 100

tstrat <- Phylostratify(stk, d[colnames(stk), ], threshold, TRUE, TRUE, FALSE)
sstrat <- Phylostratify(ssk, d[colnames(ssk), ], threshold, TRUE, TRUE, FALSE)
nstrat <- Phylostratify(sk, d[colnames(sk), ], threshold, TRUE, TRUE, FALSE)

stratdat <- data.frame(ps=unique(d$MRCA_Name), tantan=as.numeric(tstrat), seg=as.numeric(sstrat), unmasked=as.numeric(nstrat))

stratdat.melt <- melt(stratdat)
colnames(stratdat.melt) <- c('Phylostratum', 'Masker', 'Count')

stratdat.melt$Phylostratum <- factor(stratdat.melt$Phylostratum, levels=unique(d$MRCA_Name), ordered=TRUE)

ggplot(stratdat.melt, aes(x=Phylostratum, y=Count, group=Masker, colour=Masker)) + 
    geom_line() +
    theme(axis.text.x = element_text(angle=330, hjust=0))
ggsave('masker_ps.jpeg')


###### Plot control phylostratigraphs ######

fcut <- (5 * (8:20))
rcut <- (5 * (8:12))

pdf('control_strats.pdf')
PlotPhylostratigraph(srtk, d, rcut, "rtk", TRUE, TRUE, TRUE)
PlotPhylostratigraph(srsk, d, rcut, "rsk", TRUE, TRUE, TRUE)
PlotPhylostratigraph(srk, d, rcut, "rk", TRUE, TRUE, TRUE)
PlotPhylostratigraph(stk, d, fcut, "tk", TRUE, TRUE, FALSE)
PlotPhylostratigraph(ssk, d, fcut, "sk", TRUE, TRUE, FALSE)
PlotPhylostratigraph(sk, d, fcut, "k", TRUE, TRUE, FALSE)
dev.off()



###### Plot all control loci ######

ertk = rao('1000.reverse.tantan.evalue.csv', d)
ersk = rao('1000.reverse.segmasker.evalue.csv', d)
erk  = rao('1000.reverse.evalue.csv', d)
ek   = rao('1000.evalue.csv', d)
etk  = rao('1000.tantan.evalue.csv', d)
esk  = rao('1000.segmasker.evalue.csv', d)

pdf('k.pdf')
PlotMany(d, sk, ek)
dev.off()

pdf('tk.pdf')
PlotMany(d, stk, etk)
dev.off()

pdf('sk.pdf')
PlotMany(d, ssk, esk)
dev.off()



####### Alignment region visualization ######

hsp.dat <- read.csv('1000_hsp_data.csv')
s <- hsp.dat[which(hsp.dat$query == 'AT1G04680'), ]
# Sort rows by phylostratum
rn <- intersect(rownames(d), s$database)
rd <- d[rn, ] 
s <- s[as.numeric(sapply(rn, function(x) which(s$database == x))), ]
strata <- GetStrata(rd$Stratum)

q.len <- s$iteration_query_len
q.to <- s$hsp_query_to
q.from <- s$hsp_query_from
q.alen <- abs(q.to - q.from) + 1

b.low <- pmin(q.to, q.from) / q.len
b.high <- pmax(q.to, q.from) / q.len

ave.score <- s$hsp_bit_score / q.alen
ave.score[is.na(ave.score)] = 0
norm.ave.score <- ave.score / max(ave.score)


require(grDevices)
require(plotrix)
i <- (length(b.low) - 1):0
plot(c(0, 1), c(0, 63), type = "n", xlab = "", ylab = "", ann=FALSE, xaxt="n", yaxt="n")
ncol <- 10
col <- rainbow(ncol)
bg.col <- c('black', 'gray')
rect(0, i, 1, 1 + i, col=bg.col[i %% 2 + 1], border=NA)
rect(b.low, i, b.high, 1 + i, col=col[ceiling(ncol * norm.ave.score)])
color.legend(0.9, -3, 1, 0, legend=c("Low", "High"), rect.col=col, align='rb')
abline(h=(length(b.low) - c(0, strata, length(b.low))), col='orange')
