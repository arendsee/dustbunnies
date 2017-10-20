#! /usr/bin/Rscript

##########################################################################
# This script adds Most Recent Common Ancestor (MRCA) columns to
# a database csv file (or any csv file that contains a "TaxonId" column).
# This is time consuming, expect more than a second per database (for desktop).
#
# INPUT:
#       query.taxid - The taxon id of the query (e.g. 3702 for A. thaliana).
#       query.species - the spoecies name of the query
#       dbfile - File name of the database info file. This file MUST contain
#                a column titled TaxonId
#       taxdir - The directory which contains the taxon dump
#
# OUTPUT:
#       Adds three columns to the original database information file:
#
#       MRCA_Id - The taxon id of each mrca
#       MRCA_Name - The scientific name of each mrca
#       MRCA_Rank - The rank of each mrca (frequently 'not ranked')

add_mrca_columns <- function(query.taxid, query.species, dbfile, taxdir) {

    if(! any(.packages(all.available=TRUE) == "CHNOSZ")){
        cat("Please install 'CHNOSZ' library\n")
        q()
    }

    db.dat <- read.csv2(file=dbfile, header=TRUE)

    if(any(db.dat$taxid == "UNKNOWN TAXID")){
        cat("Please clean the dbinfo file (at least one unknown taxid is 
             present)\n")
        q()
    }

    # Construct mrca filename
    dbfileout <- gsub("dbinfo", paste(query.species, '_mrca', sep=""), dbfile)

    # Exit if the database info file already contains MRCA headers
    if(file.exists(dbfileout)) q()

    cat("\nWriting MRCA data to dbinfo file, this will take a few minutes\n\n")

    # Load CHNOSZ library to handle taxonomy queries
    library(CHNOSZ)
    
    names <- getnames(taxdir)
    nodes <- getnodes(taxdir)

    # The lineages are ordered from species to root
    query.lineage <- allparents(query.taxid, nodes=nodes)
 
    get.mrca <- function(id, q.lineage) {

        id.lineage <- allparents(id, nodes=nodes)
    
        # The first member of the intersect will be the mrca
        return (intersect(id.lineage, q.lineage)[1])
    }

    get.stratum <- function(mrca_id, q.lineage){
        return(which(q.lineage == mrca_id))
    }
    
    # Add mrca columns to dataframe
    db.dat$MRCA_Id <- sapply(db.dat$taxid, get.mrca, query.lineage)
    db.dat$MRCA_Name <- sciname(db.dat$MRCA_Id, names=names)
    db.dat$MRCA_Rank <- getrank(db.dat$MRCA_Id, nodes=nodes)
    db.dat$Stratum <- sapply(db.dat$MRCA_Id, get.stratum, query.lineage)

    db.dat <- db.dat[order(db.dat$Stratum),]

    # Write to new file with taxon id prefix
    write.csv2(db.dat, file=dbfileout, row.names=FALSE, quote=FALSE)
}

# Read command line arguments
args <- commandArgs(TRUE)

input.good <- function(args){
    if(length(args) < 4){
        cat("Too few arguments given to find_mrca.R\n")
        return(FALSE)
    }
    if(! file.exists(args[3])){
        cat(paste0("dbfile '", args[3], "' does not exist\n"))
        return(FALSE)
    }
    if(! file.exists(args[4])){
        cat(paste0("taxonomy file '", args[4], "' does not exist\n"))
        return(FALSE)
    }
    return(TRUE)
}

if(input.good(args)){
    # Call main function
    add_mrca_columns(args[1], args[2], args[3], args[4])
}
