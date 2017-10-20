#! /usr/bin/Rscript

#####################################################################
#
# Takes tab delimited output from a blast run (-outfmt 7) and extracts
# query-database best hit matrices. I.e. for each query-database pair, 
# finds the best hit by some metric (e.g. bitscore) and builds score
# matrices for the desired field with queries as rows and databases as
# columns.
#
# Input:
# Argument 1: name of blast report file
# Arguments 2 to n+1: fields to extract
#
# Output:
# For each field, one csv file is written
# 
# Example:
# report2mat.R "~/center/ohome/yyy.csv" bit_score evalue
#
# TODO:
# Currently there is no way for the user to change the metric by which
# the top hit is chosen. Currently the top hit chosen is always the one
# with the highest bit score. 
#
#####################################################################

args <- commandArgs(trailingOnly=TRUE)

stopifnot(length(args) > 1, file.exists(args[1]))

BuildMat <- function(d, field, filename, by="bit_score"){
  # Remove unnecessary data from the input dataframe (which may be a lot)
  rdat <- data.frame(db=d$database, qid=d$query_id, field=d[,field], by=d[,by])
  db.names <- unique(dat$database)
  query.names <- unique(dat$query_id)
  
  # Prepare score matrix output
  ld  <- length(db.names)
  lq  <- length(query.names)
  mat <- matrix(rep(0, ld*lq), nrow=lq)
  rownames(mat) <- query.names
  colnames(mat) <- db.names
  
  for(query in query.names){
    # Create reduced dataframe matching only this query
    qdat <- rdat[which(rdat$qid == query),]
    for(db in db.names){
      # Further reduce the dataframe to include only this database
      qddat <- qdat[which(qdat$db == db),]
      # Find hit that has maximum score (bit_score by default)
      best.match <- which.max(qddat$by)[1]
      mat[query,db] <- qddat$field[best.match]
    }
  }
  write.csv(mat, file=filename, quote=FALSE)
}

dat <- read.csv(file=args[1], header=TRUE, sep="\t")
stopifnot(all(args[2:length(args)] %in% colnames(dat)))

for(field in args[2:length(args)]){
  filename = gsub('.csv', paste0('_', field, '.csv'), args[1])
  BuildMat(dat, field, filename)
}