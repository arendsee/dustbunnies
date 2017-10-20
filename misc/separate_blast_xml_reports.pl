#! /usr/bin/perl

####################################################################
#
# This script takes a large file into which many BLAST reports have
# been dumped and writes them to individual files.
#
# It reads from STDIN and writes files to working directory
#
####################################################################

my $fh = \*STDIN;

# First line of first file as divider (should be shared all reports)
my $first_line = <$fh>;

my $current_file = "";
my $current_dbname = "";
while(<$fh>){
    if(/\Q$first_line/){
        if(length $current_file > 0 and length $current_dbname > 0){
            my $filename = $current_dbname;
            $filename =~ s/\..*//;
            $filename =~ s/.*\///;
            $filename .= ".xml";
            open(OUT, '>', $filename) or die $!;
            print OUT $current_file;
            close OUT;
        }
        $current_file = "";
    }
    if(/<BlastOutput_db>(.*?)<\/BlastOutput_db>/){
        $current_dbname = $1;
    }
    $current_file .= $_;
}
