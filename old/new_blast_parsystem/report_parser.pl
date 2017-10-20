#! /usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;

use vars qw{ $del $evalue $bitscore $identity $by $filename};
my $result = GetOptions (
    'f|file'        => \$filename,
    'd|delimiter=s' => \$del,
    'e|evalue'      => \$evalue,
    's|score'       => \$bitscore,
    'd|identity'    => \$identity,
    'b|by=s'        => \$by,
    );
die "Cannot parse arguments\n" unless $result;
$del = ';' unless $del;
$by = 'bitscore' unless $by;

my $fh = \*STDIN;

# Initialize the Bio::SearchIO object
my $search_obj = new Bio::SearchIO(-format => 'blastxml',
                                   -fh     => $fh);

my $result_ref = &parse_report($search_obj, 
                               \&select_by_score,
                               \&std_csv_output);

foreach my $key (keys %$result_ref){
    print "$key," . $result_ref->{$key}->evalue() . "\n";
}

# TODO pass the result hash to a function that iterates through the
# input file to correctly place each hit

# Returns a hash of top hsp objects for each query
sub parse_report {
    my($search_obj, $selection_code, $output_code) = @_;
    my $result_ref;
    my $is_first = 1;
    
    # Retrieve each Bio::Search::Result::ResultI object
    while(my $result_obj = $search_obj->next_result()){
        my $best_hsp_obj = "";
        if($result_obj->no_hits_found()){
            $output_code->($result_obj,0,$is_first);
            next;
        }

        # Retrieve each Bio::Search::Hit::HitI object
        while(my $hit_obj = $result_obj->next_hit()){
            if(not $best_hsp_obj) {
                $best_hsp_obj = $hit_obj->next_hsp();
            }
            
            # Retrieve each Bio::Search::HSP::HSPI
            while(my $hsp_obj = $hit_obj->next_hsp()){
                $best_hsp_obj = $selection_code->($best_hsp_obj, $hsp_obj);
            }
        }
        $output_code->($result_obj, $best_hsp_obj, $is_first);
        $is_first = 0;
    }
    return $result_ref;
}

sub select_by_score {
    my ($old, $new) = @_;
    my $best = ($old->score() > $new->score()) ? $old : $new;
    return $best;
}

sub std_csv_output {
    my($result, $hsp, $is_first) = @_;
    if($is_first){
        print join($del, '', 'database', 'query', 'evalue', 'bitscore',
              'conserved', 'coverage'), "\n";
    }
    my $database  = $result->database_name();
    my $query     = $result->query_name();
    my $evalue    = ($hsp) ? $hsp->evalue() : 100;
    my $bitscore  = ($hsp) ? $hsp->bits() : 0;
    my $conserved = ($hsp) ? $hsp->frac_conserved('total') : 0;
    my $coverage  = ($hsp) ? $hsp->length('total') / $hsp->length('query') : 0;
    print join($del, $database, $query, $evalue, $bitscore, $conserved,
               $coverage), "\n";
}
