#! /usr/bin/perl

=head1 NAME

=head1 DESCRIPTION

Removes internal headers, ensures all headers are equivalent, ensures all
rows have the correct number of terms.

=cut


my $fh = \*STDIN;

my $first_header = <$fh>;
my $del = substr $first_header, 0, 1;
(my $el = $first_header) =~ s/[^\Q$del\E]//g;
my $nterms = length($el) - 1;

print $first_header;

my $line_number = 1;
while(<$fh>){
    if(/^$del/){
        print STDERR "Contradictory headers, dying.\n" if not $_ eq $first_header;
        next;
    }
    print $_;
    (my $el = $_) =~ s/[^\Q$del\E]//g;
    die "Unequal term counts, dying\n" if length($el) != $nterms;
}
