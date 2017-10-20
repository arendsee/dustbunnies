#! /usr/bin/perl

while(<>){
    if(/^>/){
        if(/locus\|(.*?)\|/){
            print "$1, ";
        }
        if(/gi\|(.*?)\|/){
            print "$1\n";
        }
    }
}
