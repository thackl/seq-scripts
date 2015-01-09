#!/usr/bin/env perl
use warnings;
use strict;

use Fasta::Parser;
use Fasta::Seq;

use Kmer;

my %opt;
$opt{kmer_size} = 2000;
$opt{kmer_shift} = 1000;
$opt{min_gc} = .45;
$opt{min_stretch_length} = 5000;

my $km = Kmer->new(
	kmer_size => $opt{kmer_size},
	shift_by => $opt{kmer_shift}
);


my $fp = Fasta::Parser->new();
while(my $fa = $fp->next_seq){
    next unless length($fa->seq) > $opt{kmer_size};
    my @kmers = $km->kmerize($fa->seq);
    my $pos = 1;
    my $start;
    for (@kmers){
        if(tr/GC// / $opt{kmer_size} < $opt{min_gc}){
            if($start){
                # just continue
            }
            else{
                $start = $pos;
            }
        }else{
            if($start){
                my $end = $pos-$opt{kmer_shift}+$opt{kmer_size}-1;
                print $fa->id,"\t",$start,"\t",$end,"\n" if $end-$start >= $opt{min_stretch_length};
                $start = undef;
            }else{
                # just do nothing
            }
        }
        $pos+=$opt{kmer_shift};
    }
    #print "\n";
}



