#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use Fasta::Parser;
use Fasta::Seq;

use Kmer;

my $VERSION = 0.1.0;

my $argv = join(" ", @ARGV);

my %def = (
           win_size => 2000,
           win_shift => 1000,
           gc_min => 0,
           gc_max => 1,
           min_report_length => 5000,
          );

my %opt = ();

GetOptions(                     # use %opt (Cfg) as defaults
           \%opt, qw(
                        win_size|win-size=i
                        win_shift|win-shift=i
                        gc_min|gc-min=f
                        gc_max|gc-max=f
                        min_report_length|min-report-length=i
                        version|V!
                        debug|D!
                        help|h!
                   )
          ) or die 'Failed to "GetOptions"';

# help
$opt{help} && pod2usage(1);

# version
if ($opt{version}) {
    print "$VERSION\n"; 
    exit 0;
}

%opt = (%def, %opt);


my $km = Kmer->new(
	kmer_size => $opt{win_size},
	shift_by => $opt{win_shift}
);


my $fp = Fasta::Parser->new();
while(my $fa = $fp->next_seq){
    next unless length($fa->seq) > $opt{win_size};
    my @kmers = $km->kmerize($fa->seq);
    my @gc = map{tr/GCgc// / $opt{win_size}}@kmers;

    my $max = $opt{gc_max};
    my $min = $opt{gc_min};
    my $in_range = 0;
    my @rwin;
    my $f = $opt{win_shift};
    
    for (my $i=0; $i<@gc; $i++) {
        if ($gc[$i] <= $max && $gc[$i] >= $min) {
            unless ($in_range) {
                $in_range = 1;
                push @rwin, [$i * $f];
            }
        } else {
            if ($in_range) {
                $in_range = 0;
                $rwin[$#rwin][1] = ($i * $f) + $opt{win_size};
            }
        }
    }
    $rwin[$#rwin][1] = @gc * $f + $opt{win_size} if $in_range;


    # TODO merge close windows


    foreach (@rwin) {
        # get win gc
        my $wl = $_->[1]-$_->[0];
        next if $wl < $opt{min_report_length};
        my $s = substr($fa->seq, $_->[0], $wl);
        my $gc = $s =~ tr/GCgc//;
        print $fa->id,"\t", $_->[0],"\t", $_->[1], "\t", $gc/$wl,"\n";
    }
}


