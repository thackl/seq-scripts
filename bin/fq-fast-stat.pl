#!/usr/bin/env perl
use warnings;
use strict;

my $sc=0;
my $bp=0;
my $gc=0;
my $nc=0;
my %ls;

my $c=0;
while (<>) {
    $c++;
    if ($c == 2 ){
        chomp();
        $sc++;
        $bp+=length($_);
        $ls{length($_)}++;
        $gc+= ($_ =~ tr/GCgc//);
        $nc+= ($_ =~ tr/Nn//);
    }
    $c=0 if $c>3;
}

my $gcp = $gc/$bp*100;
my $ncp = $nc/$bp*100;

my @ls = sort{$a<=>$b}keys %ls;

foreach (@ls) {
    $ls{$_} = 0 unless exists $ls{$_};
}

my $lmin = $ls[0];
my $lmax = $ls[-1];
my $n50;
my $nbp;

for (my $i=$#ls; $i>=0; $i--) {
    my $l = $ls[$i];
    $nbp+= $ls{$l} * $l;
    if ( $nbp >= $bp/2) {
        $n50 = $l;
        last;
    }
}

print "#seqs\ttotal_bp\tmin_bp\tmax_bp\tn50_bp\tGC_%\tN_%\n";
printf "%s\t%s\t%s\t%s\t%s\t%0.2f\t%0.2f\n", $sc, $bp, $lmin, $lmax, $n50, $gcp, $ncp;

# my $bs = 10**(length(scalar @ls)-1);
#if ($bs == 1) {
    for (my $i=$#ls; $i>=0; $i--) {
        my $l = $ls[$i];
        print "##\t$l\t$ls{$l}\t" ,$ls{$l} * $l,"\n";
    }
#}else {
#    warn "binning not yet implemented\n";
#}

