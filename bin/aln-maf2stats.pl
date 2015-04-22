#!/usr/bin/env perl

use warnings;
use strict;

while (<>) {
    chomp;
    my $ah;
    my @as;
    if (substr($_, 0, 1) eq "a") {
        $ah = $_;
        while (<>) {
            last unless substr($_,0,1) eq "s";
            my @r = split(/\s+/, $_);
            push @as, $r[6];
        }
    }
    if ($ah && @as == 2) {
        my @stats = aln2stats(@as);
        print "$ah: @stats\n";
    }
}


sub aln2stats{
    my ($r,$q) = @_;

    # ref gap open/extension
    my $rg = $r =~ tr/-//;
    my $ro = $r;
    $ro =~ tr/-//s;
    my $rgo = $ro =~ tr/-//;
    my $rge = $rg - $rgo;

    # qry gap open/extension
    my $qg = $q =~ tr/-//;
    my $qo = $q;
    $qo =~ tr/-//s;
    my $qgo = $qo =~ tr/-//;
    my $qge = $qg - $qgo;

    # diffs
    my $x = $r ^ $q;
    my $mm = ($x =~ tr/\0//c) - ($rg + $qg);
    my $ma = length($r) - ($rg + $qg + $mm);
    my $idy = $ma/($rg + $qg + $mm + $ma);
    return ($idy, $ma, $mm, $rgo, $rge, $qgo, $qge);
    #return MA*$ma + MM*$mm + RGO*$rgo + RGE*$rge + QGO*$qgo + QGE*$qge;
}
