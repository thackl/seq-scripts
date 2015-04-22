#!/usr/bin/env perl

use warnings;
use strict;

use Fasta::Parser;
use Fasta::Seq;

use constant {
    ID => 0,
    OFF => 1,
    LEN => 2,
    STRAND => 3,
    RLEN => 4,
    SEQ => 5
};

my $c=0;
while (<>) {
    chomp;
    next unless $_;
    next if /^#/;
    my ($t, @r) = split(/\s+/, $_);
    if ($t eq "a") {
        $c++;
    }elsif ($t eq "s") {

        my $fa = Fasta::Seq->new(
            id => $r[ID].".$c",
            seq => $r[SEQ],
            desc => "OFF:".$r[OFF]." LEN:".$r[LEN].
                " STRAND:".$r[STRAND]." RLEN:".$r[RLEN]
            );
        $fa->{seq} =~ tr/-//d;
        print $fa->string(80);
    }
}
