#!/usr/bin/env perl
use warnings;
use strict;

use Math::Random qw(random_negative_binomial);
use Fasta::Parser;
use Fasta::Seq;

my $target_len = 5000;
my $min_len = 500;
my $len_len = $target_len - $min_len;

my $gap_len = 400;
my $ovl_len = 100;

my $rand_len = len();
my $c = 1;

my $fp = Fasta::Parser->new();

while(my $fa = $fp->next_seq){
    my $rand_sum = 0;
    my $seq_len = length($fa->seq);

    while($rand_sum < $seq_len){
        $rand_len = len();

        if( $rand_sum+$rand_len >= $seq_len){
            print Fasta::Seq->new(
                                  id => "r$c",
                                  seq => substr($fa->seq, $rand_sum)
                                 )->string(80);
            $rand_sum+= $rand_len;
        }else{
            print Fasta::Seq->new(
                                  id => "r$c",
                                  seq => substr($fa->seq, $rand_sum, $rand_len)
                                 )->string(80);
            $rand_sum+=gap();
            $rand_sum+= $rand_len;
        }
        
        $c++;
    }
}


sub len{
    return scalar Math::Random::random_negative_binomial(1000,1,1/$len_len) + $min_len;    
}

sub gap{
    return scalar Math::Random::random_negative_binomial(100,1,1/($gap_len+$ovl_len)) - $ovl_len;
}
