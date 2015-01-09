use warnings;
use strict;
  
use Fastq::Parser;
use Fastq::Seq;
  
use File::Basename;

die "3 args required" unless @ARGV == 3;

my ($ls, $f1, $f2) = @ARGV;
my @l=sort{$b<=>$a}split(",", $ls);  
my $fp1=Fastq::Parser->new(file=>$f1); 
my $fp2=Fastq::Parser->new(file=>$f2); 
 
my $op1=basename($f1,".fq");
my $op2=basename($f2,".fq");
 
my @o1;
my @o2;
my @lrev=reverse(@l);
for(my $k=0; $k<@lrev; $k++){
    my $l=$lrev[$k];
    my $to=$k < $#lrev ? "-".($lrev[$k+1]-1) : "";
    open($o1[$l], ">", $op1.".".$l.$to.".fq") or die $!;
    open($o2[$l], ">", $op2.".".$l.$to.".fq") or die $!;
}

my ($fq1, $fq2);
while( ($fq1=$fp1->next_seq) && ($fq2=$fp2->next_seq) ){
    foreach my $l (@l){
	if(length($fq2->seq) >= $l && length($fq1->seq) >= $l){
	    print {$o1[$l]} "$fq1";
	    print {$o2[$l]} "$fq2";
	    last;
	}
    }
};

close $_ for grep{defined}(@o1, @o2);
