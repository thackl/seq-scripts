#!/usr/bin/env perl
use warnings;
use strict;

use Data::Dumper;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use Getopt::Long;
use Pod::Usage;

use Bio::Tools::GFF;
use Bio::FeatureIO::gff;
use Bio::Tools::CodonTable;
use Bio::SeqIO;

use Fasta::Parser;
use Fasta::Seq;

use Cfg;

=head1 OPTIONS

=over

=item --gff

Input gff file

=item [--version] [3]

Gff version.

=item [--debug]

Enable debugging messages.

=back

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut


# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init(\<<'CFG');
	log4perl.rootLogger                 = DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [dsg] %m%n
CFG
$L->level($INFO);



my ($gff,$version,$type,$source,$aa,$debug,$help,$man) = ("",3,"","",0,0,0);

GetOptions(
           'gff=s'=>\$gff,
           'version=i'=>\$version,
           'type=s'=>\$type ,
           'source=s'=>\$source,
           'aa=s'=>\$aa,
           'debug!'=>\$debug,
           'help|?'=>\$help ,
           'man'=>\$man
          ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
$debug && $L->level($DEBUG);

my $faf = $ARGV[0];

# required
$faf || pod2usage("required: FASTA");
$gff || pod2usage("required: --gff");
# $type || pod2usage("required: --type");
# $source || pod2usage("required: --source");

Cfg->Check_binaries('samtools');

#
$L->info( qx(samtools faidx $faf) );
$? && $L->logdie($?);

# Initialise input stream
my $gffio = Bio::Tools::GFF->new(-file => $gff, -gff_version => $version) or die $!;

my $fa;
my $c=0;

# Loop over the input stream
while(my $feature = $gffio->next_feature())
{
    my $id = $feature->seq_id;
    # retrieve current seq
    if(!$fa || $fa->id ne $id){
        $fa = get_seq($id);
    }
    
    my $mod =  $feature->primary_tag();
    my $pos = $feature->location->start;
    if($mod eq "substitution"){
        my $b = substr($fa->seq, $pos, 1);
        $L->debug(" "x5,".");
        $L->debug(substr($fa->seq, $pos-5, 11));
        my $s = (grep{$_ ne $b}qw(A T G C))[int(rand(3))];
        substr($fa->{seq}, $pos, 1, $s);
        $L->debug(substr($fa->seq, $pos-5, 11));
    }elsif($mod eq "insertion"){
        $L->debug(" "x5,">");
        $L->debug(substr($fa->seq, $pos-5, 5), " ", substr($fa->seq, $pos, 5));
        my $s = (qw(A T G C))[int(rand(4))];
        substr($fa->{seq}, $pos, 1, $s);
        $L->debug(substr($fa->seq, $pos-5, 11));
    }elsif($mod eq "deletion"){
        $L->debug(" "x5,"<");
        $L->debug(substr($fa->seq, $pos-5, 11));
        substr($fa->{seq}, $pos, 1, "");
        $L->debug(substr($fa->seq, $pos-5, 5), " ", substr($fa->seq, $pos, 5));
    }else{
        $L->info("Skipping unknown modification: $mod");
    }

    
    #print Dumper($fa->id, $feature->seq_id, $primary_tag, $source_tag, $feature->location->start, $feature->location->end);	
    exit if $c++ >100;
}

# Close input stream
$gffio->close() or die $!;



sub get_seq{
    my ($id) = @_;
    open(my $fai, "samtools faidx $faf $id |") or $L->logdie;
    my $fa = Fasta::Parser->new(fh => $fai)->next_seq;
    if(! length($fa->seq)){
        $L->warn("Couldn't retrieve sequence >$id from FASTA");
        next;
    }
    close $fai;
    return $fa;
}













