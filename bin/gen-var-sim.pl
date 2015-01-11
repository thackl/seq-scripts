#!/usr/bin/env perl
use warnings;
use strict;

use Data::Dumper;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use Getopt::Long;
use Pod::Usage;

use Fasta::Parser;
use Fasta::Seq;

use Cfg;


=head1 SYNOPSIS

  gen-var-sim.pl --ref genome.fa --sdi *.sdi --var 200 --out genome-v200.fa

=head1 OPTIONS

=over

=item --ref <FA>

Reference genome file in FASTA format.

=item --sdi <SDI>

Variant annotations in SDI format, sorted by coordinates.

Each line of an sdi file consists of the columns chromosome, position, length,
reference base, consensus base.  chromosome: The same name as the chromosome id
in reference fasta file position: 1-based leftmost position length: the length
difference of the changed sequence against reference (0 for SNPs, negative for
deletions, positive for insertions) reference base: [-A-Z]+ (regular expression
range) consensus base: [-A-Z]+ ((regular expression range), IUPAC code is used
for heterozygous sites detection score: (1 detected by both IMR and DENOM, 2
detected by IMR as heterozygous SNP but detected by DENOM as homozygous SNP, 4:
detected by IMR only, 5: detected by DENOM only) HMQ coverage: the coverage from
reads with mapq >=30 SNP Phred score HMQ consensus base: the consensus base when
considering reads with high mapping quality (mapq >= 30)

=item [--out <FA>] [STDOUT]

Output file for generated sequences.

=item [--var <INT>] [MAX]

Variance rate as 1/--var. Defaults to max available rate based on given number
of available variance annotations.

=item [--debug]

Enable debugging messages.

=back

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut

# =item --snp <TSV> [<TSV ...>]
# 
# Annotation of SNPs as TSV:
# 
#   Sample  Chromosome  Position  Reference allele  New allele  [Annotation]
# 
# =item --del <TSV> [<TSV ...>]
# 
# Annotation of deletions as TSV:
# 
#   Sample  Chromosome  Begin  End  [Deletion]  [Annotation]
# 
# =item --ins <TSV> [<TSV ...>]
# 
# Annotation of insertions as TSV:
# 
#   Sample  Chromosome  Begin  End  Insertion  [Annotation]
#

## test call
# ../bin/gen-var-sim.pl --ref TAIR10_chr1.h1M.fa --sdi *.sdi --var 200 --out at-var.fa

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init(\<<'CFG');
	log4perl.rootLogger                 = DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [gvs] %m%n
CFG
$L->level($INFO);



my ($ref_fa,$out_fa,$var,$version,$debug,$help,$man) = ("","",0,0,0,0,0);
my @var_sdis;

GetOptions(
           'ref=s'=>\$ref_fa,
           'sdi=s{1,}'=>\@var_sdis,
           'out=s'=>\$out_fa,
           'var=i'=>\$var,
           'version=i'=>\$version,
           'debug!'=>\$debug,
           'help|?'=>\$help ,
           'man'=>\$man
          ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
$debug && $L->level($DEBUG);

my $faf = $ARGV[0];

# required
$ref_fa || pod2usage("--ref required");
@var_sdis || pod2usage("--sdi required");
foreach ($ref_fa, @var_sdis){$L->logdie("Cannot find $_") unless -e $_};



##----------------------------------------------------------------------------##
# analyse input

$L->info("Indexing reference: $ref_fa");

Cfg->Check_binaries('samtools');
my $sam_idx_log = qx(samtools faidx $ref_fa);
$L->info($sam_idx_log) if $sam_idx_log;
$? && $L->logdie($?);

open(FAI, $ref_fa.".fai") or $L->logdie($!);
my $ref_len = 0;
my %ref_ids;
while(my $fai = <FAI>){
    my ($id, $len) = split("\t", $fai);
    $ref_ids{$id} = 0;
    $ref_len+= $len;
}
close FAI;

$L->debug("Reference length: $ref_len");


my $var_sum = qx(wc -l @var_sdis);
$L->info("Available variances:\n$var_sum");

my ($var_tot) = $var_sum =~ /(\d+)[^\n]+$/;
$L->info("Number of available variant annotations: $var_tot");
my $var_max = sprintf("%.2f", $ref_len/$var_tot);
$L->info("Maximum possible variant rate: 1/". $var_max
         ."\n Please keep in mind that the effecitve rate of introduced SNPs\n"
         ." will be lower if annotations allow identical substitutions\n"
         ." e.g. A->A or A->R (iupac AG) or ar applied to identical positions");
my $var_sample;

if($var){
    if($var < $var_max){
        $L->logdie("1/$var variants requested, but at max 1/$var_max variant annotations available");
    }
    $var_sample = $var_max/$var;
    $L->info("Sampling ".sprintf("%.2f", $var_sample*100)."% of variants to result in estimated rate of: 1/". $var);
}


if($out_fa){
    open(FAO, '>', $out_fa) or $L->logdie($!);
}else{
    open(FAO, ">&=STDOUT") or $L->logdie($!);
}

##----------------------------------------------------------------------------##
# preping some globals

my $fa;
my $c=0;   # var counter
my $adj=0; # adjust coordinates by previously introduced indels
my $d;     # debug msg
my %fa_ids;

my %deiupac = load_deiupac();

use constant {
    SID => 0,
    POS => 1,
    LEN => 2,
    RBP => 3,
    CBP => 4
};

# SNP: LEN = 0
# INS: LEN > 0
# DEL: LEN < 0



##----------------------------------------------------------------------------##
# looping variances

open(SDI, "sort -m -k1,1 -k2,2g @var_sdis |") or $L->logdie($!);

my $var_10 = int($var_sample * $var_tot / 10);
my $var_p = 0;

$L->info("Mutating");
$|++;
printf ("%-s%% ", 0) if $L->level <= $INFO;


# Loop over the input stream
while (my $var = <SDI>) {
    # sampling
    if($var && (my $var_rand = rand(1)) > $var_sample){
        $L->debug("Skipping $var_rand\n");
        next;
    }
    
    my @var = split("\t", $var);
    # retrieve current seq
    if (!$fa || $fa->id ne $var[SID]) {
        if($fa){
            $ref_ids{$fa->id}++;
            print FAO $fa->string(80);
        }
        $fa = get_seq($var[SID]);
        $adj = 0; # reset adj if seq is loaded firston every new seq
    }

    my $pos = $var[POS]-1+$adj;      # sdi is one-based

    if ($var[LEN] == 0) { # snp
        my $b = substr($fa->seq, $pos, 1);

        if ($L->level <= $DEBUG) {
            $d = "SNP: $pos+$adj\n";
            $d.= " "x 5 .".\n";
            $d.= substr($fa->seq, $pos-5, 11)."\n";
        }

        substr($fa->{seq}, $pos, 1, deiupac($var[CBP]));

        if ($L->level >= $DEBUG){
            $d.= substr($fa->seq, $pos-5, 11)."\n";
        }

    } elsif ($var[LEN] > 0) { # ins
        if ($L->level <= $DEBUG) {
            $d = "INS: $pos+$adj\n";
            $d.= " "x 5 .">\n";
            $d.= substr($fa->seq, $pos-5, 5). " "x $var[LEN]. substr($fa->seq, $pos, 5)."\n";
        }

        substr($fa->{seq}, $pos, 0, deiupac($var[CBP]));
        $adj+=$var[LEN];

        if ($L->level >= $DEBUG){
            $d.= substr($fa->seq, $pos-5, 10+ $var[LEN])."\n";
        }
        
    } elsif ($var[LEN] < 0) { # del
        if ($L->level >= $DEBUG) {
            $d = "INS: $pos+$adj\n";
            $d.= " "x 5 ."<\n";
            $d.= substr($fa->seq, $pos-5, 10+abs($var[LEN]))."\n";
        }

        substr($fa->{seq}, $pos, abs($var[LEN]), "");
        $adj+=$var[LEN];

        if ($L->level <= $DEBUG){
            $d.= substr($fa->seq, $pos-5, 5). " "x abs($var[LEN]). substr($fa->seq, $pos, 5);
        }
    } else {
        $L->info("Check sdi file format: expecting length in 3 column");
    }

    $L->debug("@var"."$d\n");

    if($L->level <= $INFO){
        unless($c++%$var_10){
            $var_p+=10;
            printf ("%-s%% ", $var_p);
        }
    }
}
print "\n" if $L->level <= $INFO;
$|--;

close SDI;

$L->info("Finishing of output");

if($fa){
    $ref_ids{$fa->id}++;
    print FAO $fa->string(80);
}


# output non-modified refseqs
foreach my $rid (sort keys %ref_ids){
    next if $ref_ids{$rid}; # modified and outputted
    my $fa = get_seq($rid);
    print FAO $fa->string(80);
}

close FAO;



##----------------------------------------------------------------------------##

sub get_seq{
    my ($id) = @_;
    open(my $fai, "samtools faidx $ref_fa $id |") or $L->logdie;
    my $fa = Fasta::Parser->new(fh => $fai)->next_seq;
    if(! length($fa->seq)){
        $L->warn("Couldn't retrieve sequence >$id from FASTA");
        next;
    }
    close $fai;
    return $fa;
}


sub deiupac{
    join("", map{
        $deiupac{$_}()
    }split(//, $_[0]))
}


sub load_deiupac{
    return (
            A => sub{return "A"},
            T => sub{return "T"},
            C => sub{return "C"},
            G => sub{return "G"},
            R => sub{return ("A","G")[int(rand(2))]},            # A or G
            Y => sub{return ("C","T")[int(rand(2))]},            # C or T
            S => sub{return ("G","C")[int(rand(2))]},            # G or C
            W => sub{return ("A","T")[int(rand(2))]},            # A or T
            K => sub{return ("G","T")[int(rand(2))]},            # G or T
            M => sub{return ("A","C")[int(rand(2))]},            # A or C
            B => sub{return ("C","G","T")[int(rand(3))]},        # C or G or T
            D => sub{return ("A","G","T")[int(rand(3))]},        # A or G or T
            H => sub{return ("A","C","T")[int(rand(3))]},        # A or C or T
            V => sub{return ("A","C","G")[int(rand(3))]},        # A or C or G
            N => sub{return ("A","C","G","T")[int(rand(4))]},    # any base
            a => sub{return "a"},
            t => sub{return "t"},
            c => sub{return "c"},
            g => sub{return "g"},
            r => sub{return ("a","g")[int(rand(2))]},            # a or g
            y => sub{return ("c","t")[int(rand(2))]},            # c or t
            s => sub{return ("g","c")[int(rand(2))]},            # g or c
            w => sub{return ("a","t")[int(rand(2))]},            # a or t
            k => sub{return ("g","t")[int(rand(2))]},            # g or t
            m => sub{return ("a","c")[int(rand(2))]},            # a or c
            b => sub{return ("c","g","t")[int(rand(3))]},        # c or g or t
            d => sub{return ("a","g","t")[int(rand(3))]},        # a or g or t
            h => sub{return ("a","c","t")[int(rand(3))]},        # a or c or t
            v => sub{return ("a","c","g")[int(rand(3))]},        # a or c or g
            n => sub{return ("a","c","g","t")[int(rand(4))]},    # any base
           );
}
