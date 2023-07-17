#!/usr/bin/perl -w

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'fastafile|f=s','numseq|n=s','check|c=s','outseqname|o=s','help|h');

$/ = "\n>";
my $fasta = $opt{fastafile};
my $numseq = $opt{numseq};
open IN, "<$fasta" or die $!;
my $i = 0;
my $numfile = 0;
while (my $record = <IN>) {
    chomp($record);
    $record =~ s/\r\n*/\n/g;
    $record =~ s/^>//;
    $i ++;
    my ($header,@seq) = split(/\n/,$record);
    $seq = join("",@seq);
    $seq =~ s/\s+//;
    $seqname = ((split/\s+/,$header))[0];
    if ($numseq == 1 || $i % $numseq == 1) {
	close OUT;
	$numfile ++;
	$filename = "outseq_$numfile.fasta";
	if ($opt{outseqname}) {
	    $filename="$seqname.fasta";
	}
	open OUT, ">$filename" or die $!
    }
    if ($check && $check == 'NT') {
	next unless (uc($seq) =~ m/[AGTCN]+/);
    }else {
	next unless (uc($seq) =~ m/\w+/);
    }
    print OUT ">",$header,"\n",join("",@seq),"\n";
}
