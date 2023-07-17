#!/usr/bin/perl -w
################################################################################
#
#   Name:    SCIPIO - Eukaryotic Gene Identification
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2007-12-17
#   Version: 1.0
#
#
#   scipio.pl: [<options>] <target> <query> 
#
#
#   This script runs BLAT and, for each queried protein sequence, tries to find
#   *the* (correct) location. 
#   - if BLAT returns more than one location, take the highest score. 
#   - if BLAT only returns partial matches, try to assemble them by matching to
#   the end of each target.
#   
#   Also, this script tries to locate codons that are split by an intron which
#   BLAT does not do. 
#
#   The script produces output in YAML format (see documentation)
#   
#   For usage information, run scipio.pl without parameters
#
#   Counting of bps/aas: interval counting (see documentation)
#


use strict;
use File::stat;
use Getopt::Long;
use List::Util('sum', 'max', 'min');
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use YAML;
if ($YAML::VERSION gt "0.60") {
    use YAML::Dumper::Base;
}
 

my $QPFX="prot";
my $TPFX="dna";
my $TRPFX="trans";
my $QNPFX="nucl";


################################## user options ################################

my $BLAT_BIN = "blat";      # name of blat executable, dont run if empty
my $BLAT_OUTPUT = "";       # name of blat  output file
my $BLAT_TILESIZE;
my $BLAT_MIN_IDENTITY;      # minimal hit accuracy in %
my $BLAT_MIN_SCORE = 15;
my $BLAT_PARAMS = "";       # parameters passed to BLAT directly
my $FORCE_NEW = 0;          # overwrite previous output file, always run blat
my $TR_TABLE = 1;
my $VERBOSE = 0;
# my $JOIN_FS = 1;
my $SPLIT_ON_FS = 0;
# my $SPLIT_STRING = "_split:";  # split mode disabled
my $MIN_IDENTITY = 0.9;
my $MIN_BEST_SCORE = 0.3;
my $MAX_MISMATCH = 0;
my $ADD_MARGINS = 1;
my $MAX_CHECK_OVERLAP = 300;
my $DEFAULT_REGION_SIZE = 1000;
my $REGION_SIZE = $DEFAULT_REGION_SIZE;     # this much is shown in log
my $SKIP_MISMATCHES = 1;
my $SKIP_MISRATIO = 0.01;
my $SHOWALL_SCORE = 1.0;
my $ALLOW_MISSING_TARGETS = 0;

my $HIT_BLESS_LIST =
    "ID,status,reason,${QPFX}_len,${QPFX}_start,${QPFX}_end,${QPFX}_seq,".
    "target,target_len,strand,${TPFX}_start,${TPFX}_end,matches,mismatches,".
    "undetermined,unmatched,additional,score,upstream,upstream_gap,matchings,stopcodon,downstream,downstream_gap";
my %MATCH_BLESS_LISTS = ( 
    exon =>   
    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq,".
    "seqshifts,mismatchlist,undeterminedlist,inframe_stopcodons,translation,overlap",
    intron => "type,${QNPFX}_start,${TPFX}_start,${TPFX}_end,seq",
    gap =>    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq",
    any =>    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq");

my $HIDE_UNDEF = 0;         # show also undefined keys and empty lists
my $HIDE_DEFAULTS = 0;      # show also keys with default values
my $DONT_DIE = 1;
my $SHOW_BLATLINE = 0;

my %PARAMETER = ("blat_output=s" => \$BLAT_OUTPUT,
		 "force_new" => \$FORCE_NEW,
		 "verbose" => \$VERBOSE,
		 "best_size|min_score=f" => \$MIN_BEST_SCORE,
		 "min_identity=f" => \$MIN_IDENTITY,
		 "max_mismatch=i" => \$MAX_MISMATCH,
		 "transtable=i" => \$TR_TABLE,   # undocumented until BLAT supports this
		 "margins!" => \$ADD_MARGINS,    # undocumented; only needed to show feature
		 "split_on_fs" => \$SPLIT_ON_FS,
		 "region_size=i" => \$REGION_SIZE,
		 "force_score=f" => \$SHOWALL_SCORE,
		 "partial_target" => \$ALLOW_MISSING_TARGETS,
		 "show=s" => \$HIT_BLESS_LIST,
		 (map { ("show_${_}=s" => \$MATCH_BLESS_LISTS{$_}) } keys %MATCH_BLESS_LISTS),
		 "hide_undef" => \$HIDE_UNDEF,
		 "hide_defaults" => \$HIDE_DEFAULTS,

		 "blat_bin=s" => \$BLAT_BIN,
		 "blat_params=s" => \$BLAT_PARAMS,
		 "blat_tilesize=i" => \$BLAT_TILESIZE,
		 "blat_score=i" => \$BLAT_MIN_SCORE,
		 "blat_identity=f" => \$BLAT_MIN_IDENTITY,
#		 "continue_on_bug" => \$DONT_DIE,    # undocumented  
		 "show_blatline" => \$SHOW_BLATLINE  # undocumented
);

my $PRED_OVERLAP=2;    # how many codons of prediction are checked for intron start/end
my $GFF_REVERSE=1;
my $MAX_OVERLAP=5;     # how much overlap we allow for hits on different targets
my $MIN_REMAIN=4;      # how many bases have to be left after intron offsets
my $MAX_SHIFTLEN = 6;  # maximal length (in bases) of non GT..AG gap (-> seqshift)

my $MAX_GAPLEN = 6;    # maximal size (in amino acids) of gap in query that we try to close
my $HITBUFFER_SIZE = 75000;  # this much is loaded around a blat hit; 
                             # also the maximum length of an intron across different targets
my $INTERNAL_MIN_MATCH = 15;
my $IFS_MALUS = 1.5;

local $YAML::InlineSeries = 10;

####################### function headers #######################################

sub usage
{
    my $s = shift;
    print STDERR "\n $s" if (defined $s);
    print STDERR "
    usage:
    scipio.pl [<options>] <target> <query> 
   
    <target> is a DNA file
    <query> is a protein file
    both in FASTA format

    Options:    --blat_output=<filename>   name of BLAT output file (*.psl); if not specified,
                                           defaults to the newest *.psl file in the working
                                           directory, or to \"Scipio<nnnn>_blat.psl\" if none
                                           exists
                --force_new                run BLAT even if output file already exists; 
                                           (the default is never to run BLAT if the specified
                                           or default *.psl file would be overwritten)
                --verbose                  show verbose information (including progress)
                --min_score=<value>        minimal size of the largest partial hit (as
                                           a fraction of the total sequence)
                --min_identity=<value>     minimal identity in any hit (default is 0.9)
                --max_mismatch=<value>     maximal number of mismatches in any hit
                                           (0 means any number, the default)
                --split_on_fs              do not join matchings separated by sequence shifts 
                                           (the default is to consider them a single exon)
		--region_size=<value>      size of shown upstream/downstream sequence parts
                                           (defaults to $DEFAULT_REGION_SIZE, maximum is $HITBUFFER_SIZE)
		--force_score=<value>      specify a score (between 0 and 1) that forces a hit
                                           to be shown even if it is contradicting a better one
                                           (these hits are not assembled together and shown as 
                                           <queryname>_(1), etc.)
                --partial_target           accept BLAT output files containing hits referring 
                                           to nonexistent target sequences
                --show=<list>              which keys are to be shown in YAML output
                --show_intron=<list>       (for details see documentation)
                --show_exon=<list>
                --show_gap=<list>
                --hide_undef               hide keys with undefined or empty list values
                --hide_defaults            hide some keys that have default values
 
    Options passed to BLAT (ignored when BLAT is not run):
		--blat_bin=<name>          name of BLAT executable, defaults to \"blat\"
                --blat_params=<params>     see BLAT usage information
                --blat_tilesize=<..>
                --blat_score=<..> 
                --blat_identity=<..>       by default, this is set to 90% of min_identity


";
    exit 1;
}

sub run_bug;                # die on bug; continue if --continue_on_bug is set
sub revcom;                 # returns reverse complement of $1
# sub extract_offset;         # if $1 is <target>:<offset>, strips and returns offset, otherwise 0
sub with_number;            # returns ($1 $2), and adds "s" if $1 is a number not equal to 1
sub diff_str;               # compares strings $1, $2, in $3 positions, starting with position $4 aligned to $5
                            # returns "|" for matches, "X" for mismatches, and " " if "X" or "-" is found
sub str_diff_list;          # returns a list of all positions where diff_str finds mismatches
sub str_diff_count_all;        # counts the "X"s and spaces in diff_str's result
sub str_diff_count_undeterm;   # counts the spaces in diff_str's result
sub str_diff_count_mismatches; # counts the "X"s in diff_str's result
sub str_pos;                # return the occurences of pattern $1 in str

sub extend_to_full_name;    # returns full name of a target (the %full_name value given the cut targetname as key)
sub save_target_region;     # marks a region around a BLAT hit. ($1,$2,$3) = targetname, from, to
sub cut_target_regions;     # cuts regions to proposed length 
sub fill_target_region;     # fills a region with the sequence. ($1,$2,$3) = target, ref to seq, offset
sub read_targets;           # read target sequences from FASTA file, fills marked regions
sub validate_targets;       # removes empty regions, check if all sequences have been found
sub nu_to_aa;               # division by 3; rounded

# revtranslate_to_pattern ($seq)
# return a regular expression base on reverse translation of the amino acid sequence $seq.
sub revtranslate_to_pattern; 

# split_codons($aa, $dna, $backward)
# return a list of patterns modelling a splice site (dss if $backward is true, ass otherwise)
# $aa contains one amino acid, $dna the other splice site (which is known), together with 
# one half (possibly empty) of a codon. The returned list contains all patterns that complete
# the known splice site to a codon coding for $aa
sub split_codons;

# get_splicesite_patterns_forward($suffix, $dna)
# take the first two amino acids of $suffix; each of them is split with split_codons
# $dna is representing the dna sequence starting with the first possible splice point
# return a pattern that matches all possible ass's that complete $dna to a valid intron
sub get_splicesite_patterns_forward;
sub get_splicesite_patterns_backward;

# subseq ($targetname, $from, $count, $complement)
# returns subsequence of a target; 
# if $complement is true, $from must be negative to be inside the target, and -$from 
# refers to the subsequence end (so $count characters are taken beginning with -$from-$count);
# in this case the reverse complement is returned; unknown bases are given as ".", whereas
# positions outside the target are returned as spaces. It is guaranteed that the length
# of the result equals $count. If $count is negative, undef is returned.
sub subseq;    
  
# get_intron_offset ($gap_aa, $dss_seq, $ass_seq, $frameshift_mode, preferred_offset) 
# this routine concatenates a prefix of $dss_seq to a suffix of $ass_seq, where
# $dss_seq and $ass_seq are nucleotide sequences with length of 3*length($gap_aa)+2,
# and tries to align the resulting translated sequence to $gap_aa.
# Returned as offset is the length of the left part (the prefix). In frameshift mode,
# only multiples of three can be taken as offsets.
# If multiple offsets yield the minimum number of errors, then the one closest to
# preferred_offset is taken
#
# it returns a list with four values: 
#  - best offset, 
#  - best number of mismatches
#  - translated sequence
#  - did we find a splice site pattern
#
# if frameshift_mode is false, then splice side patterns are traded off against more
# errors (INTRON mode); in frameshift mode, splice side patterns are ignored
# 
sub get_intron_offset;

# calculate_diff_strings ($targetname, $queryname, 
#                         \@targetlocations, \@querylocations, $complement);
# returns a list of mismatch strings, one for each exon
sub calculate_diff_strings;

# check if successive BLAT hits are compatible to each other
sub incompatible;

# replace elements x with (from < x <= to) from list
# by elements found in $newelems
# returns the change in number of elements
sub replace_in_list;

# replace_in_last( $nexthit, $dnaseq, $protseq, $dna_add, $nsub)
# replace_in_next( $nexthit, $dnaseq, $protseq, $dna_add, $nsub)
# adds a sequence part to the first/last exon of a hit
#
# the exon matching is changed the following way:
# - the translation is set to the translation of $dna_add.$dnaseq ($dnaseq is
#   supposed to start at frame length($dna_add); it usually ends in frame 0;
#   if the translation is longer then $protseq, only the suffix of the corresponding length 
#   will be compared to $protseq; here, if dnaseq doesnt start in frame 0, 
#   $dna_add *must* be defined to match the reading frame.
# - the mismatch positions between $protseq and the translation are added to
#   the first matching
# - translation is concatenated with first matching
# - prot_start, nucl_start and dna_start are brought forward by the length of $protseq
sub replace_in_last;
sub replace_in_next;

# add initial ($hit, $type, $qlen)
# inserts an initial matching of type $type to the hit
sub add_initial;

# add_final ($hit, $type, $qlen)
# inserts a terminal matching of type $type to the hit
sub add_final;

# set_problem_field ($hit)
# sets the problem and reason keys according to the condition of the hit
sub set_problem_field;

##################################### main #####################################

my %dbase;
my %queries;
my %full_name;
my $ttable = Bio::Tools::CodonTable->new( -id => $TR_TABLE );

### Command Line
&GetOptions (%PARAMETER) or &usage;

# prepare verbose output
my $comments;
if ($VERBOSE) {
    open COMMENT, ">&", \*STDERR;
    select COMMENT;
    $| = 1;
    select STDOUT; 
} else {
    open COMMENT, '>', \$comments;
}

$REGION_SIZE=$HITBUFFER_SIZE if ($REGION_SIZE > $HITBUFFER_SIZE);
$REGION_SIZE=0 if ($REGION_SIZE < 0);

# show parameters in verbose mode
while (my ($key, $ref) = each %PARAMETER) {
    $key =~ s/[!:=].*$//;
    print COMMENT "$key=$$ref\n" if (defined $$ref);
}
print COMMENT "\n";
# parameters passed to BLAT
$BLAT_MIN_IDENTITY = 90 * $MIN_IDENTITY unless (defined $BLAT_MIN_IDENTITY);
$BLAT_PARAMS.=" -tileSize=$BLAT_TILESIZE" if (defined $BLAT_TILESIZE && $BLAT_PARAMS !~ /tileSize/);
$BLAT_PARAMS.=" -minIdentity=$BLAT_MIN_IDENTITY" unless ($BLAT_PARAMS =~ /minIdentity/);
$BLAT_PARAMS.=" -minScore=$BLAT_MIN_SCORE" unless ($BLAT_PARAMS =~ /minScore/);

&usage unless(@ARGV);

my ($dbasefile, $queryfile) = @ARGV;

&usage ("Please specify target file in FASTA format!\n") unless ($dbasefile && -f $dbasefile);
&usage ("Please specify query file in FASTA format!\n") unless ($queryfile && -f $queryfile);


### Open query file and store sequences in hash
my $querySeq = Bio::SeqIO->new(-file=>$queryfile, -format=>"fasta") 
    or usage("Queryfile did not contain sequences.\n");
my $unique_queries=1;
while (my $seq = $querySeq->next_seq()) {
    my ($newid, $newseq) = ($seq->id(), $seq->seq());
    $newseq =~ s/[-]//g;
    my $oldseq = $queries{$newid};
    if (defined $oldseq && $oldseq ne $newseq) {
	$unique_queries=0;
	my $len = length($newseq);
	$newid = "${newid}_[$len]";
	if (length($oldseq) == $len || 
	    defined $queries{$newid} && $queries{$newid} ne $newseq) {
	    print STDERR "Neither query sequence names nor sequence lengths unique. Giving up.\n";
	    exit 1;
	}
    } 
    $queries{$newid} = $newseq;
}

print STDERR "Warning: Query sequence names not unique. Trying to identify by sequence length...\n"
    unless ($unique_queries);

### Determine BLAT output file and run BLAT
unless ($FORCE_NEW) {
    unless ($BLAT_OUTPUT) {
	# choose newest file if none is specified
	# sort by mtime, newest first
	my @list = sort { (stat $b)->mtime <=> (stat $a)->mtime } <*.psl>;
	$BLAT_OUTPUT=$list[0];
    }
    if ($BLAT_OUTPUT && -f $BLAT_OUTPUT) {
	print COMMENT "Not running BLAT. Using previously created output file \"$BLAT_OUTPUT\".\n";
    } else {
	$FORCE_NEW = 1;
    }
}
if ($FORCE_NEW) {
    # do run BLAT (output file not specified, or not existing)
    unless ($BLAT_OUTPUT) {
	my ($n,$next)=("","a"); 
	do {  
	      $BLAT_OUTPUT="Scipio${$}${n}_blat.psl";
	      $n=$next; $next++;
	} while (-f $BLAT_OUTPUT);
    }
    my $showtime = time(); 
    print COMMENT 
	"Running BLAT to produce file \"$BLAT_OUTPUT\": (started at ".
	join (":", map { sprintf "%02d",$_ } (localtime($showtime))[2,1,0]).")...\n";
    my $BLAT_EXESTRING = "$BLAT_BIN -t=dnax -q=prot -noHead $BLAT_PARAMS \"$dbasefile\" \"$queryfile\" \"$BLAT_OUTPUT\"";
    print STDERR "$BLAT_EXESTRING\n" if ($SHOW_BLATLINE);
    open (BLRUN, "$BLAT_EXESTRING |");

    print COMMENT while (<BLRUN>);
    my $blat_exitstatus = $? >> 8;
    printf STDERR "Warning: blat returned exit status %d\n", $blat_exitstatus 
	if ($blat_exitstatus);
    $showtime = time() - $showtime;
    print COMMENT "...done in $showtime seconds\n";
}
open BLOUT, $BLAT_OUTPUT or die ("Could not access output file \"$BLAT_OUTPUT\"\n");

my $hitcount=0;
### Open target file and load sequence parts around hits
while (<BLOUT>)
{
    next if (/^\#/ || /^$/);
    my @result = split "\t";
    next unless @result>=21;
    my ($targetname,$targetsize,$targetstart,$targetend) = @result[13..16];
    next if (($targetstart.$targetsize.$targetend) =~ /\D/);
#    my $offset = extract_offset($targetname);
    $targetstart -= $HITBUFFER_SIZE; 
    $targetstart = 0 if ($targetstart<0);
    $targetend += $HITBUFFER_SIZE;
    &save_target_region($targetname, $targetstart, $targetend, $targetsize);
    $hitcount = $.;
}
&cut_target_regions();

print COMMENT 
    "We have to check ".with_number($hitcount," BLAT hit")." in ".
    with_number(scalar keys %dbase," target sequence").".\n";
print COMMENT "(Rerun BLAT with other options to reduce / increase this number).\n\n";
print COMMENT "Reading sequences from \"$dbasefile\":\n";

seek BLOUT, 0,0;
$.=0;
&read_targets($dbasefile);
&validate_targets;

my %hits = map { $_ => [] } (keys %queries);

###################### postprocessing of BLAT hits ######################
my $dot_stands_for = max (1,$hitcount / 100);

print COMMENT "Processing BLAT hits:\n";
BLATLINES: while(<BLOUT>)
{
    next if (/^\#/ || /^$/);  # skip comments and empty lines
    print COMMENT "." unless ($. % $dot_stands_for);
    printf COMMENT "%2.0f%%", (($.) / $hitcount)*100  unless ($. % ($dot_stands_for * 10));
    printf COMMENT "\n" unless ($. % ($dot_stands_for * 50));

    # @result[..]=
    # 0:  number of matches
    # 1:  number of mismatches
    # 2:  number of repeated matches
    # 3:  number of unspecified nucleotides
    # 4:  number of gaps in query (not including zero size gaps)
    # 5:  total size of gaps in query
    # 6:  number of gaps in target
    # 7:  total size of gaps in target
    # 8:  strand
    # 9:  name of query
    # 10: total length of query
    # 11: start of aligned part of query 
    # 12: end of aligned part of query  
    # 13: name of target
    # 14: size of complete target
    # 15: start of aligned part of target
    # 16: end of aligned part of target
    # 17: number of blocks
    # 18: size of blocks in codons/aas, comma separated
    # 19: start positions of blocks in query
    # 20: start positions of blocks in target

    my @result = split "\t";
    next unless @result>=21;  # skip lines with less than 21 entries

    # read line of BLAT output file
    # the prefix 'q' is used for amino acid positions in the query
    my ($BLAT_matchcount, $BLAT_mismatchcount, $gapcount, $strand, $queryname, $querysize, 
	$targetname, $targetsize, $block_sizes, $qFroms, $tFroms) 
	= @result[0,1,5,8,9,10,
		  13,14,18..20];


    my ($undetermined, $insertions) = (0,0);
    &extend_to_full_name($targetname);
    next unless (exists  $dbase{$targetname});

    my $queryseq = $queries{$queryname};
    if (defined $queryseq && $querysize != length($queryseq))
    {
	if ($unique_queries) {
	    print STDERR "Warning: query length mismatch. This will produce unpredictable results!\n";
	} else {
	    $queryname = "${queryname}_[$querysize]";
	    $queryseq = $queries{$queryname};
	}
    }	
    unless (defined $queryseq)
    {
	die "No query sequence '$queryname' of length $querysize found.\nNonexistent queries in psl file";
    }
    $queryseq = uc $queryseq;

    $strand =~ s/^.//;
    my $complement = ($strand =~ /-/);
 
    chomp $tFroms;
    s/\,$// foreach ($block_sizes, $qFroms, $tFroms);
    my @block_sizes = split ",",$block_sizes;
    next if (($BLAT_mismatchcount >= ($MAX_MISMATCH + @block_sizes) && $MAX_MISMATCH) ||
	     ($BLAT_mismatchcount-1) / ($BLAT_matchcount + $BLAT_mismatchcount) > 1-$MIN_IDENTITY);
    my @qFroms = split ",", $qFroms;
    my @tFroms = split ",", $tFroms;
    my @qTos=(); 
    my @tTos=();
    foreach (0..$#block_sizes)
    {
	$tFroms[$_] += $complement ? -$targetsize : 0;
	push @qTos, $qFroms[$_]+($block_sizes[$_]);
	push @tTos, $tFroms[$_]+($block_sizes[$_]*3);
    }
    my ($querystart, $queryend) = ($qFroms[0], $qTos[-1]);
    my ($targetstart, $targetend) = ($tFroms[0], $tTos[-1]);
 

    # ignore if internal partial match is too small
    next BLATLINES if ($querystart > 5 && $queryend < $querysize-5 && $BLAT_matchcount < $INTERNAL_MIN_MATCH);

    ### Process Introns
    my @targetlocations = ($targetstart);
    my @querylocations = ($querystart*3);
    my @types = ("exon");

    foreach (my $inno=0; $inno<$#block_sizes; $inno++) {
	my ($tfrom, $tto, $qfrom, $qto) = ($tTos[$inno], $tFroms[$inno+1], 
					   $qTos[$inno], $qFroms[$inno+1]);

	my $qgaplen=($qto-$qfrom);
	my $offset = 0;
	my $type = "gap";
	my $additional = $tto-$tfrom-$qgaplen*3;

	# try to add exon
	if (2 < $qgaplen && $qgaplen <= $MAX_GAPLEN  && $additional>=-2) {
	    my $gapstr=substr($queryseq, $qfrom+1, $qgaplen-2);
	    my $pattern=&revtranslate_to_pattern($gapstr);
	    if ($pattern && 
		&subseq($targetname, $tfrom, $tto-$tfrom, $complement) =~ /$pattern/) {
		$gapcount -= ($qgaplen-2);
		splice(@tFroms,      $inno+1, 0, $tfrom + $-[0]);
		splice(@tTos,        $inno+1, 0, $tfrom + $+[0]);
		splice(@qFroms,      $inno+1, 0, $qfrom+1);
		splice(@qTos,        $inno+1, 0, $qto-1);
		splice(@block_sizes, $inno+1, 0, $qgaplen-2);
		redo;
	    }
	}

	if ($qgaplen <= $MAX_GAPLEN || $tfrom > $tto) {
	    ### DEBUG PART
	    print COMMENT "WARNING. Have query overlaps.\n" if ($qgaplen < 0);
            ### END OF DEBUG PART
	    if ($additional == 0) {
		$gapcount -= $qgaplen; next; # gap was closed completely: no location needed
	    } elsif ($additional < 0) {
		### case 1: no additional bases for an intron
		### here we just fill up optimally as much of target sequence
		### as we can; this is always the case when qgaplen > MAX_GAPLEN
		$type = "seqshift" if ($qgaplen <= $MAX_GAPLEN);

		# discard overlapping parts
		my $tframe = ($tto-$tfrom) % 3;
		if ($tfrom > $tto) { 
		    my $tgap = $tfrom - $tto + $tframe;
		    $tto += $tgap;
		    $tfrom -= $tgap;
		    $qto += $tgap/3;
		    $qfrom -= $tgap/3;
		    $qgaplen = $qto-$qfrom;
		    $gapcount += $tgap/3*2;
		}
		# find optimal location for inserting frameshift
		my $tseq = &subseq($targetname, $tfrom, $tto-$tfrom, $complement);
		my $closed_gaplen = $tto-$tfrom-$tframe;
		my ($max_mismatches, $best_gappos) = ($qgaplen, 0);
#		for (my $tgappos=0; $tgappos <= $closed_gaplen; $tgappos+=3) {
		for (my $tgappos=$closed_gaplen; $tgappos >= 0; $tgappos-=3) {
		    my $gap_trans = $ttable->translate(
			substr($tseq,0,$tgappos).substr($tseq,$tgappos+$tframe));
		    my $rightlen = ($closed_gaplen-$tgappos)/3;
		    my $querypart = substr($queryseq,$qfrom,$tgappos/3).substr($queryseq,$qto-$rightlen,$rightlen);
		    my $new_mismatches = &str_diff_count_mismatches($gap_trans, $querypart);
		    if ($new_mismatches < $max_mismatches) { 
			($max_mismatches, $best_gappos) = ($new_mismatches, $tgappos);
			last if ($new_mismatches == 0);
		    }
		}
		# fill up to optimal location
		$tfrom += $best_gappos;
		$qfrom += $best_gappos/3;
		$tto = $tfrom + $tframe;
		$qto -= ($closed_gaplen-$best_gappos)/3;
		$gapcount -= $closed_gaplen/3;
	    } else { 
                ### case 2: we have a gap of nucleotides 
                ### this is an intron or additional nucleotides (e.g. frameshift)

		# calculate marginal parts to be checked; always check matching regions
		my $max_overlap = int($qto - $querylocations[-1]/3)-1;
		my $left_overlap = max(0, $max_overlap - $qgaplen);
		if ($PRED_OVERLAP <= $left_overlap) { 
		    for (my $i = 1; $i < $max_overlap; $i++) {
			if ($ttable->translate(&subseq($targetname, $tto - 3*$i, 3, $complement)) 
			    ne substr($queryseq, $qto-$i,1)) {
			    $left_overlap = max($i-$qgaplen, $PRED_OVERLAP);
			    last;
			}
		    }
		}
		my $right_overlap = $block_sizes[$inno+1]-1;
		$max_overlap = $right_overlap + $qgaplen;
		if ($PRED_OVERLAP <= $right_overlap) { 
		    for (my $i = 0; $i < $max_overlap-1; $i++) {
			if ($ttable->translate(&subseq($targetname, $tfrom + 3*$i, 3, $complement)) 
			    ne substr($queryseq, $qfrom+$i,1)) {
			    $right_overlap = max($i-$qgaplen+1, $PRED_OVERLAP);
			    last;
			}
		    }
		}
		my $gapstart=($qfrom-$left_overlap);
		my $qgaplenplus=$qgaplen + $left_overlap + $right_overlap;
		### DEBUG PART
		if ($qgaplenplus < 1) {
		    print COMMENT "ERROR. There should not be query overlaps. Skipped this hit.\n";
		    next BLATLINES;
		}
		### END OF DEBUG PART

		# find optimal location for inserting intron
		my ($dss_seq, $ass_seq, $gap_aa) = 
		    (&subseq($targetname, $tfrom-3*$left_overlap, $qgaplenplus*3+2, $complement),
		     &subseq($targetname, $tto-($qgaplenplus-$right_overlap)*3-2, 
			     $qgaplenplus*3+2, $complement),
		     substr($queryseq."*", $gapstart, $qgaplenplus));
		### DEBUG PART
		if ($gap_aa eq "*") {
		    &run_bug("Internal error. This shouldn't happen!");
		}
		### END OF DEBUG PART
		my ($get_offset, $mismatchadd, $gap_trans, $splicesite_score) = 
		    &get_intron_offset($gap_aa, $dss_seq, $ass_seq,
				       $additional <= $MAX_SHIFTLEN,
				       3*($qgaplenplus-$right_overlap)   );
		$get_offset -= 3*$left_overlap;
		$mismatchadd=$qgaplen if ($qgaplen < $mismatchadd);
		
		# close the gap in query, if possible with few mismatches; otherwise keep the gap
		if ($mismatchadd <= 2 && $tfrom+$get_offset >= $targetlocations[-1]+2)
		{   
		    $qto=$qfrom;
		    $tto -= $qgaplen*3;
		    $gapcount -= $qgaplen; 
		    $offset = $get_offset;
		    $type = $additional <= $MAX_SHIFTLEN ? "seqshift" : $splicesite_score>0.5 ? "intron" : "intron?";
		}
	    }
	}
	push @types, ($type, "exon"); 
	push @targetlocations, ($tfrom+$offset, $tto+$offset);
	push @querylocations, ($qfrom*3+$offset, $qto*3+$offset);
    } # end process introns


    # try to add suffix that wasn't found by BLAT
    if ($querysize-$MAX_GAPLEN <= $queryend && $queryend < $querysize && $ADD_MARGINS)
    {
	my $dss =  &subseq($targetname, $targetend-3, 7, $complement);
	my $suffix = substr($queryseq, $queryend-1)."*";
	my $pattern = &get_splicesite_patterns_forward($suffix,$dss);
	my $target = &subseq($targetname, $targetend, $HITBUFFER_SIZE-$REGION_SIZE, $complement);
	if ($pattern && $target =~ /$pattern/)
	{
	    my $qnto = $queryend*3;
	    my ($tfrom, $tto) = ($-[0]+2, $+[0]-3);
	    $_ += $targetend foreach ($tfrom, $tto);
	    my $offset = $querysize*3 - $qnto - ($tto-$tfrom);
	    push @querylocations, ($qnto + $offset) x 2;
	    push @targetlocations, ($targetend + $offset, $tfrom);
	    push @types, ("intron", "exon");
	    $targetend = $tto;
	    $queryend = $querysize;
	}
    }

    # try to add prefix that wasn't found by BLAT
    if (0 < $querystart && $querystart <= $MAX_GAPLEN && $ADD_MARGINS)
    {
	my $ass = &subseq($targetname, $targetstart-4, 7, $complement);
	my $prefix = substr($queryseq, 0, $querystart+1);
	my $searchstart = $targetstart - $HITBUFFER_SIZE + $REGION_SIZE;
	$searchstart=0 if (!$complement && $searchstart < 0);
	my $target = &subseq($targetname, $searchstart, $targetstart-$searchstart-1, $complement);
	my $pattern = &get_splicesite_patterns_backward($prefix, $ass);
	if ($pattern && $target =~ /^.*($pattern)/)
	{
	    my $qfromref = \$querylocations[0];
	    my ($tfrom, $tto) = ($-[1], $+[1]-2);
	    $_ += $searchstart foreach ($tfrom, $tto);
	    my $offset = $$qfromref - ($tto-$tfrom);
	    $targetlocations[0] -= $offset;
	    $$qfromref = $tto-$tfrom; 
	    unshift @querylocations, (0, $$qfromref);
	    unshift @targetlocations, ($tfrom, $tto);
	    unshift @types, ("exon", "intron");
	    $targetstart = $tfrom;
	    $querystart = 0;
	}
    }
    push @targetlocations, $targetend;
    push @querylocations, $queryend*3;

    # if last match reaches end of query, check if there is a stop codon
    my $stopcodon;
    if ($queryend == $querysize)
    {
	$stopcodon = &subseq($targetname, $targetend, 3, $complement);
	$stopcodon = undef unless($ttable->is_ter_codon($stopcodon));
    } 

    # matchsize is the length of the aligned part of the query
    # matchsize differs from querysize only when partial hits are found
    my $matchsize = ($queryend-$querystart);
	
    my @diff_strings = 
	&calculate_diff_strings($targetname, $queryname, \@targetlocations, \@querylocations, $complement);
    
    ### debug part
    if (length(join "", map { /:(.*)$/; $1 } @diff_strings) + $gapcount != $matchsize)
    {
	&run_bug("Incorrect calculation of unmatched aa's in line $.!\n");
    }
    ### end of debug part

    my @matchings = ();
    my @mismatchpos = ();
    my $mismatchcount = 0;

    ### Construct matchings hash
    foreach (@types)
    {
	my $previous = $matchings[-1];
	my ($tfrom, $qnfrom) = (shift @targetlocations, shift @querylocations);
	my ($tto, $qnto) =  ($targetlocations[0], $querylocations[0]);
	my ($qfrom, $qto) = map ( &nu_to_aa($_), ($qnfrom, $qnto) );
	my $current = { "${TPFX}_start" => $tfrom,
			"${TPFX}_end" => $tto,
			"${QNPFX}_start" =>  $qnfrom,
			"${QPFX}_start" => $qfrom,
			"${QNPFX}_end" => $qnto,
			"${QPFX}_end" => $qto };
	if (/seqshift/) {
	    my $translation = $ttable->translate(&subseq($targetname, $tfrom, $tto-$tfrom, $complement)."--");
	    if ($SPLIT_ON_FS) {
		$current->{translation} = $translation;
	    } else {
		my $translen = length($translation);
		if ($qnto == $qnfrom) {
		    $insertions += $translen;
		} else {
		    $gapcount -= $translen;
		    $undetermined += $translen;
		    push @{$previous->{undeterminedlist}}, ($qfrom+1 .. $qfrom + $translen);
		    push @{$previous->{gaplist}}, ($qfrom + $translen +1 .. $qto);
		}
		$previous->{translation} .= $translation;
		push @{$previous->{seqshifts}}, $current;
		while ($translation =~ /[*]/g) {
		    push @{$previous->{inframe_stopcodons}}, $qfrom;
		}
		next;
	    }
	}
	if (/exon/) 
	{
	    my $diffstr;
	    (shift @diff_strings) =~ /^(.*):(.*)$/;
	    ($current->{translation}, $diffstr) = ($1, $2);
	    $current->{$_} = [] 
		foreach("seqshifts","mismatchlist","gaplist",
			"undeterminedlist","inframe_stopcodons");
	    while ($diffstr =~ /([X ])/g) {
		my $qpos =  $-[0] + $qfrom + 1;
		if ($1 eq 'X') {
		    push @{$current->{mismatchlist}}, $qpos;
		    ### DEBUG PART ###
		    push @mismatchpos, $qpos;
		    ### END OF DEBUG PART ###
		    $mismatchcount++;
		} elsif ($1 eq ' ') {
		    push @{$current->{undeterminedlist}}, $qpos;
		    $undetermined++;
		}
	    }
	    while ($current->{translation} =~ /[*]/g) {
		push @{$current->{inframe_stopcodons}}, $-[0] + $qfrom + 1;
	    }
	    if ($previous && $previous->{type} eq "exon")
	    {   # add to previous exon...
		@$previous{"${TPFX}_end","${QNPFX}_end","${QPFX}_end"} = ($tto, $qnto, $qto);
		push (@{$previous->{$_}}, @{$current->{$_}} )
		    foreach ("mismatchlist","undeterminedlist","inframe_stopcodons");
		$previous->{translation} .= $current->{translation};
		# and don't put on list
		next;
	    }
	} elsif (/intron/) {
	    @$current{"${QNPFX}_pos", "${QPFX}_pos"} = ($qnfrom, $qfrom) ;
	}
	$current->{type} = $_;
	push @matchings, $current;
    } # end foreach (@types)

    ### debug part
    my @mismatchcounts = map { scalar @{$_->{mismatchlist}} } grep { $_->{type} eq "exon" }  @matchings;
    if ($mismatchcount != sum (@mismatchcounts,0))
    {
	print STDERR "ID: $.\nMismatchcount: $mismatchcount\n".
	    scalar(@mismatchpos),"\naufgeschlüsselt: (".join("+",@mismatchcounts).")\n";
	print STDERR "Exons: ".join("/", @targetlocations)."\n";
	&run_bug("Incorrect calculation of mismatches!");
    }
    &run_bug("Schlimm!") unless defined ($qTos[-1]);
    ### end of debug part

    # matchcount is the number of matched positions of the query sequence
    my $matchcount = $matchsize - $mismatchcount - $gapcount - $undetermined; 

    # the score is the fraction of matches minus mismatches; 
    # any query should have at least one hit with a score of $MIN_BEST_SCORE (=0.3)
    my $score = ($matchcount-$mismatchcount)/$querysize;
    
    ### Collect all data for this BLAT hit
    my $hitref = { "${QPFX}_len" => $querysize, "${QPFX}_start" =>  $querystart, "${QPFX}_end" => $queryend,  
		   "target" => $targetname, "target_len" => $dbase{$targetname}{length},
		   "${TPFX}_start" => $targetstart, "${TPFX}_end" => $targetend, 
		   "strand" => $complement?"-":"+", "complement" => $complement, 
		   "stopcodon" => $stopcodon,
		   "score" => $score, "unmatched" => $gapcount,
		   "undetermined" => $undetermined, "additional" => $insertions,
		   "matches" => $matchcount, "mismatches" => $mismatchcount, "mismatchcounts" => \@mismatchcounts,
		   "ID" => $., "matchings" => \@matchings };
    next if (($mismatchcount >= $MAX_MISMATCH && $MAX_MISMATCH) || 
	     $matchcount/($mismatchcount+$matchcount) < $MIN_IDENTITY);
#            $matchcount/$matchsize
    push @{$hits{$queryname}}, $hitref;
    $dbase{$targetname}{users}{$.} = 1;
} # end BLATLINES
print COMMENT "done\n";
print COMMENT sum (map {scalar @{$hits{$_}}} keys %hits)." of $hitcount hits saved.\n";

# free memory of unused sequences
my $delcount = 0;
foreach (keys %dbase) {
    unless (keys %{$dbase{$_}{users}}) { 
#	print COMMENT "Deleting sequence $_ - not needed anymore.\n";
	$delcount++;
	delete $dbase{$_};
    } 
} 


print COMMENT "Discarded $delcount unused sequences.\n" if ($delcount);
print COMMENT "Now assembling BLAT hits...\n";
#print COMMENT "[";
$delcount = 0;

### Sort BLAT hits
# before postprocessing is completed, all results are stored in %hits which is a hash 
# indexed by query names whose values are references to lists, each list entry
# referring to the information about one hit for the query

# after postprocessing, only one hitref (or sequence of consecutive hits) will remain
# for each query, and sequences are added

my @missing = ();
my @hitkeyqueries = keys %hits;

foreach my $queryname (@hitkeyqueries) {
    my @usercounts = map { scalar keys %{$_->{users}} } values %dbase;
#    print COMMENT "$queryname: (".join(",",@usercounts).") => ".sum(@usercounts);
    my $hitlist = $hits{$queryname};
    my $is_copy = $queryname =~ s/_\(\d+\)$//;

    if (@$hitlist) {
	@$hitlist = sort { $b->{"score"} <=> $a->{"score"} } @$hitlist;
	if (grep { ! defined $_->{"${QPFX}_len"} } @$hitlist) {
	    print STDERR "${QPFX}_len undefined!\n";
	}
#	unless (grep { ($_->{"${QPFX}_end"} - $_->{"${QPFX}_start"}) / $_->{"${QPFX}_len"} >= $MIN_BEST_SCORE } @$hitlist) {
	if ($hitlist->[0]{score} < $MIN_BEST_SCORE) {
 	    @$hitlist = ();
	    printf COMMENT
 		"Dropping $queryname because largest piece has score less than %f\n", $MIN_BEST_SCORE;
	}
    }
    unless (@$hitlist)
    {
	push @missing, $queryname; 
	delete $hits{$queryname};
	next;
    }

    my $i=1;
    $i++ while ($i < @$hitlist && $hitlist->[$i]{score} >= $SHOWALL_SCORE);
    if ($i>1) { # we found an additional high scoring hit
	my @alt_names = map("${queryname}_($_)", (1..$i-1));
	@hits{@alt_names} = map [ $_ ], splice(@$hitlist, 1, $i-1);
	push @hitkeyqueries, @alt_names; # for blessing 
	print COMMENT 
	    "Added ".with_number($i-1," additional high scoring hit")." as queries ${queryname}_(1)".
	    ($i==2 ? "" : $i==3 ? ",_(2)" : ",...,_(".($i-1).")").".\n";
	$i=1;
    }
  CAND: 
    while ($i<@$hitlist) {
	foreach (@$hitlist[0..($i-1)])	{
	    if (&incompatible($_, $hitlist->[$i])) {
		remove_user(@{$hitlist->[$i]}{"target","ID"});
		splice(@$hitlist, $i, 1);
		next CAND;
	    }
	}
	$i++;
    } # end CAND

    @$hitlist = sort { $a->{"${QPFX}_start"} <=> $b->{"${QPFX}_start"} } @$hitlist;
    my $queryseq = $queries{$queryname};
    

    # adjusting with neighbouring hits
    foreach (0..$#$hitlist) {
	my $current = $hitlist->[$_];
	my ($targetname, $complement, $tfrom, $tto, $tlen, $qto, $qlen) 
	    = @{$current}{"target", "complement", "${TPFX}_start", "${TPFX}_end", 
			  "target_len", "${QPFX}_end", "${QPFX}_len"};
	if ($_==0) {
	    my $upstream_size = min ($REGION_SIZE, $tfrom + ($complement ? $tlen : 0));
	    $current->{upstream} = &subseq($targetname, $tfrom - $upstream_size, $upstream_size, $complement);
	    my $qfrom = $current->{"${QPFX}_start"};
	    push @missing, "${queryname}[1..$qfrom]" if ($qfrom > 0 && !$is_copy);
	} 
	my $prevnucs = ($complement ? 0 : $tlen) - $tto;
	
	if ($_< $#$hitlist) {
	    my $nexthit = $hitlist->[$_+1];
	    my ($nexttfrom, $nextqfrom, $nexttlen) = @{$nexthit}{"${TPFX}_start", "${QPFX}_start", "target_len"};

	    my $t_o = $current->{target_overlaps};                   # copy to variable to prevent autovivification of key
	    my $has_overlap = grep { $_ == $nexthit->{ID} } @$t_o;   # true if ID of next hit is among overlaps

	    my $postnucs = ($nexthit->{complement}? $nexttlen : 0) + $nexttfrom;
	    my $qgaplen = $nextqfrom-$qto;
	    my $additional = $prevnucs+$postnucs-$qgaplen*3;
#	    if ($qgaplen <= $MAX_GAPLEN || $has_overlap) {
		if ($additional <= 0 || $has_overlap) { # not an intron, so here just extend the exons of both parts
		    # all dna sequences here have length divisible by 3, so dna_add is not needed
		    my ($prevrem, $postrem) = 
			(&subseq($targetname, $tto, $prevnucs, $complement),
			 &subseq($nexthit->{target}, $nexttfrom-$postnucs, $postnucs, $nexthit->{complement}));
		    my $dna_add = "n" x (-$postnucs%3);
		    if ($has_overlap) {
			$current->{matchings}[-1]{overlap} = $additional;
			$additional-=$prevnucs;
			&replace_in_last($current, "", "", $postrem, $additional);
			$qto-=&nu_to_aa($additional);
			$prevnucs = 0;
			if ($postnucs % 3 == 2) { 
			    $dna_add = &subseq($targetname, $tto-$additional-1,1,$complement);
			}
		    } elsif ($prevnucs <= 4) {
			&replace_in_last($current, $prevrem, substr($queryseq, $qto, &nu_to_aa($prevnucs)),
					 ($postnucs <= 4 && $additional==0) ? $postrem : "nn");
			$qto += &nu_to_aa($prevnucs);
			$prevnucs = 0;
			if ($additional == 0) { $dna_add = $prevrem; }
		    } 
		    if ($postnucs <= 4 || $has_overlap) {
			my $postaas = &nu_to_aa($postnucs);
			&replace_in_next($nexthit, $postrem, substr($queryseq, $nextqfrom-$postaas, $postaas),
					 $dna_add);
			$nextqfrom -= $postaas;
			$postnucs = 0;
		    }
		} elsif ($qgaplen <= $MAX_GAPLEN) { 
		    my $left_border = $current->{matchings}[-1]{"${QPFX}_start"} + $MIN_REMAIN;
		    my $right_border = $nexthit->{matchings}[0]{"${QPFX}_end"} - $MIN_REMAIN;
		    my $left_overlap = min ( max ( $PRED_OVERLAP, $PRED_OVERLAP-$qgaplen ), max ($qto - $left_border, 0) );
		    my $right_overlap = min ( max ( $PRED_OVERLAP, $PRED_OVERLAP-$qgaplen ), max ($right_border - $nextqfrom, 0));
		    my $qgaplenplus = $qgaplen + $left_overlap + $right_overlap;
		    if ($left_border > $qto || $right_border < $nextqfrom || $qgaplenplus < 0) {
			&run_bug("This is a bug and shouldn't happen. Please report.");
		    }
		    my $asslen = min ( $postnucs+3*$right_overlap, 3*$qgaplenplus +2);
		    my ($dss_seq, $ass_seq, $gap_aa) =
			(&subseq($targetname, $tto-3*$left_overlap, 
				 min ($prevnucs+3*$left_overlap, 3*$qgaplenplus +2),
				 $complement),
			 &subseq($nexthit->{target}, $nexttfrom+3*$right_overlap-$asslen, $asslen, $nexthit->{complement}),
			 substr($queryseq, $qto-$left_overlap, $qgaplenplus));
		    my $preferred_offset = 3*int($qgaplenplus/2);
		    my ($get_offset, $mismatchadd, $gap_trans, $splicesite_score) =
			#always look for splice sites
			&get_intron_offset($gap_aa, $dss_seq."xx", "xx".$ass_seq, 0 , $preferred_offset); # todo: $additional <= $MAX_SHIFTLEN
		    if (($mismatchadd <= 2 && $splicesite_score > 0.5) || $qgaplen<=0 ) { # only add sure introns in the end
			my $asslen = 3*$qgaplenplus - $get_offset;
			my ($dss, $ass) = (substr($dss_seq."xx",0,$get_offset+2),substr("xx".$ass_seq,-$asslen-2));
			my $dss_inpart = substr($dss,-2,2,"");
			my $ass_inpart = substr($ass,0,2,"");
			
			# for replace_in_last we just set $dna_add to length 2
			&replace_in_last($current, $dss,
					 substr($gap_aa, 0, &nu_to_aa($get_offset)),
					 substr($ass,0,2), 3*$left_overlap );
#			$prevnucs += ($get_offset - 3*$left_overlap);
			# for replace_in_next we set $dna_add such that its length equals the reading frame
			&replace_in_next($nexthit, $ass,
					 substr($gap_aa, &nu_to_aa($get_offset)),
					 substr($dss,-($get_offset%3), $get_offset%3),
 					 3*$right_overlap);
#			$postnucs += (3*$qgaplen+ 3*$left_overlap-$get_offset);
			$qto += &nu_to_aa($get_offset) - $left_overlap;
			$nextqfrom -= ($qgaplenplus - $right_overlap - &nu_to_aa($get_offset));
			## DEBUG PART
			if ($qto != $nextqfrom) {
			    &run_bug("Internal error. Serious bug with intron offset calculation. Please report!");
			}
			## END OF DEBUG PART
			my ($dsstype, $asstype) = ("intron", "intron");
			if ($splicesite_score <= 0.5) {
			    $dsstype .= "?" unless ($dss_inpart =~ /^g[ct]$/);
			    $asstype .= "?" unless ($ass_inpart eq "ag");
			}
			&add_final($current,$dsstype); $prevnucs=0;
			&add_initial($nexthit,$asstype); $postnucs=0;
		    } # mismatchadd <= 2
		    ## DEBUG PART
		    elsif ($qto == $nextqfrom && $mismatchadd > 2) {
			print STDERR "current: ".$current->{ID}." next hit: ".$nexthit->{ID}."\n";
			print STDERR "prevnucs=$prevnucs postnucs=$postnucs mismatchadd=$mismatchadd\n";
 			substr($dss_seq,$left_overlap*3,0)=".";
 			substr($ass_seq,-$right_overlap*3,0) = "." if ($right_overlap);
			print STDERR "to=nextfrom=$qto additional=$additional length=".length($queryseq)."\n";
			print STDERR "dss=$dss_seq ass=$ass_seq match=".substr($queryseq, $qto-$left_overlap, $left_overlap)."()".substr($queryseq,$nextqfrom,$right_overlap)."\n";
			print STDERR "offset=$get_offset\n";
			&run_bug("Internal error. Adding no aas should not add errors. Please report!");
		    }
		    ## END OF DEBUG PART
		}  # additional > 0
#	    } # gaplen <= MAX_GAPLEN
	    if ($qto < $nextqfrom) {
		# add gap sequences in all cases in which we have a gap; 
		# empty gap sequences are needed later to determine that there is a gap to neighbouring hits
		$current->{downstream_gap} = &subseq($targetname, $tto, min($REGION_SIZE, $prevnucs), $complement);
		my $next_upstream_size = min ($REGION_SIZE, $postnucs);
		$nexthit->{upstream_gap} = &subseq($nexthit->{target}, $nexttfrom-$next_upstream_size, $next_upstream_size, $nexthit->{complement});
		push @missing, "${queryname}[".($qto+1)."..$nextqfrom]";
	    } elsif ($prevnucs > 0 || $postnucs > 0)  {
		### DEBUG PART
		&run_bug ("Internal error: gap is closed - we should not have regions to show here. Please report!");
		### END OF DEBUG PART
	    }
	} else { # $_==$#$hitlist  or no assembling
	    $current->{downstream} = &subseq($targetname, $tto, min($REGION_SIZE, $prevnucs), $complement);
	    push @missing, "${queryname}[".($qto+1)."..$qlen]" if ($qto < $qlen && !$is_copy)
	}

	foreach (@{$current->{matchings}}) {
	    my ($reftfrom, $reftto, $refqnfrom, $refqnto) = \@$_{"${TPFX}_start","${TPFX}_end", "${QNPFX}_start", "${QNPFX}_end"};
	    my ($qfrom, $qto) = @$_{"${QPFX}_start","${QPFX}_end"} = map ( &nu_to_aa($_), ($$refqnfrom, $$refqnto) );
	    $_->{seq} = &subseq($targetname, $$reftfrom, $$reftto-$$reftfrom, $complement);
	    if ($_->{type} eq "intron") {
		if ($_->{seq} =~ s/\.{3}(\.+)/.../) { $$reftto -= length($1); }
		elsif ($_->{seq} =~ s/^(\.+)\.{3}/.../) { $$reftfrom += length($1); }
	    } elsif ($_->{type} eq "exon" || $_->{type} eq "gap" ) {
		$_->{"${QPFX}_seq"} = substr($queryseq, $qfrom, $qto-$qfrom);
	    }
	}
	&remove_user($targetname,$current->{ID});
	my $qfrom = $current->{ "${QPFX}_start" };
	## DEBUG PART
	if ($qto != $current->{"${QPFX}_end"}) {
	    print STDERR "";
	    &run_bug("BUG: coordinates mismatch. id:".
		$current->{ID}." qto=$qto prot_end=".$current->{"${QPFX}_end"}.
		". Please report!");
	}
	## END OF DEBUG PART
	$current->{"${QPFX}_seq"} = substr($queryseq, $qfrom, $qto-$qfrom);
	$current->{score} = sprintf("%5.3f",($current->{matches} - $current->{mismatches})/$qlen);
	&set_problem_field($current);

	# choose keys for YAML output
	my @bless_list = split ",",$HIT_BLESS_LIST;
	my %default_values = ( 
	    matches => $current->{"${QPFX}_len"}, 
	    mismatches => 0,
	    unmatched => 0,
	    undetermined => 0,
	    additional => 0,
	    upstream_gap => "",
	    downstream_gap => "",
	    "${QPFX}_start" => 0,
	    "${QPFX}_end" => $current->{"${QPFX}_len"} );
	if ($HIDE_UNDEF) {
	    @bless_list = grep { defined $current->{$_} } @bless_list;
	}
	if ($HIDE_DEFAULTS) {
	    foreach (@bless_list) {
		print STDERR "\$current->$_\n" unless defined $current->{$_};
		print STDERR "\$default_values{$_}\n" if (exists $default_values{$_} && !defined $default_values{$_});
	    }
	    @bless_list = grep { ! (exists $default_values{$_} && $current->{$_} eq $default_values{$_}) } @bless_list;
	}
	YAML::Bless($current)->keys(\@bless_list);
	foreach my $match (@{$current->{matchings}}) {
	    my $type = $match->{type};
	    $type = "intron" if ($type =~ /intron/);
	    $type = "any" unless (exists $MATCH_BLESS_LISTS{$type});
	    my @bless_list = split ",",$MATCH_BLESS_LISTS{$type};
	    if ($type eq "exon") {
		foreach my $fs (@{$match->{seqshifts}}) {
		    YAML::Bless($fs)->keys( [ grep { exists $fs->{$_} } @bless_list ] );
		}
	    }
	    if ($HIDE_UNDEF) {
		foreach (values %$match) {
		    $_ = undef if ( ref($_) eq 'ARRAY' && !@$_);
		}
		@bless_list = grep { defined $match->{$_} } @bless_list;
	    }
	    YAML::Bless($match)->keys(\@bless_list);
	}
    } # end foreach (@hitlist)
} # end foreach (@hitkeyqueries)


###################### output of results ################################

print "### Scipio v1.0 output\n"; 
print "# query file    $queryfile\n";
print "# target file   $dbasefile\n";
print "# BLAT output   $BLAT_OUTPUT\n\n";

my @querybless = grep { defined $hits{$_} } sort keys %hits;
YAML::Bless(\%hits)->keys(\@querybless);
print YAML::Dump(\%hits);


if (@missing) {
     print "### not found: ".(join ",",@missing)."\n\n";
}
print "### end of Scipio output\n";
exit scalar @missing;


####################### function definitions ####################################

sub run_bug {
    my $errmess = shift;
    if ($DONT_DIE) { print STDERR "$errmess\n" }
    else { die $errmess }
}

sub revcom {
    my $result = reverse shift;
    $result =~ tr/tcyawmhgksbrdvnTCYAWMHGKSBRDVN/agrtwkdcmsvyhbnAGRTWKDCMSVYHBN/;
    return $result;
}

sub with_number {
    my ($n, $s, $alt) = @_;

    $alt = "${s}s" unless (defined $alt);
    return ($n==1) ? "$n$s" : "$n$alt";
}

sub diff_str {
    my ($s1, $s2, $count, $s1from, $s2from) = @_;
    $s1from = 0 unless (defined $s1from);
    $count = length($s1)-$s1from unless (defined $count);

    $s2from = $s1from unless (defined $s2from);
    if (length($s1)<$s1from+$count || length($s2)<$s2from+$count)  {
	print STDERR "diff_str returns undef: input is\n" .
	    "s1: $s1\ns2: $s2\ncount: $count\nfrom1 :$s1from\nfrom2: $s2from\n";
    }
    return undef if (length($s1)<$s1from+$count || length($s2)<$s2from+$count);
    my @s = split //, substr($s1, $s1from, $count);
    my @t = split //, substr($s2, $s2from, $count);
    
    my $result=""; 
    foreach (@s)  {  # $_: translation <-> $comp: query
	my $comp = shift(@t);
	### DEBUG part
	&run_bug ("'-' in query or translation!!!") if ($comp eq '-' || $_ eq '-');
	### end of DEBUG part
	# count 'X' as undetermined
	# others as match/mismatch
	$result .= ($comp =~ /^[X-]$/ || /^[X-]$/) ?  " " : ( $comp eq $_ ? "|" : "X" );
	
    }
    return $result;
}

sub str_diff_list {
    my @result=(); my $pattern=shift;
    local $_ = &diff_str(@_);
    push @result, $-[0] while(/$pattern/g);
    return @result;
}

sub str_diff_count_all {
    local $_ = &diff_str(@_);
    s/\|//g;
    return length;
}

sub str_diff_count_undeterm {
    local $_ = &diff_str(@_);
    s/[^ ]//g;
    return length;
}

sub str_diff_count_mismatches {
    local $_ = &diff_str(@_);
    s/[^X]//g;
    return length;
}

sub str_pos { 
    my @result=(); my $pattern=shift;
    local $_ = shift;
    push @result, $-[0] while(/$pattern/g);
    return @result;
}

sub extend_to_full_name {
    my $result = $full_name{$_[0]};
    $_[0] = $result if (defined $result);
}

sub save_target_region {
    my ($targetname, $from, $to, $total_length) = @_;
    $from = 0 if ($from<0);
    return if $to<=$from;

    unless (exists $dbase{$targetname})  {
	$dbase{$targetname} = { "length" => $total_length, "users" => {} };
    }
    my $target = $dbase{$targetname};
    $target->{length} = $total_length if ($target->{length} < $total_length);

    foreach (keys %$target)  {
	my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	next unless (defined $pieceto);
	return if ($piecefrom <= $from && $to <= $pieceto);
	if ($from <= $piecefrom && $pieceto <= $to) {
	    delete $target->{$_};
	    next;
	}
	if ($piecefrom < $from && $from <= $pieceto) {
	    delete $target->{$_};
	    $from = $piecefrom;
	}
	elsif ($piecefrom <= $to && $to < $pieceto) {
	    delete $target->{$_};
	    $to = $pieceto;
	}
    }
    $target->{"$from,$to"} = "";
}

sub cut_target_regions {
    foreach my $target (values %dbase) {
	my $length = $target->{length};
	foreach (keys %$target) {
	    my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	    next unless (defined $pieceto && $pieceto>$length);
	    my $piece = $target->{$_};
	    delete $target->{$_};
	    $target->{"$piecefrom,$length"}="" if ($piecefrom < $length);
	}
    }
}
	
sub fill_target_region {
    my ($targetname, $sref, $offset) = @_;
    my $target = $dbase{$targetname};
    my $end = $offset+length($$sref);
    
    foreach (keys %$target) {
	my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	next unless (defined $pieceto);
	next if ($offset >= $pieceto || $piecefrom >= $end);
	if ($target->{$_}) { # we already have this target; something is wrong
	    die "Error: multiple target sequences with the same name '$targetname' found! Aborting...\n"
	}
	if ($piecefrom < $offset) {
	    delete $target->{$_};
	    $target->{"$piecefrom,$offset"} = "";
	    $piecefrom=$offset;
	}
	if ($pieceto > $end) {
	    delete $target->{$_};
	    $target->{"$end,$pieceto"} = "";
	    $pieceto=$end;
	}
	$target->{"$piecefrom,$pieceto"}=substr($$sref,$piecefrom-$offset,$pieceto-$piecefrom);
    }
}	

sub read_targets {
    my ($targetname, $offset);
    my $seq = "";
    open FASTAFILE, shift;
    my $cutoff_warning_no=0;
    my $seqcount=0;

    while (<FASTAFILE>) {
	if (s/^>\s*//) { # new target sequence
	    if ($seq) {  # save old sequence
		&fill_target_region($targetname, \$seq, $offset);
		$seq = "";
	    }
	    chomp;
	    if (/^(\S*)(\s+.*)$/) { # new target name contains white spaces
		$cutoff_warning_no++;
		$full_name{$1}=$_;
		$dbase{$_}=$dbase{$1};
		delete $dbase{$1};
	    }
	    $offset = 0;
	    $targetname = $_;
	    if (defined $dbase{$targetname} && $seqcount <= 80) {
		print COMMENT ($seqcount == 80) ? "[...] " : ".";
		$seqcount++;
	    }
	} elsif (defined $dbase{$targetname}) { # ignore sequences not among BLAT results
	    chomp;
	    $seq.=lc;
	    ### don't let $seq become larger than 20MBases
  	    if (length($seq) > 20000000) {
  		&fill_target_region($targetname, \$seq, $offset);
  		$offset += length($seq);
  		$seq = "";
 	    }
	}
    }
    if ($seq) {
	&fill_target_region($targetname, \$seq, $offset);
    }
    print COMMENT "done\n";	
    if ($cutoff_warning_no) {
	print STDERR "Warning: Some target sequences stripped at white spaces.\n";
    }
}

sub validate_targets {
    while (my ($targetname, $target) = each %dbase) {
	foreach (keys %$target) {
	    if (/^(\d+),(\d+)$/ && ($2 > $1) && $target->{$_} eq "")
	    {
		die "Missing or incomplete target sequence '$targetname'.\n".
		    "Nonexistent targets in psl file" unless ($ALLOW_MISSING_TARGETS);
		delete $dbase{$targetname};
		last;
	    }
	}
    }
}

sub nu_to_aa {
    # if changing the assignment of nucleotides, check in particular: postnucs(l. 863)
    my $arg = shift;
    my $residue = ($arg+1) % 3 -1;   # -1, 0, 1
    return ($arg-$residue)/3;
}

sub revtranslate_to_pattern {
    my $seq = uc shift;
    return $seq =~ /[^ABCDEFGHIKLMNPQRSTVWYZ*]/? undef :  
	join ("", map { "(".join("|",$ttable->revtranslate($_)).")" } split(//, $seq));
}

sub split_codons {
    my ($aa, $dna, $backward) = @_;
    my @result=();
    my @codons = $ttable->revtranslate($aa);
    foreach (@codons) {
	/^(.)(.)(.)$/;
	my ($c1,$c2,$c3) = ($1,$2,$3);
	if ($backward) {
	    if ($dna =~ /ag$c3$/) { 
		push @result, "$c1${c2}g[ct]"; 
	    } else {
		if ($dna =~ /ag$c2$c3$/) { push @result, "${c1}g[ct]"; }
		if ($dna =~ /ag$/) { push @result, "${_}g[ct]"; }
	    }
	} else {
	    if ($dna =~ /^$1g[ct]/) {
		push @result, "ag$c2$c3"; 
	    } else {
		if ($dna =~ /^$c1${c2}g[ct]/) { push @result, "ag$c3"; }
		if ($dna =~ /^g[ct]/) { push @result, "ag$_"; }
	    }
	}
    }
    return @result;
}

sub get_splicesite_patterns_forward {
    my ($suffix, $dna) = @_;
    $suffix =~ s/^(.)(.)//;
    my $longpatt = &revtranslate_to_pattern($suffix);
    my @result = &split_codons($1, $dna, 0);
    @result = ( "(".join("|", @result).")".&revtranslate_to_pattern($2) ) if @result;
    push @result, &split_codons($2, substr($dna,3),0);
    ### debug part
    print STDERR "Invalid suffix: $suffix" unless (defined $longpatt);
    ### end of debug part
    return undef unless (@result && defined $longpatt);
    return "(".join("|", @result).")".$longpatt;
}
   
sub get_splicesite_patterns_backward {
    my ($prefix, $dna) = @_;
    $prefix =~ s/(.)(.)$//;
    my $longpatt = &revtranslate_to_pattern($prefix);
    my @result = &split_codons($2, $dna, 1);
    @result = ( &revtranslate_to_pattern($1)."(".join("|",@result).")" ) if @result;
    substr($dna,-3)="";
    push @result, &split_codons($1, $dna, 1);
    return undef unless (@result && defined $longpatt);
    return &revtranslate_to_pattern($prefix)."(".join("|", @result).")";
}

sub subseq {
    my ($targetname, $from, $count, $complement) = @_;
 
    return undef if ($count<0);
    return "" if ($count==0);
    
    my $target = $dbase{$targetname};
    return undef unless (defined $target);
    $from=-$from-$count if ($complement);
    
    my $prelap = max ( 0, -$from);
    my $postlap = max ( 0, $from+$count-$target->{length} );
    if ($count < $prelap+$postlap) {
	return " " x $count;
    }
    my $result = (" " x $prelap).("." x ($count-$prelap-$postlap)).(" " x $postlap);
    foreach (keys %$target) {
	next unless (/^(\d+),(\d+)$/);
	my ($piecefrom,$pieceto) = ($1, $2);
	next unless ($from<$pieceto && $from+$count>$piecefrom);
	    
	my ($pieceoffset,$offset,$piececount) = ($from-$piecefrom,0,$count);
	if ($from<$piecefrom) {
	    $offset = $piecefrom-$from;
	    $pieceoffset = 0;
	    $piececount -= $offset;
	}
	if ($from+$count>$pieceto) {
	    $piececount = $pieceto-$piecefrom-$pieceoffset;
	}
 	substr($result, $offset, $piececount) = substr($$target{$_},$pieceoffset,$piececount);
    }
    return $complement ? 
	&revcom($result) : $result;
}

sub remove_user {
    my ($targetname, $ID) = @_;
#    my $oldkeycount = keys %{$dbase{$targetname}{users}};
    delete $dbase{$targetname}{users}{$ID};
    my $newkeycount = keys %{$dbase{$targetname}{users}};
    unless ($newkeycount) {
	delete $dbase{$targetname};
#	print COMMENT "Sequence '$targetname' discarded - not needed anymore.\n";
    } # elsif ($oldkeycount == $newkeycount) {
# 	print COMMENT "No user removed for '$targetname'\n";
#     } else {
# 	print COMMENT "$newkeycount users remain for '$targetname'\n";
#     }
}
    
sub get_intron_offset {   
    my ($gap_aa, $dss_seq, $ass_seq, $frameshift_mode, $preferred_offset) = @_;
    my $gaplen=length($gap_aa)*3;
    my ($bestoffset,$bestmismatchcount,$gap_trans,$splicesite_score) = (0,length($gap_aa),"-" x $gaplen,0); 
    
    my $ssb = $frameshift_mode ? 0 : 2; # in intron mode, account for splice sites
    
    my ($leftend, $rightend) = ($gaplen+$ssb-length($ass_seq), length($dss_seq)-$ssb);
    $leftend=0 if ($leftend < 0);
    $rightend=$gaplen if ($rightend > $gaplen);
    foreach my $offset  (reverse $leftend..$rightend) {
	next if ($frameshift_mode && $offset % 3);
	my $dna_comp = substr($dss_seq, 0, $offset);
	my $iscore_cand;
	$dna_comp .= substr($ass_seq, $offset-$gaplen) if ($offset < $gaplen);
	$dna_comp =  $ttable->translate($dna_comp);


	unless ($frameshift_mode) {
	    my $ipattern_cand = substr($dss_seq, $offset, 2).substr($ass_seq, $offset-$gaplen-2, 2);
	    $iscore_cand 
		= { "gtag" => 2.0, "gcag" => 1.8, "atac" => 0.9, "gaag" => 0.8, "ggag" => 0.7 }->{$ipattern_cand};
	    unless (defined $iscore_cand) {
		$iscore_cand = $ipattern_cand =~ /ag$/ ? 0.5 :
		    $ipattern_cand =~ /^g[tc]/ ? 0.3 : 0;
	    }
	}

	my $mismatchcount = &str_diff_count_mismatches($dna_comp, $gap_aa); #no of mismatches between dna_comp and gap
	$mismatchcount += ($dna_comp =~ /[*]/) * $IFS_MALUS;  # punishment for stop codons
	my $optimiser = $mismatchcount - $bestmismatchcount;
	$optimiser += ($splicesite_score - $iscore_cand) unless ($frameshift_mode);
	if ($optimiser < 0 || ($optimiser == 0 && $offset >= 2*$preferred_offset - $bestoffset)) {
	    ($bestoffset,$bestmismatchcount,$gap_trans,$splicesite_score)
		=($offset,$mismatchcount,$dna_comp,$iscore_cand);
	    last if ($mismatchcount==0 && ($frameshift_mode || $iscore_cand >= 2.0) && $offset <= $preferred_offset);
	}
    }
   
    return ($bestoffset, $bestmismatchcount, $gap_trans, $splicesite_score);
}  ### end sub &get_intron_offset

sub calculate_diff_strings {
    my ($targetname, $queryname, $targetlocations, $querylocations, $complement) = @_;
    my $transcript = ""; 
    for (my $i=0; $i<@$targetlocations; $i+=2) {
	my ($from, $to) = @$targetlocations[$i, $i+1]; 
	$transcript .= &subseq ($targetname, $from, $to-$from, $complement);
    }
    my $product = $ttable->translate($transcript);
    my @result = ();
    for (my $i=0; $i<@$querylocations; $i+=2) {
	my ($from, $to) = map { &nu_to_aa($_) } @$querylocations[$i, $i+1];
	my $partial_product = substr($product, 0, $to-$from, "");
	my $diff_str = &diff_str($partial_product, $queries{$queryname}, $to-$from, 0, $from);
	### DEBUG part ###
	unless (defined $diff_str) {
	    print STDERR "$targetname:$transcript\n";
	    &run_bug("diff_str returned undef!");
	}
	### end of DEBUG part ###
	push @result, "$partial_product:$diff_str";
    }
    return @result;
}
    

# Two hits on different targets are incompatible 
# - if they overlap more than MAX_OVERLAP
# Two hits on the same target are incompatible
# - if they overlap, or are in reversed order or on different strands
# 
sub incompatible {
    my ($hitref1, $hitref2) = @_;
    my $q1from = $hitref1->{"${QPFX}_start"};
    my $q2from = $hitref2->{"${QPFX}_start"};
    ($q2from, $hitref1, $hitref2) = ($q1from, $hitref2, $hitref1)
	if ($q1from > $q2from);
    my ($t1from, $sq1);
    ($q1from, $t1from, $sq1) = @{$hitref1->{matchings}[-1]}{"${QPFX}_start", "${TPFX}_start", "seqshifts"};
    my ($q2to,$t2to, $sq2) = @{$hitref2->{matchings}[0]}{"${QPFX}_end","${TPFX}_end", "seqshifts"};
    my ($q1to, $t1to, $target1, $comp1) = @$hitref1{"${QPFX}_end", "${TPFX}_end", "target", "complement"};
    my ($t2from, $target2, $comp2) = @$hitref2{"${TPFX}_start","target", "complement"};

    my ($length1, $length2) = map { $_->{length}} @dbase{$target1, $target2};
 
    # hits on the same target must not overlap and must be in correct order
    return 1 if (($target1 eq $target2) && ($q1to > $q2from  || $comp1 ne $comp2 || $t1to > $t2from));
 
    if ($q1to > $q2from) {
	foreach (@$sq1) { return 1 if ($_->{"${QPFX}_end"} >= $q2from) }
	foreach (@$sq2) { return 1 if ($_->{"${QPFX}_start"} <= $q1to) }
	my ($t1start, $t1end) = $comp1 ? (-$length1,0) : (0,$length1);
	my ($t2start, $t2end) = $comp2 ? (-$length2,0) : (0,$length2);
	my $offset =  3*($q1to-$q2from) - $t1to + $t2from;
	$t1start = max ( $t1start, $t2start - $offset);
	$t2start = $t1start + $offset; # max ( $t2start, $t1start + $offset);
    
	my $count = min ( $t1end - $t1start , $t2end - $t2start );
	if ($count <= $MAX_CHECK_OVERLAP) {
	    my $acceptable_errors = min ($count / 15, $count / 100 +1);
	    my ($subseq1, $subseq2) = (&subseq($target1, $t1start, $count, $comp1), &subseq($target2, $t2start, $count, $comp2));
	    if ( &str_diff_count_all($subseq1, $subseq2) <= $acceptable_errors  ) {
		push @{$hitref1->{target_overlaps}}, $hitref2->{ID};
		push @{$hitref2->{target_overlaps}}, $hitref1->{ID};
		$t1start++ unless ($comp1);
		$t2start++ unless ($comp2);
		print COMMENT 
		    "We found an overlap of target sequenzes:\n A: ${target1}\n B:$target2\n".
		     " A[".abs($t1start)."..".abs($t1start+$count-1)."] ".($comp1?"(comp)":"")." equals\n".
		     " B[".abs($t2start)."..".abs($t2start+$count-1)."] ".($comp2?"(comp)":"")."\n".
		return 0;
	    }
	}    
    }

    # no target overlap: ok if query overlap is no more than MAX_OVERLAP, 
    # and non-overlapping parts have at least the minimum length
    return 
	( $q1to  >  $q2from+$MAX_OVERLAP ||
	  $q1to > $q2to - $MIN_REMAIN ||
	  $q1from > $q2from - $MIN_REMAIN );

# 	if (  $q1to  <=  $q2from+$MAX_OVERLAP && 
# 	      $q1to <= $q2to - $MIN_REMAIN && 
# 	      $q1from <= $q2from - $MIN_REMAIN  );


}

sub replace_in_list {
    my ($from, $to, $listref, $newelems) = @_;
    my $oldsize = @$listref;
    @$listref = ( defined $from ? grep { $_ <= $from } @$listref : (), 
		  @$newelems,
		  defined $to ? grep { $_ > $to } @$listref : ());
    return @$listref - $oldsize;
}
    

sub replace_in_last {
    my ($prevhit, $dnaseq, $protseq, $dna_add, $nsub) = @_;
    $dna_add="" unless (defined $dna_add);
    $nsub=0 unless (defined $nsub);

    my $prevmatch = $prevhit->{matchings}[-1];
    my $oldend = $prevmatch->{"${QNPFX}_end"};
    my $newend = $oldend - $nsub;
    my $psub = &nu_to_aa($oldend) - &nu_to_aa($newend);
    
    my ($ncount, $pcount) = map  length ,($dnaseq, $protseq);
    my $translation = substr($ttable->translate($dnaseq.$dna_add),0, $pcount);
    $prevhit->{"${QPFX}_end"} += ($pcount - $psub);
    map { $_ += $ncount-$nsub } (@{$prevmatch}{"${QNPFX}_end","${TPFX}_end"}, 
				 $prevhit->{"${TPFX}_end"});
    substr(substr($prevmatch->{translation},-$psub-1),1) = $translation;
    my %newlists = ( mismatchlist => [ &str_diff_list('X',$translation, $protseq) ],
		     undeterminedlist => [ &str_diff_list(' ',$translation, $protseq) ],
		     inframe_stopcodons => [ &str_pos("[*]", $translation) ] );
    map { $_ += &nu_to_aa($newend) + 1 } @$_ foreach values %newlists;
    my (undef, $mismatchadd, $undetermadd) = 	map {
	&replace_in_list(&nu_to_aa($newend), undef, $prevmatch->{$_}, $newlists{$_}); 
    } sort keys %newlists;

    $prevhit->{mismatches} += $mismatchadd; 
    $prevhit->{undetermined} += $undetermadd;
    $prevhit->{matches} += ($pcount-$psub-$mismatchadd-$undetermadd);
}
    
sub replace_in_next {
    my ($nexthit, $dnaseq, $protseq, $dna_add, $nsub) = @_;
    $dna_add="" unless (defined $dna_add);
    $nsub=0 unless (defined $nsub);

    my $nextmatch = $nexthit->{matchings}[0];
    my $oldstart = $nextmatch->{"${QNPFX}_start"};
    my $newstart = $oldstart + $nsub;
    my $psub = &nu_to_aa($newstart) - &nu_to_aa($oldstart);
    
    my ($ncount, $pcount) = map  length ,($dnaseq, $protseq);
    my $translation = substr($ttable->translate($dna_add.$dnaseq), -$pcount, $pcount);

    $nexthit->{"${QPFX}_start"} -= ($pcount -$psub);
    map { $_ -= $ncount-$nsub } (@{$nextmatch}{"${QNPFX}_start","${TPFX}_start"}, 
				 $nexthit->{"${TPFX}_start"});
    substr($nextmatch->{translation},0,$psub)=$translation;
    my %newlists = ( mismatchlist => [ &str_diff_list('X',$translation, $protseq) ],
		     undeterminedlist => [ &str_diff_list(' ',$translation, $protseq) ],
		     inframe_stopcodons => [ &str_pos("[*]", $translation) ] );
    map { $_ += &nu_to_aa($newstart) - $pcount + 1 } @$_ foreach values %newlists;
    my (undef, $mismatchadd, $undetermadd) = 	map {
	&replace_in_list(undef, &nu_to_aa($newstart),  $nextmatch->{$_}, $newlists{$_}); 
    } sort keys %newlists;

    $nexthit->{mismatches} += $mismatchadd;
    $nexthit->{undetermined} += $undetermadd;
    $nexthit->{matches} += ($pcount-$psub-$mismatchadd-$undetermadd);
}

sub add_initial {
    my ($hit, $type) = @_;
    my $first_match = $hit->{matchings}[0];
    my $targetname = $hit->{target};
    my $complement = $hit->{complement};
    $hit->{"${TPFX}_start"} = $complement ? -$dbase{$targetname}{length} : 0;
    unshift @{$hit->{matchings}}, { "${TPFX}_start" => $hit->{"${TPFX}_start"},
				    "${TPFX}_end" => $first_match->{"${TPFX}_start"},
				    type => $type,
				    "${QNPFX}_start" => $first_match->{"${QNPFX}_start"},
				    "${QNPFX}_end" => $first_match->{"${QNPFX}_start"} }
}

sub add_final {
    my ($hit, $type) = @_;
    my $last_match = $hit->{matchings}[-1];
    my $targetname = $hit->{target};
    my $complement = $hit->{complement};
    $hit->{"${TPFX}_end"} = $complement ? 0 : $dbase{$targetname}{length};
    push @{$hit->{matchings}}, { "${TPFX}_start" => $last_match->{"${TPFX}_end"},
				 "${TPFX}_end" => $hit->{"${TPFX}_end"},
				 type => $type,
				 "${QNPFX}_start" => $last_match->{"${QNPFX}_end"},
				 "${QNPFX}_end" => $last_match->{"${QNPFX}_end"} }
}

sub set_problem_field {
    my $current = shift;
    my @problems;
    my ($qfrom, $qto, $qlen) = @$current{"${QPFX}_start","${QPFX}_end","${QPFX}_len"};
    foreach my $bad_type  (grep { $_ ne "exon" && $_ ne "intron" } 
			    map { $_->{type} } @{$current->{matchings}}) {
	$bad_type = $bad_type eq "intron?" ? "bad intron" : $bad_type;
	push @problems, $bad_type unless (grep { $_ eq $bad_type } @problems);
    }
    push @problems, "sequence shift" 
	if (grep  { (defined $_->{seqshifts} && @{$_->{seqshifts}}) } 
	    @{$current->{matchings}}); 
    push @problems, "mismatches" if ($current->{mismatches});
#   unnecessary because this implies a frameshift or a gap
#    push @problems, "unmatched" if ($current->{unmatched});
    push @problems, "missing stopcodon" if ($qto == $qlen && !$current->{stopcodon});
    push @problems, "in-frame stopcodon" 
	if (grep  { (defined $_->{inframe_stopcodons} && @{$_->{inframe_stopcodons}}) } 
	    @{$current->{matchings}}); 
    push @problems, "gap to querystart" if ($current->{upstream} && $qfrom>0);
    push @problems, "gap to previous hit" if (defined $current->{upstream_gap});
    push @problems, "gap to next hit" if (defined $current->{downstream_gap});
    push @problems, "gap to queryend" if ($current->{downstream} && $qto<$qlen);
#   add unmatched only unless there are other problems
#   in most cases unmatched implies a frameshift or a gap
    unless (@problems || !$current->{unmatched}) {
	push @problems, "unmatched" if ($current->{unmatched});
    }
    if (@problems) {
	$current->{status} = "incomplete";
	$current->{reason} = join "/",@problems;
    } elsif ($qfrom>0 || $qto < $qlen) {
	$current->{status} = "partial";
    } else {
	$current->{status} = "auto";
    }
}

################################################################################
