#!/usr/bin/perl -w
################################################################################
#
#   Name:    SCIPIO - Eukaryotic Gene Identification
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2007-11-29
#   Version: accompanying SCIPIO 1.0
#
#
#   yaml2log.pl
#
#
#   This script transforms the yaml output of Scipio into a log file
#   showing the alignment
#
#
use strict;
use List::Util('sum', 'max', 'min');
use YAML;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Getopt::Long;

my $WIDTH=98;
my $TR_TABLE=1;
### debug mode
my $SEQSHIFTNAME="seqshifts";
### end of debug mode

# my $QUERYFILE;


my $QPFX="prot";
my $TPFX="dna";
my $TRPFX="trans";
my $QNPFX="nucl";

my $SHOW_GAPTRANS=0;

&GetOptions ("transtable=i" => \$TR_TABLE,
	     "seqshifts=s" => \$SEQSHIFTNAME,
	     "gaptrans" => \$SHOW_GAPTRANS,
	     "width=i" => \$WIDTH);

# my %queries;
# if (defined $QUERYFILE) {
#     my $querySeq;
#     unless (-f $QUERYFILE && ($querySeq = Bio::SeqIO->new(-file=>$QUERYFILE, -format=>"fasta"))) {
# 	print STDERR "Please specify query file in FASTA format!\n";
# 	exit 1;
#     }
#     while (my $seq = $querySeq->next_seq()) {
# 	$queries{$seq->id()} = $seq->seq();
#     }
# }
my %examplelist = ( 
    1 =>  "mismatch",
    2 =>  "undetermined query",
    3 =>  "undetermined target",
    4 =>  "additional codon in target",
    5 =>  "unmatched query",
    6 =>  "frameshift (+1) target only",
    7 =>  "frameshift (+1) target/query",
    8 =>  "frameshift (+2) target only",
    9 =>  "frameshift (+2) target/query",
    10 => "frameshift (-2) target only",
    11 => "frameshift (-2) target/query",
    12 => "frameshift (-1) target only",
    13 => "frameshift (-1) target/query",
    14 => "stopcodon target/query",
    15 => "stopcodon, target only",
    16 => "stopcodon, undetermined query",
    17 => "additional stopcodon");


my $ttable = Bio::Tools::CodonTable->new( -id => $TR_TABLE );

sub nu_to_aa {
    my $arg = shift;
    return undef unless (defined $arg);
    my $residue = ($arg+1) % 3 -1;   # -1, 0, 1
    return ($arg-$residue)/3;
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
	# count 'X' in translation always as unmatched
	# ('X' in query can be mismatch but never match)
	$result .= ($comp =~ /^[X-]$/ || /^[X-]$/) ?  " " : ( $comp eq $_ ? "|" : "X" );
	
    }
    return $result;
}

# sub escape {
#     $_=shift if (@_);
#     s/([^a-zA-Z0-9.:^*$@!+_?-|])/sprintf("%%%x",ord($1))/e;
# }

sub alignments
{
    my ($exon, $total_protseq) = @_;
    my ($tseq, $translation, $toffset, $qnoffset) 
	= @$exon{"seq", "translation", "${TPFX}_start", "${QNPFX}_start"};
    my ($qfrom, $qto) = map { &nu_to_aa($_); } @$exon{"${QNPFX}_start", "${QNPFX}_end"};
    my $qseq = $exon->{"${QPFX}_seq"};
    if (!defined $qseq && defined $total_protseq) {
	$qseq = substr($total_protseq, $qfrom, $qto-$qfrom);
    }
    $qseq = "-" x ($qto-$qfrom) unless (defined $qseq);
    my $trans_delta = length($translation)-($qto-$qfrom);
    $tseq = uc $tseq;
    if (defined $exon->{$SEQSHIFTNAME}) {
	foreach (reverse @{$exon->{$SEQSHIFTNAME}}) {
	    my ($tfrom, $qnfrom, $tto, $qnto) = @$_{"${TPFX}_start","${QNPFX}_start", "${TPFX}_end", "${QNPFX}_end"};
	    my $tshift = $tto - $tfrom;
	    my $qlen = $qnto - $qnfrom;
	    my $tgap = -$tshift % 3;
	    my $translen = ($tshift+$tgap)/3;
	    my $qpos = &nu_to_aa($qnfrom) - &nu_to_aa($qnoffset);
	    $trans_delta -= ($translen-$qlen/3);
	    if ($qlen < $tshift) {   # extra nucleotides in DNA
		if ($qlen) {
		    print STDERR "Invalid YAML entry: Invalid Frameshift!\n"; return (undef) x 4;
		}
		substr($qseq, $qpos, 0) = "-" x $translen;
	    } else	{   # missing nucleotides in DNA
		$tgap = $qlen - $tshift;
		substr($translation, $qpos + $trans_delta + $translen, 0) = "-" x ($qlen/3-$translen);
	    }
	    substr($tseq, $tfrom - $toffset, $tshift) = lc substr($tseq, $tfrom-$toffset, $tshift).("-" x $tgap);
	}
    }
    if ($trans_delta) {
	print STDERR "$translation\n$qseq\n";
	print STDERR "Invalid YAML entry: Frameshifts do not correspond to translation length\n";
	return (undef) x 4;
    }
    substr($tseq,0,0) = " " x ($qnoffset % 3);
    my $diffstr = &diff_str($translation, $qseq);
    if ($qnoffset % 3 == 2) { foreach ($translation, $diffstr, $qseq) { substr($_,0,0) = " "; } }
    s/(.)/ $1 /g foreach ($diffstr, $translation, $qseq);
    return ($tseq, $translation, $diffstr, $qseq);
}

sub pretty_sequence {
    my ($seq, $from, $linepref, $show_trans) = @_;
    my $result = "";
    $linepref = "" unless($linepref);
    my $indent = $from % 10;
    my $count = ($WIDTH-19-length($linepref))/11;
# next line enables last index counting mode for reverse strand
#    $from-- if ($from<0);

    my @trans = map { $show_trans && length($seq) > $_ ? 
			  $ttable->translate(substr("$seq",$_)) : "" } (0,1,2);
    substr($trans[2],0,0)=" ";
    while ($seq) {
    	my $line = substr($seq,0,60,"");
	my @trline = $show_trans ? map { substr($trans[$_],0,20,"") } (0,1,2) : ();
	if (@trline) {
	    s/(.)/$1  /g foreach (@trline); 
	    substr($trline[$_-1],0,0) = " " x $_ foreach (1,2);
	}
	$from += length($line);
	foreach my $i (60, 50, 40, 30, 20, 10) {
	    substr($line,$i-$indent,0)=" " if ($i-$indent<length($line));
	    foreach (@trline) {
		substr($_, $i-$indent,0)=" " if ($i-$indent<length);
	    }
	}
	$result .= "$linepref$line".(" " x (66 - length($line))).sprintf('%10d', abs($from))."\n";
	$result .= "$linepref$_\n" foreach (@trline);
    }
    return $result;
}

sub pretty_alignment {
# to do: use variable WIDTH here

    my ($tseq, $translation, $diffstr, $qseq, $tfrom, $qfrom, $linepref) = @_;
    return undef unless defined $tseq;
    my $result="";
    $linepref="" unless ($linepref);
# next line enables last index counting mode for reverse strand
#    $tfrom-- if ($tfrom<0);

    while ($tseq) {
	my @lines = 
	    map { sprintf('%-66s',substr($_,0,66,"")); } ($tseq, $translation, $diffstr, $qseq);
	$result .= $linepref.$lines[0]; 
	$lines[0] =~ s/[ \-]//g;
	$tfrom += length($lines[0]);
	$result .= sprintf('%10d', abs($tfrom))."\n$linepref";
	$result .= join((" "x9)."|\n$linepref", @lines[1..3]);
	$lines[3] =~ s/[ \-]//g;
	$qfrom += length($lines[3]);
	$result .= sprintf('%10d', $qfrom)."\n$linepref\n";
    }
    return $result;
}


sub pretty_print_matchings {
    my ($queryname, $hitref, $queryseq) = @_;
    my ($targetname, $strand, $upfrom, $downfrom, 
	$protseq, $protstart, $protlen)
	= @$hitref{"target", "strand", "${TPFX}_start", "${TPFX}_end", 
		   "${QPFX}_seq",  "${QPFX}_start", "${QPFX}_len"};
    
    my $result = "";
    $protstart = 0 unless defined $protstart;
    foreach ("upstream","upstream_gap") {
	my $seq = $hitref->{$_};
	$result.="$_\n".&pretty_sequence($seq, $upfrom-length($seq), "  ") if (defined $seq);
    }
    foreach (@{$hitref->{matchings}}) {
	my ($tfrom, $tto) = @$_{"${TPFX}_start","${TPFX}_end"};
	if ($_->{type} eq "exon" || $_->{type} eq "gap" ) {
	    my ($qfrom, $qto) = 
		defined $_->{"${QPFX}_start"} ?  
		@$_{"${QPFX}_start", "${QPFX}_end"} : 
		map { &nu_to_aa($_) } @$_{"${QNPFX}_start", "${QNPFX}_end"};
	    $qto = $qfrom unless (defined $qto);
	    if (defined $_->{protseq}) {
		substr($queryseq, $qfrom, $qto-$qfrom) = $_->{protseq};
	    }
	    if ($_->{type} eq "exon") {
		my $s = &pretty_alignment(&alignments($_, $queryseq), $tfrom, $qfrom, "  ");
		return "[invalid]\n" unless defined $s;
		$result.="exon\n$s";
	    } else { # type eq "gap"
		$result.="gap\n".&pretty_sequence($_->{seq}, $tfrom,"  ", $SHOW_GAPTRANS);
		if ($SHOW_GAPTRANS) {
		    my $qseq = substr($queryseq, $qfrom, $qto-$qfrom);
		    while ($qseq) {
			my $s = substr($qseq,0,66,"");
			$qfrom += length($s);
			$result.=sprintf("  %-66s%10d\n",$s, $qfrom);
		    }
		}
	    } 
	} else { # intron / frameshift
	    $result .= $_->{type}."\n".&pretty_sequence($_->{seq}, $tfrom, "  ");
	}
    }
    my $overlap = $hitref->{matchings}[-1]{overlap};
    if (defined $overlap) {
	$result.= "... (cut off overlap of $overlap bps; continued on next hit)\n";
    }
    foreach ("downstream","downstream_gap") {
	my $seq = $hitref->{$_};
	$result.="$_\n".&pretty_sequence($seq, $downfrom, "  ") if (defined $seq);
    }
    
    return $result;
}
###
  
sub querypos2dnapos {
    ## $qnpos must not be inside a seqshift!
    ## $qnpos is absolute (starting with query start), result is
    ## relative to matching start;
    my ($qnpos, $matching) = @_;
    return $qnpos - $matching->{"${QNPFX}_start"} + 
	sum ( 0, map { ($_->{"${TPFX}_end"}-$_->{"${TPFX}_start"}) 
			- ($_->{"${QNPFX}_end"}-$_->{"${QNPFX}_start"}) }
	      grep { $qnpos >= $_->{"${QNPFX}_end"} } @{$matching->{seqshifts}});
}

sub dnapos2transpos {
    ## tpos must not be inside a seqshift!
    ## tpos and result are relative to matching start
    my ($tpos, $matching) = @_;
    my $offset = 0;
    return &nu_to_aa($tpos + 
		     sum ( 0, map { ($_->{"${TPFX}_start"} - $_->{"${TPFX}_end"})%3 }
			   grep { $tpos + $matching->{"${TPFX}_start"} >= $_->{"${TPFX}_end"} }
			   @{$matching->{seqshifts}}) );
}

sub get_problem_overview {
    my ($queryname, $hitlist) = @_;
    my @problems = ();
    my $result = {};

    foreach (@$hitlist) {
	my ($queryseq, $querystart) = @{$_}{"${QPFX}_seq", "${QPFX}_start"};
	$querystart=0 unless defined $querystart;
	foreach (@{$_->{matchings}}) {
	    my ($qnfrom, $qnto, $seq, 
		$tfrom,  $translation, $type) = 
		( @{$_}{"${QNPFX}_start","${QNPFX}_end", "seq", 
			"${TPFX}_start", "translation","type"} );
	    my ($qfrom, $qto) = map ( &nu_to_aa($_), ($qnfrom, $qnto));
	    if ($type eq "gap") {
		$result->{$qfrom+1} = { 'U' => ($qto - $qfrom), case => 0 };
	    }
	    next unless ($type eq "exon");
	    $_->{seqshifts} = [] unless (defined $_->{seqshifts});

	    ### prepare list of mismatches, inframe stopcodons and undetermined
	    my @positions = map { 
		defined $_ ? @$_ : () 
	    } @{$_}{"mismatchlist","inframe_stopcodons","undeterminedlist"};
	    push @positions,
		map {
		    my ($qnfrom, $qnto) = @{$_}{"${QNPFX}_start","${QNPFX}_end"};
		    ($qnfrom==$qnto)?($qnfrom/3):();
		} @{$_->{seqshifts}};
	    foreach my $querypos (@positions)
	    {
		my $dnapos = &querypos2dnapos(($querypos-1)*3, $_);
		my $transpos = &dnapos2transpos($dnapos, $_);
		$result->{$querypos} = {
		    trans => substr($translation,$transpos,1),
		    seq => substr("..$seq..", $dnapos+2, 3),
		    qseq => substr($queryseq, $querypos-$querystart-1, 1)
		}
	    }
	    @{$result->{$_}}{"M","case"}=(1,1) foreach @{$_->{mismatchlist}};
	    @{$result->{$_}}{"U","case"}=(1, $result->{$_}{qseq} eq 'X' ? 2 : 3) 
		foreach @{$_->{undeterminedlist}};
	    foreach (map { $result->{$_} } @{$_->{inframe_stopcodons}}) {
		$_->{"*"}++;
		$_->{case} += 14;
	    }

	    ### prepare list of sequence shifts
	    foreach my $fs (@{$_->{seqshifts}}) {
		my ($qnfrom, $qnto, $fsstart) = @{$fs}{"${QNPFX}_start", "${QNPFX}_end", "${TPFX}_start"};
		my $qfrom = $qnfrom/3;
		my $fslen = $fs->{"${TPFX}_end"} - $fsstart;
		$fsstart -= $tfrom;
		my $transpos = &dnapos2transpos($fsstart, $_);
		my $current = $result->{$qfrom};
		if ($qnfrom != $qnto)  {
		    @{$result->{$qfrom+1}}{"trans","seq","qseq"} = ("") x3;
		    $current = $result->{$qfrom+1};
		}
		my @cases = ($current->{case});
		@cases = () unless (defined $cases[0]);
		my $queryX = (@cases && $cases[0]==2);
		@cases = () if ($qnfrom != $qnto);
		$current->{trans} .= substr($translation, $transpos, ($fslen+2)/3);
		$current->{seq} .= substr($seq, $fsstart, $fslen);
		$current->{qseq} .= substr($queryseq, $qfrom-$querystart, $qnto/3-$qfrom);
		my $additional = int(($fslen+2)/3) - ($qnto/3-$qfrom);
		if ($additional > 0) { $current->{"+"} = $additional }
		elsif ($additional < 0) { $current->{"-"} = -$additional }
		if ($fslen%3) { 
		    $current->{"F"} = $additional ? ($fslen %3) : -(-$fslen %3); 
		}
			      
		if ($fslen>=3) {
		    push @cases, $current->{"*"} ? 17 : 4 if ($fslen>=3);
		}
		if ($fslen % 3) {
		    push @cases, ($additional>0 ? 4 : 8) + ($fslen % 3)*2;
		    $cases[-1]++ if ($queryX);
		} 
		push @cases, 5 if ($additional < 0);
		$current->{case} = join ",",@cases;
	    }
	}
    }
    foreach my $current (values %$result) {
	$current->{counted_as} = 
	    join ",",  map { 
		$current->{$_} != 1 ? "$_(".$current->{$_}.")" : $_ 
	    } grep($current->{$_}, 
		   "M","-","U","+","*");
    }
    return $result;
}



sub output {
    # to make clear which coordinates we are referring to, we use prefixes:
    # 'q' is amino acid positions in the query
    # 'qn' is nucleotide postitions in the query
    # positions without prefix are always nucleotide positions

    my ($queryname, $hitlist) = @_;  
    my $overview = &get_problem_overview($queryname, $hitlist);
    my $result = "\n".("#"x 78)."\nquery".(" "x11)."$queryname\n";
    my @statuslist = map { $_->{status} } (@$hitlist);
    my $total_hits = @$hitlist;
    my $status = (grep(/incomplete/, map($_->{status},@$hitlist))?"in":"")."complete";
    $result .= "status".(" "x10)."$status\n";

    my $count=0;
    my ($qlen, $q0from, $q0to) = @{$hitlist->[0]}{"${QPFX}_len","${QPFX}_start","${QPFX}_end"};
    $q0to = $qlen unless defined $q0to;
    my $queryseq = "." x $qlen;
    foreach ( @$hitlist) {
	my ($protseq, $protstart) = @{$_}{"${QPFX}_seq", "${QPFX}_start"};
	$protstart = 0 unless defined $protstart;
	$protseq = "" unless defined $protseq;
	substr($queryseq, $protstart, length($protseq)) = $protseq;
    }    
    if ($q0from || $q0to < $qlen ) {
	$result .= ("-" x 78)."\n";
	$result .= "partial hit".(($total_hits>1)?"s on $total_hits different targets":"").":\n";
	my $qto = 0;
	foreach (@$hitlist) {
	    $count++;
	    my ($qfrom, $targetname) = @{$_}{"${QPFX}_start","target"};
	    $qfrom = 0 unless defined $qfrom;
	    if ($qfrom > $qto) {
		$result .= 
		    sprintf("%9d", $qto + 1) .
		    ($qfrom > $qto +1 ? sprintf("..%-5d", $qfrom) : " " x 7).
		    "missing\n";
	    }
	    $qto = $_->{"${QPFX}_end"};
	    $qto = $_->{"${QPFX}_len"} unless defined $qto;
	    $result .= sprintf("(#$count)%5d..%-5dfound on target $targetname\n", $qfrom+1, $qto);
	}
	my $qend = $hitlist->[0]{"${QPFX}_len"};
	if ($qend > $qto) {
	    $result .= 
		    sprintf("%9d", $qto + 1) .
		    ($qend > $qto +1 ? sprintf("..%-5d", $qend) : " " x 7).
		    "missing\n";
	}
    }
    $result.=("-"x78)."\nmismatches and sequence shifts:\n" if (keys %$overview);
    foreach (sort {$a <=> $b} keys %$overview) {
	if ($overview->{$_}{case} eq "0") {
	    $result.=sprintf("%9d..%-4d: missing, %s\n", 
			     $_, $_ + $overview->{$_}{U} -1, $overview->{$_}{counted_as});
	} else {
	    my ($qseq, $trans, $seq, $fs, $count_as, $cases) = 
		@{$overview->{$_}}{"qseq","trans","seq","F", "counted_as","case"};
	    $seq =~ s/(...)/$1 /g;
	    $seq =~ s/ *$/]/;
	    $trans = "-" unless $trans;
	    $fs = $fs ? "F:$fs," : "";
	    my $description = join (" / ", @examplelist{split ( ",", $cases)});
	    $result.=sprintf ("  at%5d: expected:%-4sfound:%-4s[%-9s %-12s example:%s, %s\n",
			      $_, $qseq, $trans, $seq, $fs.$count_as, $cases, $description);
	}
    }
	
    $count=0;
    foreach my $hitref (@$hitlist) {
	$count++;
	
	my ($hit_id, $matchings, $strand, $targetname, $targetlength, $score, 
	    $gapcount, $undetermined, $insertions, $stopcodon, $status, $reason, $tfrom, $tto) = 
	    @$hitref{"ID", "matchings", "strand", "target", "target_len", "score",
		     "unmatched", "undetermined", "additional", "stopcodon", "status", "reason",
		     "${TPFX}_start", "${TPFX}_end"};
	my @problems = defined $reason ? split /\//, $reason : ();
	my $complement = ($strand eq "-");
	my @mismatches = map { @{$_->{mismatchlist}} } grep { defined $_->{mismatchlist} } @$matchings;
	my $total_mismatches = @mismatches;
	
	my ($querystart, $queryend, $total_length) = @$hitref{"${QPFX}_start", "${QPFX}_end", "${QPFX}_len"};
	$querystart = 0 unless (defined $querystart);
	$queryend = $total_length unless (defined $queryend);
	$gapcount = $hitref->{"gaps"} unless (defined $gapcount);
	$gapcount = 0 unless (defined $gapcount);
	$undetermined = 0 unless (defined $undetermined);
	my $matchsize = $queryend-$querystart;
	my $matchcount = $matchsize - $gapcount - $undetermined - $total_mismatches;
	my @intronlocs = map ( [ @{$_}{"${QNPFX}_start","${TPFX}_start","${TPFX}_end"} ] ,
			       grep { $_->{type} =~ /^intron/ } @$matchings );
	my $introntlocs = join(",", map (  $complement ? 
					   (0-$_->[1])."..".(1-$_->[2]) : ($_->[1]+1)."..".($_->[2]), 
					   @intronlocs));
	my $intronqlocs = join(",", map { &nu_to_aa($_->[0]) } @intronlocs);
	$intronqlocs = "introns at      $intronqlocs\n" if ($intronqlocs);
	$introntlocs = "introns at      $introntlocs\n" if ($introntlocs);
	
	$result .= ("#"x78)."\n"
	    ."ID              $hit_id\n"
	    ."status        "
	    .($status eq "incomplete" ? "!" : " ")." $status\n"
	    .($status eq "incomplete" ? "reason        ! ".join("\n              ! ",@problems)."\n" : "" )
	    ."query length    $total_length\n"
	    .($matchsize == $total_length ? "" : 
	      "query location  ".($querystart+1)."..$queryend".($total_hits>1?" (#$count/$total_hits)\n" : "\n"))
	    .$intronqlocs
	    ."target          $targetname\n"
	    ."target length   $targetlength\n"
	    ."strand          $strand\n"
	    ."no of exons     ".(scalar grep { $_->{type} eq "exon" } @$matchings)."\n"
	    ."target location ".( $complement ? (0-$tfrom)."..".(1-$tto) : ($tfrom+1)."..$tto" )."\n"
	    .$introntlocs
	    ."hit length      $matchsize\n"
 	    .($matchcount == $matchsize ? "" :
	      "- matches     ! $matchcount\n")
	    .($total_mismatches == 0 ? "" :
	      "- mismatches  ! $total_mismatches\n")
	    .($gapcount == 0 ? "" :
	      "- unmatched   ! $gapcount\n")
	    .($undetermined == 0 ? "" :
	      "- undeterm.   ! $undetermined\n")
	    .($insertions == 0 ? "" :
	      "insertions    ! $insertions\n")
	    .($matchcount == $matchsize ? "" :
	      "identity      ! ".sprintf("%04.1f%%\n", $matchcount / $matchsize * 100))
	    ."score           $score\n"
	    .("-"x78)."\n";
	$result.=&pretty_print_matchings ($queryname, $hitref, $queryseq);
    }
    print $result;
}    
  


my @input  = split("\n---",join("",<>));
my $header = shift @input if (@input && $input[0]!~/^---/);

print "### yaml header:\n$header\n" if (defined $header);
foreach (@input) {
    s/^/---/;
    s/$/\n/;
    my $hits = YAML::Load($_);
   my @keys = sort keys %$hits;
#    my @keys = keys %$hits;
    foreach (@keys) {
	&output($_, $hits->{$_});
    }
}



