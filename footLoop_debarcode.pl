#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use File::Basename qw(dirname); use Cwd qw(abs_path);
use vars qw($opt_v $opt_b $opt_i $opt_d $opt_x $opt_l $opt_f $opt_o);
getopts("vb:i:d:xl:fo:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
   print "\n- Pushed $libPath into perl lib path INC\n";
}

use myFootLib; use FAlite; use footPeakAddon;

my $md5script = `which md5` =~ /md5/ ? "md5" : "md5sum";
my $homedir = $ENV{"HOME"};
my $footLoopScriptsFolder = dirname(dirname abs_path $0) . "/footLoop";
my @version = `$footLoopScriptsFolder/check_software.pl | tail -n 12`;
my $version = join("", @version);
if (defined $opt_v) {
   print "$version\n";
   exit;
}
my ($version_small) = "vUNKNOWN";
foreach my $versionz (@version[0..@version-1]) {
   ($version_small) = $versionz =~ /^(v?\d+\.\d+\w*)$/ if $versionz =~ /^v?\d+\.\d+\w*$/;
}

my $usage = "
-----------------
$YW $0 $version_small $N
-----------------

Usage: $YW$0$N ${LGN}[-x: Dry run]$N $CY-i<ccs_fq.fq.gz>$N $LPR-b <barcode.csv>$N $LGN-l<id number e.g. 1234>$N 

$LRD IMPORTANT!!!! $LCY
Label (-l) will be used as file name for$YW ALL resulting files$N, so make sure it's$YW unique and distinct$LCY from other samples.
We recommend to use sequencing date for ease of organization, e.g. -l 190204
Or use the first few digits of sequencing ID from fq file, e.g. m180515_233839_42145, then -l 180515 or -l 180515233$N

Optional: 
-o <output dir> [default: same as -l]
-d <name for debug e.g. 203>

";

print $usage unless defined $opt_i and -e $opt_i and defined $opt_b;
print "$LRD ERROR$N -i not defined\n" if not defined $opt_i;
print "$LRD ERROR$N -b not defined\n" if not defined $opt_b;
die "\n" if not defined $opt_i or not defined $opt_b or (defined $opt_i and not -e $opt_i) or (defined $opt_b and not -e $opt_b);
#die "$usage\n\nDoes not exist -b ${LPR}opt_b$N!\n" if not defined $opt_b;
#die "$usage\n\nDoes not exist -i $LCY$opt_i$N!\n" if defined $opt_i and not -e $opt_i;
die "\n\n" if not defined $opt_x and (not defined $opt_i or not defined $opt_b or (defined $opt_i and not -e $opt_i));
die "\n\n" if defined $opt_x and (not defined $opt_b);


my $input1 = $opt_i;

my $debugname = $opt_d;
$debugname = "NA" if not defined $opt_d;
my $barcodeFile = $opt_b;
my @barcodeFiles = split(",", $barcodeFile);

if (defined $opt_x and not defined $opt_i) {
	system("touch dummy.temp") == 0 or die "Failed to create $LGN dummy.temp$N: $!\n";
	$opt_i = "dummy.temp";
	$input1 = $opt_i;
}
my $label = defined $opt_l ? $opt_l : getFilename($input1);
my @labels = split("", $label);
$label = "";
for (my $i = 0; $i < @labels; $i++) {
	$label .= $labels[$i] if $labels[$i] =~ /^[0-9]$/;
}
my $date = myFootLib::getDate_simple();
if ($label eq "") {
	$label = $date;
}
my $num = 0;
$label = "PCB$label";
#my $label2 = "PCB$label";
#while (-d $label2) {
#	$num ++;
#	$label2 = "PCB$label" . $num;
#}
#$label = $label2;
my $dbfolder = defined $opt_o ? $opt_o : "$label";
my $dbfolderFull = myFootLib::getFullpath($dbfolder);
print "$LRD WARNING$N Output dir already exist (-o $dbfolder), printing there anyway\n\n" if -d $dbfolder;
mkdir $dbfolder if not defined $opt_x;
my $inputName = $label;
print "-b Barcode       : $LPR$opt_b$N\n" if defined $opt_b;
print "-i Input fastq   : $LCY$opt_i$N\n" if defined $opt_i;
print "-l Label         : $LRD$label$N\n";
print "-o Output dir    : $dbfolderFull\n";
my %outfq;

my $outLog;
open ($outLog, ">", "$inputName\_logfile.txt") or die "Cannot write to $inputName\_logfile.txt: $!\n" if not defined $opt_x;
# parse barcode
my (%bc, %out);
my @line = `cat @barcodeFiles`; my $linecount = -1;
LOG($outLog, "\nOutput:\n");
foreach my $line (@line[0..@line-1]) {
	chomp($line);
	$linecount ++;
	next if $line =~ /^\n$/;
	next if $line !~ /[a-zA-Z0-9]+/;
	next if $line =~ /^(\#|No)/i or $line =~ /\ttarget/i;
	#my $delim = $line =~ /\t/ ? "\\t" : $line =~ /,/ ? "," : $line =~ /[ ]/ ? " " : die "Cannot determine delimiter at line=$LCY\n$line\n$N\n";
	#LOG($outLog, "Delimiter = $delim\n");
	my @arr = split("[\t ]+", $line);
	my ($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;
	if (@arr == 2) {
		($bc, $bcseq, $desc, $plasmid) = @arr;
	}
	else {
		($no, $target, $size, $conc, $vol, $bc, $bcseq, $desc, $plasmid) = @arr;
	}
	$target 	= $bc if not defined $target;
	$desc 		= $bc if not defined $desc;
	$plasmid 	= $bc if not defined $plasmid;
	($desc, $bcseq, $bc, $plasmid, $target) = fixformat($desc, $bcseq, $bc, $plasmid, $target);
	($desc, $bc, $plasmid) = (uc($desc), uc($bc), uc($plasmid));
	my $bcbefore = $bc;
	$bc = "bc$bc\_plasmid$plasmid\_desc$desc";
	die "\n${LRD}FATAL ERROR$N: Multiple barcode name for the same barcode sequence.\n$LCY" . $bc{uc($bcseq)}{bc} . "$N\ncurrent=$LPR$bc$N\nbarcode=$LGN" . uc($bcseq) . "$N)\n\n" if defined $bc{uc($bcseq)};
	if (@arr == 2) {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		$bc{uc($bcseq)}{bc} = $bc;
		LOG($outLog, "  $linecount.\t$LCY$dbfolder$N/$LGN$inputName\_$bc.fq.gz$N\tbc=$LPR$bcbefore$N,plasmid=$LCY$plasmid$N,desc=$YW$desc$N\n");
#		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fq$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	else {
		$bc{uc($bcseq)}{target} = $target;
		$bc{uc($bcseq)}{bc} = $bc;
		$bc{uc($bcseq)}{total} = 0;
		$bc{uc($bcseq)}{line} = $linecount;
		LOG($outLog, "  - $linecount.\t$LCY$dbfolder$N/$LGN$inputName\_$bc.fq.gz$N\tbc=$LPR$bcbefore$N,plasmid=$LCY$plasmid$N,desc=$YW$desc$N\n");
#		print STDERR "\t$LGN$linecount\t$dbfolder/$inputName\_$bc.fq.gz$N\tbc=$LPR$bcbefore$N, plasmid=$LCY$plasmid$N, desc=$YW$desc$N\n";
	}
	if (not defined $opt_x) {
		open ($out{$bc}, ">", "$dbfolder/$inputName\_$bc.fq") or die "Cannot write to $dbfolder/$inputName\_$bc.fq: $!\n";
	}
	$outfq{"$dbfolder/$inputName\_$bc.fq"} = 1;
}
LOG($outLog, "  - $linecount\t$LCY$dbfolder$N/$LGN$label\_UNKNOWN.fq$N\tbc=${LPR}NA$N,plasmid=${LCY}NA$N,desc=${YW}NA$N\n");
$outfq{"$dbfolder/$label\_UNKNOWN.fq"} = 1;

if (defined $opt_x) {
	LOG($outLog, "\nDry run ${LGN}SUCCESS$N\n\n");
	exit;
}
open (my $outnobc, ">", "$dbfolder/$inputName\_UNKNOWN.fq") or die "Cannot write to $dbfolder/$inputName\_UNKNOWN.fq: $!\n";
my ($seq_nobc, $seq_dblbc, $seq_totbc, $readCount) = (0,0,0,0);
my (%data, %bad, @len, $in1);
my ($folder1, $fileName1) = getFilename($input1, "folder");
LOG($outLog, "${LGN}1. -b $LCY$barcodeFile$N was parsed successfully!\n$N\n");
# parse fq
my $isDir = -d $input1 ? "(is a directory!)" : -e $input1 ? "" : "($input1 does not exist!)";
die "$LRD FATAL ERROR$N: -i $LCY$input1$N is not a fq file! $LRD$isDir$N\n" if -d $input1 or not -e $input1;
open ($in1, "<", $input1) or die "Cannot read from $input1: $!\n" if $input1 !~ /.gz$/i;
open ($in1, "zcat < $input1|") or die "Cannot read from $input1: $!\n" if $input1 =~ /.gz$/i;
while (my $line = <$in1>) { chomp($line);
	# 1. name: m171023_184836_42145_c101402632550000001823289101251862_s1_p0/0/ccs
	$readCount ++;
	my ($readName) = $line;
	$line = <$in1>; chomp($line);
	# 2. seqs: ACTGTCG....
	my $len = int(length($line));
	my ($seq) = $line;
	
	$line = <$in1>; chomp($line);
	# 3. junk (+)

	$line = <$in1>; chomp($line);
	# 4. Quality
	my $qua = $line;

	my ($bc, $bcs, $seq2, $qua2) = infer_barcode($seq, $qua, $readName);
	print $outnobc "$readName\n$seq\n\+\n$qua\n" if not defined $seq2;
	print $outLog "\nPrev Seq:\n$seq\n$qua\n" if $readName eq $debugname;
	print {$out{$bc}} "$readName\n$seq2\n\+\n$qua2\n" if defined $seq2;

	LOG($outLog, date() . "\t$YW$input1$N: Done $LGN$readCount$N\n") if $readCount % 1000 == 0;
}
close $in1;
LOG($outLog, "\n\n${LGN}2. -i $LCY$input1$N was parsed successfully!\n$N\n");
my $double_bc = (keys %bad);
my $seq_goodbc = $seq_totbc - $seq_nobc - $double_bc;
my $seq_goodbcPerc = int($seq_goodbc / $seq_totbc * 1000)/10;
LOG($outLog, "

Total seq : $seq_totbc
Has bc    : $seq_goodbc ($seq_goodbcPerc \%)
Double bc : $double_bc
No bc     : $seq_nobc

");
foreach my $bcseq (sort {$bc{$a}{line} <=> $bc{$b}{line}} keys %bc) {
	print $outLog "\t$bc{$bcseq}{bc}\t$bcseq\t$bc{$bcseq}{total}\n";
	print STDERR "\t$bc{$bcseq}{bc}\t$bcseq\t$bc{$bcseq}{total}\n";
}

LOG($outLog, "\n\n${LGN}3. gzipping all fq files!$N\n\n");
foreach my $file (keys %outfq) {
	my $cmd = "gzip -f $file";
	LOG($outLog, "$LGN$cmd$N\n");
	system($cmd) == 0 or LOG($outLog, "Failed to$LGN $cmd$N: $!\n");
}

sub fixformat {
	my (@strings) = @_;
	for (my $i = 0; $i < @strings; $i++) {
		my @string = split("", $strings[$i]);
		$strings[$i] = "";
		for (my $j = 0; $j < @string; $j++) {
			$strings[$i] .= $string[$j] if $string[$j] =~ /^[0-9A-Za-z]$/;
		}
	}
	return @strings;
}

sub infer_barcode {
	my ($seq, $qua, $name) = @_;
	my $barcode; my $barcode2; my $type;
	my $seq1; my $seq2;
	my $seq_length = length($seq);
	my ($beg, $end) = 0;
	for (my $l = 0; $l < 4; $l++) {
		foreach my $bcseq2 (sort {$bc{$a}{bc} cmp $bc{$b}{bc}} keys %bc) {
			my $bcseq = substr($bcseq2, $l, length($bcseq2)-$l+1);
			my $length = length($bcseq);
			my $bcseqrev = $bcseq;
			$bcseqrev =~ tr/ACGT/TGCA/;
			$bcseqrev = reverse($bcseqrev);
			if (($seq =~ /^(.){0,10}($bcseq|$bcseqrev)/ or $seq =~ /($bcseq|$bcseqrev).{0,10}$/)) {
				my $pos = 0 if $seq =~ /^($bcseq|$bcseqrev)/;
				   $pos = $seq_length if $seq =~ /($bcseq|$bcseqrev)$/;
				my ($pos1, $pos2) = $seq =~ /^(.+)($bcseq|$bcseqrev)(.+)$/ if not defined $pos and $seq =~ /($bcseq|$bcseqrev)/;
					$pos = defined $pos ? $pos : length($pos1) <= 10 ? length($pos1) : length($pos2);
				my $type = $pos < 0.5 * $seq_length ? 1 : 2;
				die "$name: barcode=$barcode, bcseq=$bcseq, pos1=$pos1 pos2=$pos2\n" if not defined $pos;
				$barcode2 = $seq =~ /$bcseq/ ? "$type,$bcseq" : "$type,$bcseqrev";
				if (defined $barcode) {
					$bad{$name} = 1;
					$seq2 = $bcseq;
					print $outLog "$LCY$name\tDBL\t$pos/$seq_length\t$barcode\t$bc{$bcseq2}{bc}\t$seq1\t$seq2$N\n";
				}
				else {
					$barcode = $bc{$bcseq2}{bc};
					$seq1 = $bcseq;
					$data{$name} = $barcode;
					$bc{$bcseq2}{total} ++;
					print $outLog "$LGN$name\tGOOD\t$pos/$seq_length\t$barcode$N\n";
				}
			}
			die if $l > 3;
		}
		$seq_totbc ++ if $l == 0;
		if (defined $barcode and not defined $bad{$name}) {
			my $bcseq    = $seq1;
			my $bcseqrev = $bcseq;
			$bcseqrev =~ tr/ACGT/TGCA/;
			$bcseqrev = reverse($bcseqrev);
			my ($bad1, $good1) = $seq =~ /^(.{0,14}$bcseq|.{0,14}$bcseqrev)(.+)$/ if $seq =~ /^(.{0,14}$bcseq|.{0,14}$bcseqrev)/;
			my ($good2, $bad2) = $seq =~ /^(.+)($bcseq.{0,14}|$bcseqrev.{0,14})$/ if $seq =~ /^(.+)($bcseq.{0,14}|$bcseqrev.{0,14})$/;
			if (defined $bad1) {
				my $len = length($bad1);
				$seq =~ s/^.{$len}//;
				$qua =~ s/^.{$len}//;
			}
			if (defined $bad2) {
				my $len = length($bad2);
				$seq =~ s/.{$len}$//;
				$qua =~ s/.{$len}$//;
			}
			$bad1 = "" if not defined $bad1;
			$bad2 = "" if not defined $bad2;
			return ($barcode, $barcode2, $seq, $qua);
		}
	}
	print $outLog "$LRD$name\tNO BARCODE!$N\n" if not defined $barcode;
	$seq_nobc ++ if not defined($barcode);
	return;
}

__END__
