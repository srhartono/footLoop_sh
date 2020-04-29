#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use colorz;
use vars qw($opt_n);
getopts("n:");

my ($footPeak_folder) = ($opt_n);
die "\nusage: $YW$0$N -n $CY<foorPeak folder>$N\n\n" unless defined $opt_n and -d $opt_n;

my %coor;
die "Cannot find footPeak_logFile.txt from $LGN$footPeak_folder$N\n\n" if not -e "$footPeak_folder/footPeak_logFile.txt";
my @coors = `grep "def=" $footPeak_folder/footPeak_logFile.txt`;
foreach my $coor (@coors) {
	chomp($coor);
	my ($chr, $beg, $end, $gene, $val, $strand) = $coor =~ /^def\=.+, coor=(.+), (\d+), (\d+), (.+), (\d+), ([\+\-\.])$/;
	die "can't parse chr beg end etc from coor=$LCY$coor$N at file=\n$LPR$footPeak_folder/footPeak_logFile.txt$N\n\n" if not defined $chr;
	$coor{uc($gene)} = $strand;
}
my @peakFiles = <$footPeak_folder/PEAKS_LOCAL/*bed>;

foreach my $peakFile (sort @peakFiles) {
	my $data = parse_peakFile($peakFile);
	parse_callFile($peakFile, $data);

}
my %data;
sub parse_peakFile {
	my ($input1) = @_;
	my ($fileName1) = $input1 =~ /PEAKS_LOCAL\/(.+).local.bed$/; die "Can't parse fileName from input1=$LCY$input1$N\n\n" if not defined $fileName1;
	my $peakFile = "$footPeak_folder/PEAKS_LOCAL/$fileName1.local.bed";
	die "Cannot find peak local file $LCY$peakFile$N!\n\n" if not -e $peakFile;
	open (my $in1, "<", $peakFile) or die "Cannot read from $peakFile: $!\n";	
	my $data;
	while (my $line = <$in1>) {	
		chomp($line);	
		next if $line =~ /^#/;
		my ($readname, $beg, $end) = split("\t", $line);
		$data->{$readname}{"$beg,$end"} = 1;
	}
	close $in1;
	return $data;
}
sub parse_callFile {
	my ($input1, $data) = @_;
	system("mkdir -p $footPeak_folder/CONV/") if not -d "$footPeak_folder/CONV/";
	my ($fileName1) = $input1 =~ /PEAKS_LOCAL\/(.+)\.\w+\.local.bed$/; die "Can't parse fileName from input1=$LCY$input1$N\n\n" if not defined $fileName1;
	my $callFile = "$footPeak_folder/.CALL/$fileName1.PEAK.out";
	my ($name) = $fileName1 =~ /^(PCB.+desc\w+_gene\w+)_(Pos|Neg|Unk).+$/; $name = $fileName1 if not defined $name;
	$input1 = $callFile;
	my ($PCB, $BC, $plasmid, $desc, $gene, $strand, $window, $thres, $cpg) = $fileName1 =~ /(PCB\d+)_bcBC(\d+)_plasmid(.+)_desc(.+)_gene(.+)_(Pos|Neg|Unk)_(\d+\.?\d*)_(\d+\.?\d*)_([A-Z][A-Z][A-Z]?)\.?/;
	$thres = int(100*$thres);
	my $type = "w${window}t$thres";
	my $goodstrand = $coor{uc($gene)}; die "Cannot get goodstrand from gene=$gene fileName=$LGN$fileName1$N!\n\n" if not defined $goodstrand;
	$goodstrand = $goodstrand eq "+" ? "Pos" : $goodstrand eq "-" ? "Neg" : die "goodstrand isn't + or - ($LGN$goodstrand$N for gene=$LCY$gene$N file=$LPR$input1$N\n\n";
	my $cpgtype = "UNK";
	if ($goodstrand eq "Pos") {
		if ($strand eq "Pos") {
			$cpgtype = "PEAK" if $cpg eq "CH";
			$cpgtype = "PEAK_C" if $cpg eq "CG";
			$cpgtype = "PEAK_RCONV" if $cpg eq "GH";
			$cpgtype = "PEAK_RCONV_C" if $cpg eq "GC";
		}
		elsif ($strand eq "Neg") {
			$cpgtype = "PEAK_TEMP_RCONV" if $cpg eq "CH";
			$cpgtype = "PEAK_TEMP_RCONV_C" if $cpg eq "CG";
			$cpgtype = "PEAK_TEMP" if $cpg eq "GH";
			$cpgtype = "PEAK_TEMP_C" if $cpg eq "GC";
		}
	}
	elsif ($goodstrand eq "Neg") {
		if ($strand eq "Neg") {
			$cpgtype = "PEAK" if $cpg eq "GH";
			$cpgtype = "PEAK_C" if $cpg eq "GC";
			$cpgtype = "PEAK_RCONV" if $cpg eq "CH";
			$cpgtype = "PEAK_RCONV_C" if $cpg eq "CG";
		}
		elsif ($strand eq "Pos") {
			$cpgtype = "PEAK_TEMP_RCONV" if $cpg eq "GH";
			$cpgtype = "PEAK_TEMP_RCONV_C" if $cpg eq "GC";
			$cpgtype = "PEAK_TEMP" if $cpg eq "CH";
			$cpgtype = "PEAK_TEMP_C" if $cpg eq "CG";
		}
	}
	my $output = "$footPeak_folder/CONV/$fileName1.PEAK.out.$cpgtype.readconv";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";	
	open (my $out1, ">", "$output") or die "Cannot write to $output: $!\n";
#	print $out1 "#PCB\tBC\tplasmid\tdesc\tgene\treadname\tCpG\tperc_conv\tconv\ttotal_c\tpeak_length\n" if defined $PCB;
#	print $out1 "#filename\t.\t.\t.\t.\treadname\tCpG\tperc_conv\tconv\ttotal_c\tpeak_length\n" if not defined $PCB;
	print $out1 "#filename\ttype\tcpgtype\treadname\tperc_conv\tconv\ttotal_c\tpeak_length\n";
#	my $cpg = $fileName1 =~ /(CG|GC)/ ? "Y": "N";
	while (my $line = <$in1>) {	
		chomp($line);	
		next if $line =~ /^#/;
		my ($readname, @conv) = split("\t", $line);
		my $length = @conv;
		my $print = "";
		my $peakcoors = $data->{$readname};
		die "Cannot find peakcoor for read $readname!\n\n" if not defined $peakcoors;
		my %res = ("con"=>0,"conp"=>0,"notp"=>0,"not"=>0,"tot"=>0,"all"=>0,"len"=>0);
		foreach my $peakcoor (sort keys %{$peakcoors}) {
			my ($con, $not, $tot, $all) = (0,0,0);
			my ($beg, $end) = split(",", $peakcoor);
			my $length2 = $end - $beg;
			$print .= "$readname\tlength=$length, beg=$beg, end=$end\t";
			for (my $i = $beg; $i < $end; $i++) {
				$print .= "$conv[$i]";
				$con ++ if $conv[$i] =~ /^[89]$/;
				$not ++ if $conv[$i] =~ /^[45]$/;
				$tot ++ if $conv[$i] =~ /^[4589]$/;
				$all ++ if $conv[$i] !~ /^[4589]$/;
				die "Thre is a conv, but not peak (6/7) in $print!\n\n" if $conv[$i] =~ /^[67]$/;
			}
			my $conp = int($con/$tot*1000+0.5)/10;
			my $notp = int($not/$tot*1000+0.5)/10;
			$print .= "\ncon=$LGN$con$N $conp \%, not=$LCY$not$N $notp \%, tot=$LGN$tot$N, all=$YW$all$N, length=$LPR$length$N, lengthchunk=$LGN$length2$N\n";
#			print "$print\n";
			die "length chunk $length2 isn't the same as total+non C ($tot+$all=" . ($tot+$all) . ")\n\n" if $length2 != $tot+$all;
			$res{con} += $con;
			$res{not} += $not;
			$res{tot} += $tot;
			$res{all} += $all;
			$res{len} += $length2;
		}
		$res{conp} = $res{tot} == 0 ? 0 : int($res{con}/$res{tot}*1000+0.5)/10;
		$res{notp} = $res{tot} == 0 ? 0 : int($res{not}/$res{tot}*1000+0.5)/10;
##		print $out1 "$PCB\t$BC\t$plasmid\t$desc\t$gene\t$readname\t$cpg\t$res{conp}\t$res{con}\t$res{tot}\t$res{len}\n" if defined $PCB;
#		print $out1 "$fileName1\t.\t.\t.\t.\t$readname\t$cpg\t$res{conp}\t$res{con}\t$res{tot}\t$res{len}\n" if not defined $PCB;
		print $out1 "$name\t$type\t$cpgtype\t\"$readname\"\t$res{conp}\t$res{con}\t$res{tot}\t$res{len}\n";
	}
	close $in1;
	close $out1;
	print "OUTPUT:\n$LGN$output$N\n";
}
__END__
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

# 0 is bad or not data (6)
# 1 is non C/G
# 4 is CH non conv
# 5 is CG non conv
# 6 is CH conv
# 7 is CG conv
# 8 is CH peak
# 9 is CG peak
