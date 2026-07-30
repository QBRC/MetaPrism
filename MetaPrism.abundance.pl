#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use IPC::Open2;
use List::Util qw(sum);
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

my $defaultBaseAbundanceGenes = 'K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965';
my $defaultNumberOfBaseAbundanceGenes = 10;

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 1000),
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
	'B=f' => \(my $baseAbundance = ''),
	'b=s' => \(my $baseAbundanceGenes = ''),
	'e=i' => \(my $ignoreOneSideOutlierNumberOfBaseAbundanceGenes = 2),
	'use_merged_bam_file' => \(my $use_merged_bam_file = ''),
	'use_samtools_depth' => \(my $use_samtools_depth = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism.abundance.pl [options] MetaPrism.txt genome.sorted.bam [...] > MetaPrism.abundance.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -f INT   include flag [$includeFlag]
         -F INT   exclude flag [$excludeFlag]
         -S STR   stranded, "f" or "r"
         -B FLOAT base abundance
         -b STR   base abundance genes or number of base abundance genes [$defaultBaseAbundanceGenes or $defaultNumberOfBaseAbundanceGenes if $dataPath/single_copy_gene.count.txt is available]
         -e STR   ignore one-side outlier number of base abundance genes [$ignoreOneSideOutlierNumberOfBaseAbundanceGenes]
         --use_merged_bam_file use a temporary merged bam file
         --use_samtools_depth use samtools depth instead of parsing cigar strings

EOF
}
if($baseAbundanceGenes eq '') {
	$baseAbundanceGenes = $defaultBaseAbundanceGenes;
	$baseAbundanceGenes = $defaultNumberOfBaseAbundanceGenes if(-s "$dataPath/single_copy_gene.count.txt" and -r "$dataPath/single_copy_gene.count.txt");
}
{
	my $parentPid = $$;
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryDirectory/fork.$hostname.$parentPid.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryDirectory/fork.$hostname.$parentPid.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryDirectory/fork.$hostname.$parentPid.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my $temporaryPrefix = "$temporaryDirectory/MetaPrism.$hostname.$$";
if(-r $baseAbundanceGenes) {
	my @baseAbundanceGeneList = ();
	open(my $reader, $baseAbundanceGenes);
	while(my $line = <$reader>) {
		chomp($line);
		push(@baseAbundanceGeneList, $line);
	}
	close($reader);
	$baseAbundanceGenes = join(',', @baseAbundanceGeneList);
} elsif($baseAbundanceGenes =~ /^[0-9]+$/) {
	my @baseAbundanceGeneList = ();
	open(my $reader, "$dataPath/single_copy_gene.count.txt") or die "Can't open '$dataPath/single_copy_gene.count.txt': $!";
	while(scalar(@baseAbundanceGeneList) < $baseAbundanceGenes) {
		chomp(my $line = <$reader>);
		my ($baseAbundanceGene) = split(/\t/, $line, -1);
		push(@baseAbundanceGeneList, $baseAbundanceGene);
	}
	close($reader);
	$baseAbundanceGenes = join(',', @baseAbundanceGeneList);
}
my ($metaPrismFile, @bamFileList) = @ARGV;
if($use_merged_bam_file || $use_samtools_depth) {
	if(scalar(@bamFileList) == 1) {
		system("ln -sfr @bamFileList $temporaryPrefix.bam");
	} else {
		system("samtools merge --threads $threads $temporaryPrefix.bam @bamFileList");
	}
	if($use_samtools_depth) {
		if($minimumMappingQuality == 0 && $includeFlag == 0 && $excludeFlag == 0) {
			system("ln -sfr $temporaryPrefix.bam $temporaryPrefix.filtered.bam");
		} else {
			system("samtools view -b -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $temporaryPrefix.bam > $temporaryPrefix.filtered.bam");
		}
		if($stranded) {
			open(my $reader, "samtools view -h $temporaryPrefix.filtered.bam |");
			open(my $writerF, "| samtools view -b - > $temporaryPrefix.filtered.f.bam");
			open(my $writerR, "| samtools view -b - > $temporaryPrefix.filtered.r.bam");
			my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ /^\@/) {
					print $writerF $line, "\n";
					print $writerR $line, "\n";
					next;
				}
				my %tokenHash = ();
				(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
				$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
				print $writerF $line, "\n" if(grep {($tokenHash{'flag'} & 253) == $_} 97, 145, 0);
				print $writerR $line, "\n" if(grep {($tokenHash{'flag'} & 253) == $_} 81, 161, 16);
			}
			close($reader);
			close($writerF);
			close($writerR);
			system("samtools index $temporaryPrefix.filtered.f.bam");
			system("samtools index $temporaryPrefix.filtered.r.bam");
		} else {
			system("samtools index $temporaryPrefix.filtered.bam");
		}
	} else {
		system("samtools index $temporaryPrefix.bam");
		@bamFileList = ("$temporaryPrefix.bam");
	}
}
my $use_sam_module = eval {require Bio::DB::Sam; 1;};
foreach my $bamFile (@bamFileList) {
	if(-r "$bamFile.bai") {
	} else {
		(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
		if(-r $baiFile) {
			system("ln -s $baiFile $bamFile.bai");
		} else {
			system("samtools index $bamFile");
		}
	}
}
my $pid = open2(my $reader, my $writer, 'sort -n');
forkPrintParentWriter($writer);
my @columnList = ();
{
	open(my $reader, $metaPrismFile);
	chomp(my $line = <$reader>);
	@columnList = split(/\t/, $line, -1);
	my $order = 0;
	my @tokenListList = ();
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		push(@tokenListList, [$order += 1, @tokenList]);
		if(scalar(@tokenListList) >= $numberPerThread) {
			if($threads == 1) {
				printAbundance(@tokenListList);
			} else {
				forkPrintSubroutine(\&printAbundance, @tokenListList);
			}
			@tokenListList = ();
		}
	}
	if(@tokenListList) {
		if($threads == 1) {
			printAbundance(@tokenListList);
		} else {
			forkPrintSubroutine(\&printAbundance, @tokenListList);
		}
	}
	close($reader);
	forkPrintWait();
}
forkPrintParentWriter();
close($writer);
{
	my @tokenHashList = ();
	my %geneAbundanceHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(my $order, @tokenHash{@columnList, 'abundance'}) = split(/\t/, $line, -1);
		push(@tokenHashList, \%tokenHash);
		$geneAbundanceHash{$tokenHash{'gene'}} += $tokenHash{'abundance'};
	}
	if($baseAbundanceGenes ne '') {
		my @abundanceList = map {defined($_) ? $_ : 0} @geneAbundanceHash{split(/,/, $baseAbundanceGenes)};
		if($ignoreOneSideOutlierNumberOfBaseAbundanceGenes > 0) {
			@abundanceList = sort {$a <=> $b} @abundanceList;
			@abundanceList = @abundanceList[$ignoreOneSideOutlierNumberOfBaseAbundanceGenes .. $#abundanceList - $ignoreOneSideOutlierNumberOfBaseAbundanceGenes];
		}
		$baseAbundance = sum(@abundanceList) / scalar(@abundanceList);
	}
	print join("\t", @columnList, 'abundance'), "\n";
	if($baseAbundance eq '') {
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			print join("\t", @tokenHash{@columnList, 'abundance'}), "\n";
		}
	} elsif($baseAbundance > 0) {
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			$tokenHash{'abundance'} = $tokenHash{'abundance'} / $baseAbundance;
			print join("\t", @tokenHash{@columnList, 'abundance'}), "\n";
		}
	} else {
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			$tokenHash{'abundance'} = "$tokenHash{'abundance'}/$baseAbundance";
			print join("\t", @tokenHash{@columnList, 'abundance'}), "\n";
		}
	}
}
close($reader);
waitpid($pid, 0);
if($use_merged_bam_file || $use_samtools_depth) {
	system("rm $temporaryPrefix.bam");
	if($use_samtools_depth) {
		system("rm $temporaryPrefix.filtered.bam");
		if($stranded) {
			system("rm $temporaryPrefix.filtered.f.bam");
			system("rm $temporaryPrefix.filtered.f.bam.bai");
			system("rm $temporaryPrefix.filtered.r.bam");
			system("rm $temporaryPrefix.filtered.r.bam.bai");
		} else {
			system("rm $temporaryPrefix.filtered.bam.bai");
		}
	} else {
		system("rm $temporaryPrefix.bam.bai");
	}
}

sub printAbundance {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} @bamFileList if($use_sam_module);
	foreach(@_) {
		my ($order, @tokenList) = @$_;
		my %tokenHash = ();
		@tokenHash{@columnList} = @tokenList;
		my ($chromosome, $start, $end, $strand) = @tokenHash{'chromosome', 'start', 'end', 'strand'};
		my $depth = 0;
		if($use_samtools_depth) {
			if($stranded eq '') {
				$depth = getDepth("$temporaryPrefix.filtered.bam", $chromosome, $start, $end);
			} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
				$depth = getDepth("$temporaryPrefix.filtered.f.bam", $chromosome, $start, $end);
			} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
				$depth = getDepth("$temporaryPrefix.filtered.r.bam", $chromosome, $start, $end);
			}
		} elsif($use_sam_module) {
			foreach my $sam (@samList) {
				foreach my $alignment ($sam->get_features_by_location(-seq_id => $chromosome, -start => $start, -end => $end)) {
					next if($alignment->qual < $minimumMappingQuality);
					next if(($alignment->flag & int($includeFlag)) != $includeFlag);
					next if($alignment->flag & int($excludeFlag));
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($alignment->flag & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($alignment->flag & 253) == $_} 81, 161, 16);
					}
					my @positionList = getPositionList($alignment->start, $alignment->cigar_str);
					@positionList = grep {$_ ne ''} @positionList;
					@positionList = grep {$start <= $_ && $_ <= $end} @positionList;
					$depth += scalar(@positionList);
				}
			}
		} else {
			foreach my $bamFile (@bamFileList) {
				my $reader;
				if($stranded eq '') {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | samtools depth -d 0 - |");
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@tokenList = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$tokenList[1] & 253) == \$_} 97, 145, 0);' | samtools depth -d 0 - |");
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@tokenList = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$tokenList[1] & 253) == \$_} 81, 161, 16);' | samtools depth -d 0 - |");
				}
				while(my $line = <$reader>) {
					chomp($line);
					my ($chromosome, $position, $depth1) = split(/\t/, $line, -1);
					$depth += $depth1 if($start <= $position && $position <= $end);
				}
				close($reader);
			}
		}
		my $length = $end - $start + 1;
		my $meanDepth = $depth / $length;
		my $abundance = $meanDepth;
		forkPrint(join("\t", $order, @tokenList, $abundance), "\n");
	}
}

sub getDepth {
	my ($bamFile, $chromosome, $start, $end) = @_;
	my $depth = 0;
	open(my $reader, "samtools depth -a -a -r $chromosome:$start-$end $bamFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'chromosome', 'position', 'depth'} = split(/\t/, $line, -1);
		$depth += $tokenHash{'depth'};
	}
	close($reader);
	return $depth;
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}
