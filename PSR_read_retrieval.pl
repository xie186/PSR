#!/usr/bin/perl -w
# PSR_read_retrieval.pl
# AUTHOR: Concita Cantarella, Nunzio D'Agostino
# LAST REVISED: September 2015


use strict;
use POSIX;
use Getopt::Std;
use Text::CSV;
use Data::Dumper;
use Bio::DB::Sam;
use File::Basename;


my $tabSep="\t";
my $newline="\n";

my $usage = "$newline usage: $0 [-m <tabdel file> -d <directory with BAM > -o <output directory> [-t minimum read depth] [-p polymorphism percentage depth] $newline".
            "-m tab-del		tab-delimited file with detected SSR  $newline".
            "-d dir		directory path containing bam files $newline".
            "-o 		output dir $newline".
            "-t 		Minimum read depth at a position to make a call (default =10) $newline".
            "-p 		Minimum supporting reads at a position to call variants (default=30%) $newline$newline";

my %opt; 
getopts("m:d:o:t:p:",\%opt) ; 

if ((!defined $opt{'m'}) || (!defined  $opt{'d'}) || (!defined  $opt{'o'}))  {print $usage; exit;} # mandatory arguments

my $fileMisa=$opt{'m'};
my $directory=$opt{'d'}; 
my $outdir=$opt{'o'};
my $logFilename='PSR_reads_retrieval.log';
my $tabFilename='PSR_counts.txt';

my $threshold=10; # min number of reads; 
my $polyThreshold=30; #polymorphism percentage 

if (defined $opt{'t'}) {$threshold=$opt{'t'};}
if (defined $opt{'p'}) {$polyThreshold=$opt{'p'};}

print $newline.'Minimum read depth ='.$threshold.$newline;
print $newline.'Minimum percentage to call variants ='.$polyThreshold.$newline;

print $newline.$newline.'wait.... be patient!...'.$newline.$newline;

my $outfileExt='.fa';
my $logFile=$outdir.$logFilename;
my $tabFile=$outdir.$tabFilename;

open(my $lfh, '>', $logFile) or die "Could not open file $logFile $!";
open(my $tabfh, '>', $tabFile) or die "Could not open file $tabFile $!";

my @misaRows;
my $csv = Text::CSV->new ( { sep_char => $tabSep } )  
	or die "Cannot use CSV: ".Text::CSV->error_diag ();

open(my $ifh, "<:encoding(UTF-8)", $fileMisa) 
	or die "can't open UTF-8 encoded filename: $!"; 

#header file of detected SSRs - it must contain 7 columns 
my @header=($csv->getline( $ifh ));
my $ssrID_label	=$header[0][0];  
my $ssrNR_label	=$header[0][1];  
my $ssrType_label 	=$header[0][2];
my $ssr_label		=$header[0][3];
my $ssrSize_label	=$header[0][4];
my $ssrStart_label	=$header[0][5];
my $ssrEnd_label	=$header[0][6];

my $typePrefix='p'; #accepted microsatellites type

$csv->column_names(@header);

while ( my $row = $csv->getline_hr( $ifh ) ) {
	chomp $row;
	push @misaRows, $row;
}

my $misasize=$#misaRows+1;
print $lfh "SSR total = $misasize $newline";
print $tabfh 'Seq_ID'.$tabSep.'SSR'.$tabSep.'start'.$tabSep.'stop'.$tabSep.'# repeat unit'.$tabSep.'# reads'.$newline;

$csv->eof or $csv->error_diag();
close $ifh;
$csv -> DESTROY;

my %unitRepCounts;
my %reads;

	
for (my $i = 0; $i < $misasize; ++$i){
	my $transcrp_id	=  $misaRows[$i]{$ssrID_label};
	my $ssrNR	= $misaRows[$i]{$ssrNR_label};
	my $SSRtype	=$misaRows[$i]{$ssrType_label};
	# select only perfect SSR	
	if (substr($SSRtype,0,1) ne $typePrefix) {print $lfh $tabSep."$SSRtype isn't a perfect microsatellites $newline"; next;}
	elsif (substr($SSRtype,0,1) eq $typePrefix){ 
		my $ssrUnitLength	= substr($SSRtype,length($typePrefix));
		my $SSR		= $misaRows[$i]{$ssr_label};
		my $ssr		= substr($SSR,1,$ssrUnitLength);
		my $ssr_size	= $misaRows[$i]{$ssrSize_label};
		my $ssr_start	= $misaRows[$i]{$ssrStart_label};
		my $ssr_stop	= $misaRows[$i]{$ssrEnd_label};
		my $unitReps	= $ssr_size/$ssrUnitLength;

		my $minUnitRep =  4;
		
		# misa.ini
		# definition(unit_size,min_repeats):  1-10 2-6 3-5 4-5 5-5 6-5
		
		SWITCH:{
			$ssrUnitLength ==1  && do {$minUnitRep =  8;  last SWITCH; };
			$ssrUnitLength ==2  && do {$minUnitRep =  4;  last SWITCH; };
			$ssrUnitLength ==3  && do {$minUnitRep =  3;  last SWITCH; };
			$ssrUnitLength ==4  && do {$minUnitRep =  3;  last SWITCH; };
			$ssrUnitLength ==5  && do {$minUnitRep =  3;  last SWITCH; };
			$ssrUnitLength ==6  && do {$minUnitRep =  3;  last SWITCH; };
		}		

		my $summaryhead=$transcrp_id.$tabSep.$SSR.$tabSep.$ssr_start.$tabSep.$ssr_stop;
		print $lfh $summaryhead.$newline;		
		print $lfh "$tabSep Input files: $newline";

		opendir(my $dir, $directory) or die $!;
		my @bams = grep { /\.bam$/ && -f "$directory/$_" } readdir($dir);
		closedir($dir);
		
		#retrieve reads aligning the current SSR   
		foreach my $bamfile (@bams) {
			my $sam = Bio::DB::Sam->new(-bam => $directory.$bamfile);
			my @alignments 	= $sam->get_features_by_location(-seq_id => $transcrp_id,-start => $ssr_start,-end => $ssr_stop);
														 
			my $alignsize=$#alignments +1;
			print $lfh "$tabSep$tabSep $bamfile $newline";
			print $lfh "$tabSep$tabSep num of reads aligning SSR = $alignsize $newline";

			my @todelete; #reads aligned >1 times
				
			for my $a (@alignments) {
				my $refId			= $a->seq_id;
				my $start			= $a->start;
				my $end			= $a->end;
				my $strand		= $a->strand;
				my $cigar		= $a->cigar_str;
				my $queryName	= $a->qname;
				my $query_seq	= $a->query->dna;
				my $maxUnitRep  = int(length($query_seq)/$ssrUnitLength);
				my @alignInfo;
				if ( exists $reads{$queryName}) {
					push (@todelete,$queryName); 
					print $lfh $tabSep.$tabSep.$tabSep.$queryName.' was discarded due to multiple alignment in the same transcript/genomic locus'.$newline;
					next;
				}
				else {
					if ($start < $ssr_start and $end > $ssr_stop ){
						my $ssrSeqInRead;
						my $i=0;
						while ($query_seq =~ m/(($ssr){$minUnitRep,$maxUnitRep})/g){
							$ssrSeqInRead = $1; 
							my $ssrSeqStart=index($query_seq,$ssrSeqInRead,$i);
							my $ssrSeqStop=$ssrSeqStart+length($ssrSeqInRead);
							if (($ssrSeqStart > 0) and ($ssrSeqStop < length($query_seq))) {
								my $length5p=$ssrSeqStart;
								my $overhang5p=substr($query_seq,0,$length5p);
								my $length3p=length($query_seq)-$ssrSeqStop;
								my $overhang3p=substr($query_seq,-$length3p);
								my $ssrUnitNum;
								if (length($overhang5p) < length ($ssr)) {
									if ($overhang5p eq substr($ssr,-length($overhang5p))){
										print $lfh $tabSep.$tabSep.$tabSep.$queryName." was discarded due to undefined 5' limit".$newline;
										$i=$ssrSeqStop;
										next;
									}	
									elsif (length($overhang3p) < length ($ssr)){
										if ($overhang3p eq substr($ssr,0,length($overhang3p))) {
											print $lfh $tabSep.$tabSep.$tabSep.$queryName." was discarded due to undefined 3' limit".$newline;
											$i=$ssrSeqStop;
											next;
										}
										else {
											$ssrUnitNum = int(length($ssrSeqInRead)/$ssrUnitLength);
											@alignInfo=($query_seq,$ssrUnitNum);
											$reads{$queryName}=[@alignInfo];
										}
									}
								}
								else {
										$ssrUnitNum = int(length($ssrSeqInRead)/$ssrUnitLength);
										@alignInfo=($query_seq,$ssrUnitNum);
										$reads{$queryName}=[@alignInfo];
								}
							}
							else { print $lfh $tabSep.$tabSep.$tabSep.$queryName.' was discarded because the SSR is to one of two extremes '.$query_seq.$newline;}
							$i=$ssrSeqStop;
						}
					}	
					else {
						print $lfh $tabSep.$tabSep.$tabSep.$queryName.' was discarded due to its start e stop positions '.$start.$tabSep.$end.$newline;
					}
				}
			}

			foreach my $repetitiveRead (@todelete) {
				delete $reads{$repetitiveRead};
			}
			undef @alignments;
			undef @todelete;
			undef $sam;
		}

		my @readCount=keys(%reads);
		my $size=$#readCount +1;
		my $sumOfReads =$size;
		my %stats; # hash array (ssrUnitReps -> number of reads)
		my %polyReads;
		print $lfh "$tabSep total useful reads = ".$sumOfReads.$newline;
		print $lfh "$tabSep Minimum read depth = ".$threshold.$newline;
		if ($sumOfReads >= $threshold ){
			foreach my $readID (keys %reads) {
				my $readSequence= $reads{$readID}[0];
				my $unitNum = $reads{$readID}[1];
				push @{ $polyReads{$unitNum} }, '>'.$readID.$newline.$readSequence;  
				if (not exists $stats{$unitNum}) {$stats{$unitNum}=1;}
				else {$stats{$unitNum}++;}
			}

			my @sortedkeys;
			my @sortedvals;
			# sort numerically descending
			@sortedkeys = sort { $stats{$b} <=> $stats{$a} } keys(%stats); 
			@sortedvals = @stats{@sortedkeys};
		
			my $n=0;
			my $polymorfism;

			foreach $polymorfism (@sortedkeys){
				my $polyreadsNum=$sortedvals[$n];
				my $result;
				print $lfh "$tabSep Number of reads with $polymorfism repeated unit = $polyreadsNum $newline";
				if ($polyreadsNum >= ($polyThreshold *$sumOfReads)/100 and $polyreadsNum >=$threshold ){
					print $tabfh $transcrp_id.$tabSep.$SSR.$tabSep.$ssr_start.$tabSep.$ssr_stop.$tabSep.$polymorfism.$tabSep.$polyreadsNum.$newline;
					my $outfilename='PSR_'.$transcrp_id.'_'.$ssr.'_'.$ssr_start."-".$ssr_stop.'_'.$polymorfism.'rU'.$outfileExt;											 
					$outfilename =~ s/\//_/g;
					my $outfile=$outdir.$outfilename;
					open(my $ofh, '>', $outfile) or die "Could not open file $outfile $! $newline";
					my @fileContent=@{$polyReads{$polymorfism}};
					$result=join($newline,@fileContent);
					print $ofh $result; 
					close $ofh;
					print $lfh "$tabSep $outfile written! $newline";
					undef $result;
					undef $outfile;
				}
				$n++;
			}
			undef @sortedkeys;
			undef @sortedvals;
			undef $polymorfism;
		}
		undef %reads;
		undef $sumOfReads;
		undef @bams;
	}
}
close $tabfh;
close $lfh;	
print "$logFile written! $newline";
print "$tabFile written! $newline";
