#!/usr/bin/perl -w
# PSR_poly_finder.pl
# AUTHOR: Concita Cantarella, Nunzio D'Agostino
# LAST REVISED: September 2015
				

use strict;
use POSIX;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use DBI;
use DBD::mysql;
use File::Copy;
use DateTime;

my $tabSep="\t";
my $newline="\n";

my $usage = $newline.$newline."usage: $0 [-i <input folder> -o <output directory >]".$newline.
            "Reads MISA output".$newline.
            "-i dir 	directory containing  PSR_counts.txt files, one for each sample".$newline.$newline.
            "-o dir	output dir".$newline;

my %opt; 
getopts("i:o:",\%opt);
if (!defined $opt{'i'})  {die $usage}; # mandatory arguments

my $date = DateTime->now;
my $ymd    = $date->ymd;
my $hms    = $date->hms('_');
my $formatDate=$ymd.'_'.$hms;

my $mysqldir='/var/lib/mysql/PSR_DB/';
my $outdir='/tmp/'; 


my $PSRcountsFolder=$opt{'i'};
if (defined $opt{'o'}) {$outdir=$opt{'o'};}

my $logFilename='PSR_poly_finder.log';
my $outfilename='PSR_polySSR_'.$formatDate.'.txt';
my $dbname='PSR_DB';
my @tablesname;
my @fieldsname=('Seq_ID','SSR','start','stop','Num_of_repeatUnit','Num_of_reads','ID');

my $logfilepath=$outdir.$logFilename;
my $outfilepath=$outdir.$outfilename;

my $host = 'localhost';       # MySQL Host
my $user = 'psr';             # your username
my $password ='';
my $mysqlfilename=$mysqldir.$outfilename;
my $unionTable='PSR_allPolymorphysm';
my $compressSuffix='_compress';
my $polyPrefix='poly_';

open(my $lfh, '>', $logfilepath) or die "Could not open file $logfilepath $!";
print $newline."log file = $logfilepath".$newline;
open(my $outfh, '>', $outfilepath) or die "Could not open file $outfilepath $!";

# connect to PSR_DB
my $dbh = DBI->connect("dbi:mysql:database=$dbname;host=$host", $user, $password)
    or die "Couldn't connect to database: $DBI::errstr".$newline;


$dbh->do("USE $dbname;");

my @suffixlist=('.txt'); 
opendir(my $indir, $PSRcountsFolder) or die $!;

while (my $countfile = readdir($indir)){
	next if ($countfile =~ m/^\./);
	print $lfh $newline.$countfile.$newline;
	my $tablename = basename($countfile,@suffixlist);
	push @tablesname, $tablename;
	createTable($tablename);
	loadData(($tablename,$countfile));
	compressTable($tablename);
}

createUnionTable(@tablesname);
alterUnionTable(@tablesname);
updateUnionTable(@tablesname);

getResult(@tablesname);

move($mysqlfilename,$outfilepath);
print $newline."output file $outfilepath written!".$newline.$newline;
dropTables(@tablesname);


close ($lfh);
$dbh->disconnect();
closedir($indir);


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#			SUBROUTINES
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sub createTable{
	my $tablename=$_[0];
	my $query="CREATE TABLE $dbname.$tablename ($fieldsname[0] varchar(150),$fieldsname[1] varchar(10),$fieldsname[2] int(12),$fieldsname[3] int(12),$fieldsname[4] int(7),$fieldsname[5] int(7),$fieldsname[6] int(7) auto_increment,PRIMARY KEY (ID)); ";
	print $lfh $newline.$query.$newline;
	my $sth = $dbh->prepare($query);
	$sth->execute();
}


sub loadData{
	my $tablename=$_[0];
	my $countfile=$_[1];
	my $query= "LOAD DATA LOCAL INFILE '$PSRcountsFolder$countfile' INTO TABLE $dbname.$tablename IGNORE 1 LINES;";
	print $lfh $newline.$query.$newline;
	my $sth = $dbh->prepare($query);
	$sth->execute();
}


sub compressTable{
	my $sampleTable=$_[0];
	my $query='';
	$query.="CREATE TABLE $dbname.$sampleTable$compressSuffix ";
	$query.="SELECT $dbname.$sampleTable.$fieldsname[0],$dbname.$sampleTable.$fieldsname[1],$dbname.$sampleTable.$fieldsname[2],$dbname.$sampleTable.$fieldsname[3], ";
	$query.="GROUP_CONCAT(DISTINCT $dbname.$sampleTable.$fieldsname[4] ORDER BY $dbname.$sampleTable.$fieldsname[4] DESC SEPARATOR ';') ";
	$query.="AS $polyPrefix$sampleTable FROM $dbname.$sampleTable GROUP BY $dbname.$sampleTable.$fieldsname[0];";
	print $lfh $newline.$query.$newline;	
	my $sth = $dbh->prepare($query);
	$sth->execute();
}


sub createUnionTable{
	my $query='';
	$query.="CREATE TABLE $dbname.$unionTable AS SELECT ";
	my $i=1;
	foreach my $sampleTable (@_) {
		$query.="$dbname.$sampleTable$compressSuffix.$fieldsname[0] FROM $dbname.$sampleTable$compressSuffix ";
		if ($i < $#_+1){ $query.=" UNION SELECT ";}
		$i ++;	
	}
	$query.="GROUP BY $fieldsname[0];";
	print $lfh $newline.$query.$newline;
	my $sth = $dbh->prepare($query);
	$sth->execute();
}


sub alterUnionTable{	
	my $query='';
	$query.="ALTER TABLE $dbname.$unionTable ADD COLUMN $fieldsname[1] varchar(10) NOT NULL DEFAULT 0 AFTER $fieldsname[0], ADD COLUMN $fieldsname[2] int(12) NOT NULL DEFAULT 0 AFTER $fieldsname[1], ADD COLUMN $fieldsname[3] int(12) NOT NULL DEFAULT 0 AFTER $fieldsname[2]";
	my $j=0;
	foreach my $sampleTable (@_) {
		$query.=", ADD COLUMN $polyPrefix$sampleTable varchar(30) NOT NULL DEFAULT 0 AFTER";
		if ($j==0) {$query.=	" $fieldsname[3] ";}
		elsif ($j > 0) {
			my $previousField=$_[$j-1];
			$query.= " $polyPrefix$previousField ";
		}
		$j ++;			
	}
	$query.= ';';
	print $lfh $newline.$query.$newline;
	my $sth = $dbh->prepare($query);
	$sth->execute();
}	


sub updateUnionTable{
	foreach my $sampleTable (@_) {
		my $query='';
		$query.="UPDATE $dbname.$unionTable INNER JOIN $dbname.$sampleTable$compressSuffix ON $dbname.$unionTable.$fieldsname[0] = $dbname.$sampleTable$compressSuffix.$fieldsname[0]";
		$query.=" SET ";
		$query.="$dbname.$unionTable.$fieldsname[1] = $dbname.$sampleTable$compressSuffix.$fieldsname[1],";
		$query.="$dbname.$unionTable.$fieldsname[2] = $dbname.$sampleTable$compressSuffix.$fieldsname[2],";
		$query.="$dbname.$unionTable.$fieldsname[3] = $dbname.$sampleTable$compressSuffix.$fieldsname[3],";
		$query.="$dbname.$unionTable.$polyPrefix$sampleTable = $dbname.$sampleTable$compressSuffix.$polyPrefix$sampleTable;";
		print $lfh $newline.$newline.$query.$newline;	
		my $sth = $dbh->prepare($query);
		$sth->execute();
	}
}	


sub getResult{
	my $query='';
	my $i=0;
	$query.="SELECT '$fieldsname[0]','$fieldsname[1]','$fieldsname[2]', '$fieldsname[3]',";
	foreach my $sampleTable (@_) {
		$query.="'$polyPrefix$sampleTable'";
		if ($i < $#_) {$query.=',';}
		$i++;
	}
	$query.=" UNION ALL ";
	$query.="SELECT * from $dbname.$unionTable INTO OUTFILE '$mysqlfilename'";
	print $lfh $newline.$query.$newline;
	my $sth = $dbh->prepare($query);
	$sth->execute();
}


sub dropTables{
	my $query='';
	$query.="DROP TABLE";
	foreach my $table (@_) {
		$query.=" $dbname.$table,";
		$query.=" $dbname.$table$compressSuffix,";
	}
	$query.=" $dbname.$unionTable;";
	print $lfh $newline.$query.$newline;	
	my $sth = $dbh->prepare($query);
	$sth->execute();
}
