#!perl -w
use strict;

my $similarity = $ARGV[0];

my %individuals;
my $id=0;
open INPUT, "snp.ped";
while( my $line=<INPUT> ){
	if( $line=~/^(\S+)/ ){
		$individuals{$id}=$1;
		$id++;
	}
}
close INPUT;

my %unliableNumber;
open INPUT, "../uncertain";
while( my $line=<INPUT> ){
	if( $line=~/(\S+)\s+(\S+)/ ){
		$unliableNumber{$2}=$1;
	}
}
close INPUT;

my $lineNumber=0;
open INPUT, "plink.mibs";
while( my $line=<INPUT> ){
	my @cells=split(" ", $line);
	my $lineSize = @cells;

	my $i=0;
	while( $i < $lineSize ){
		if( $cells[$i]>$similarity && $lineNumber<$i ){
			if( $unliableNumber{$individuals{$lineNumber}} > $unliableNumber{$individuals{$i}} ){
				print "$individuals{$lineNumber}\t$individuals{$lineNumber}\n";
			}else{
				print "$individuals{$i}\t$individuals{$i}\n";
			}
		}
		$i++;
	}
	$lineNumber++;
}
close INPUT;
