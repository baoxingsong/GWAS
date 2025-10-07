#!perl -w
use strict;

my %individuals;
my $id=0;
open INPUT, "snp_matrix.fam";
while( my $line=<INPUT> ){
	if( $line=~/^(\S+)/ ){
		$individuals{$id}=$1;
		$id++;
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
		if( $cells[$i]>0.95 && $lineNumber<$i ){
			print "$individuals{$lineNumber}\t$individuals{$lineNumber}\n";
		}
		$i++;
	}
	$lineNumber++;
}
close INPUT;
