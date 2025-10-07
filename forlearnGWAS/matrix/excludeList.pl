#!perl -w
use strict;

my %avaliableNumber;
open INPUT, "plink.lmiss";
my $linenumber=0;
while( my $line=<INPUT> ){
	if( $line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ && $linenumber >0){
		my $number=$4-$3;
		$avaliableNumber{$2}=$number;
	}
	$linenumber++;
}
close INPUT;

$linenumber=0;
open OUTPUT, ">excludeSnps";
open INPUT, "plink.frq";
while( my $line=<INPUT> ){
	if( $line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ && $linenumber>0){
		if( $5 eq "NA" ){
			print OUTPUT "$2\n";
		}elsif( $avaliableNumber{$2}*$5 <5 ){
			print OUTPUT "$2\n";
		}elsif( $1 eq '0' ){
			print OUTPUT "$2\n";
		}elsif( $1 eq '26'){
			print OUTPUT "$2\n";
		}
	}
	$linenumber++;
}
close INPUT;
close OUTPUT;

