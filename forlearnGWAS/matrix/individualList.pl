#!perl -w
use strict;
my $lineNumber=0;
open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
	if( $line=~/^(\S+)\s/ ){
		print "lineNumber\t$lineNumber\t$1\n";
	}
	$lineNumber++;
}
close INPUT;