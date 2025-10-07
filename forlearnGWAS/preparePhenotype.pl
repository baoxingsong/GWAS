#！perl -w
use strict;

my %pheno;
open INPUT, "yiqiaoPheno.txt";
while( my $line=<INPUT> ){
	if ($line=~/^(\S+)\s(.*)$/){
		$pheno{$1} = $2;
	}
}
close INPUT;

open INPUT, "hwc.snp.chr.mis9_maf5.ibs.fam.backup";
while( my $line=<INPUT> ){
	my @elements = split('\s', $line);
	if (exists $pheno{$elements[0]} ){
		print "$elements[0]\t$elements[0]\t0\t0\t0\t$pheno{$elements[0]}\n";
	}else{
		print "$elements[0]\t$elements[0]\t0\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
}
close INPUT;
