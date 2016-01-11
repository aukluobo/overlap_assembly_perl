use strict;use warnings;

die "perl $0 <relation.lst>" if(@ARGV==0);

my $bac_num=100;

my $nn=0;my $cal=0;
open OUT,">test.$nn" or die;
open IN,$ARGV[0] or die;
while(<IN>){
	if(/bac/){
		print OUT $_;
		$cal++;
	}else{
		if($cal >=$bac_num){
			$cal=0;$nn++;
			close OUT;
			open OUT,">test.$nn" or die;
			print OUT $_;
		}else{
			print OUT $_;
		}
	}
}
close IN;
close OUT;




