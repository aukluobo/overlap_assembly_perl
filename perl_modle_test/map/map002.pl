use strict;
use warnings;

#huangguodong@genomics.cn 20140924

die "perl $0 <out_form.m8.ff.changbactoctg.cM.83.freq.map>" if(@ARGV==0);

my ($file)=@ARGV;

my %all;
open IN,$file or die;
while(<IN>){
	my @aa=split;
	$all{$aa[1]}{$aa[0]}=[@aa];
}
close IN;

foreach my $cm(sort {$a<=>$b} keys %all){
	my @first;my @last;my @midd;
	foreach my $ctg(keys $all{$cm}){
		my @bb=@{$all{$cm}{$ctg}};
		#three region
		if(@bb==7){
			push @last,[@bb];
		}elsif(@bb==5){
			#two region
			push @midd,[@bb];
		}else{
			#one region
			push @first,[@bb];
		}
	}
	if(@first>0){
		foreach my $key(@first){
			my @cc=@{$key};
			foreach my $key_c(@cc){
				print "$key_c\t";
			}
			print "\n";
		}
	}
	if(@midd){
		foreach my $key(@midd){
			my @cc=@{$key};
			foreach my $key_c(@cc){
				print "$key_c\t";
			}
			print "\n";
		}
	}
	if(@last){
		foreach my $key(@last){
			my @cc=@{$key};
			foreach my $key_c(@cc){
				print "$key_c\t";
			}
			print "\n";
		}
	}
}




