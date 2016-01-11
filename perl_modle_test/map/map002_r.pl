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

foreach my $cm(sort {$b<=>$a} keys %all){
	my @first;my @last;my @midd;
	foreach my $ctg(keys $all{$cm}){
		my @cc=@{$all{$cm}{$ctg}};
		my $cc_ctg=shift @cc;
		my @bb;my %tmp;
		push @bb,$cc_ctg;
		for(my $i=0;$i<@cc;$i+=2){
			$tmp{$cc[$i]}=$cc[$i+1];
		}
		foreach my $cc_cm(sort {$b<=>$a} keys %tmp){
			push @bb,$cc_cm,$tmp{$cc_cm};
		}
		#three region
		if(@bb==7){
			push @first,[@bb];
		}elsif(@bb==5){
			#two region
			push @midd,[@bb];
		}else{
			#one region
			push @last,[@bb];
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




