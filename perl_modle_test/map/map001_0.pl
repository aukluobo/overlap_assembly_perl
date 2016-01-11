use strict;
use warnings;

#huangguodong@genomics.cn 20141104

# this program is try to get the best record for the marker

die "perl $0 <lst.cm>" if(@ARGV==0);

my ($lst)=@ARGV;

open IN,$lst or die;
my %all;
while(<IN>){
	my @aa=split;my $rec=$_;chomp $rec;
	#4.4     3AS_3432506     CTG1579 15      951     9       936     ---     100
	my $len=$aa[4]-$aa[3]+1;
	push @{$all{$aa[1]}{$len}},$rec;
}
close IN;

foreach my $mar(keys %all){
	my @key=sort {$a<=>$b} keys $all{$mar};
	my $longest=$key[-1];
	my @key1=@{$all{$mar}{$longest}};
	my %tmp;my $smal=100000000000;
	if(@key1>1){
		foreach my $rec(@key1){
			my @aa=split /\s+/,$rec;
			if($aa[-1]<$smal){$smal=$aa[-1];}
		}
		my @tmp;
		foreach my $rec(@key1){
			my @aa=split /\s+/,$rec;
			if($aa[-1]==$smal){push @tmp,$rec;}
			elsif($aa[-1]<=10){push @tmp,$rec;}
		}
		my $block=10000000;
		if(@tmp>1){
			#check block
			foreach my $rec(@tmp){
				my @aa=split /\s+/,$rec;
				my @bb=split //,$aa[-2];
				if(@bb<$block){$block=@bb;}
			}
			my @last;
			foreach my $rec(@tmp){
				my @aa=split /\s+/,$rec;
				my @bb=split //,$aa[-2];
				if(@bb==$block){
					#print $rec,"\n";
					push @last,$rec;
				}
			}
			if(@last<=2){
				foreach my $re(@last){print "$re\n";}
			}
		}else{
			#only one best record
			foreach my $rec(@tmp){
				my @aa=split /\s+/,$rec;
				print $rec,"\n";
			}
		}
	}else{
		print "$key1[0]\n";

	}
}



