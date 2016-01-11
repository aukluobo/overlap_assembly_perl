use strict;use warnings;

#huangguodong@genomics.cn 20141014

die "perl $0 <s2.o> <oo> <raw_map1> <raw_map2>" if(@ARGV==0);

my ($result,$oo,$map1,$map2)=@ARGV;

my %map1;my %map2;my %map1_ctg_sort;my %map2_ctg_sort;
read_map($map1,\%map1,\%map1_ctg_sort);
read_map($map2,\%map2,\%map2_ctg_sort);

my %block;
open IN,$oo or die;
while(<IN>){
	if(/^blo/ && /CTG/){
		my @aa=split;
		my $block=shift @aa;
		foreach my $ctg(@aa){
			$block{$block}=[@aa];
		}
	}
}
close IN;

open IN,$result or die;
while(<IN>){
	if(/block/){
		my @aa=split;
		my $block=$aa[0];
		foreach my $ctg(@{$block{$block}}){
			print "$map1{$ctg}[1]\t$map2{$ctg}[1]\t$block\n";
		}
	}else{print;}
}
close IN;




sub read_map{
	my ($map,$hash,$hash1)=@_;
	my $num=1;
	open IN,$map or die;
	while(<IN>){
		my @aa=split;chomp;
		$hash->{$aa[0]}=[$num,$_];
		$hash1->{$num}=$aa[0];
		$num++;
	}
	close IN;
}


