use strict;use warnings;

#huangguodong@genomics.cn 20141014

die "perl $0 <s2.o> <oo> <raw_map1> <raw_map2> <map_sort:1/-1>" if(@ARGV==0);

my ($result,$oo,$map1,$map2,$map_sort)=@ARGV;

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
my %check;my $num=1;my %all;my $check_num=1;my %cm_type;
while(<IN>){
	my @aa=split;
	$all{$num}=[@aa];
	$num++;
	my @bb;if($aa[0]=~/CTG/){@bb=grep /CTG/,@aa;}else{@bb=grep /block/,@aa;}
	if(@bb>1){
		$check{$check_num}=[@aa];$check_num++;$cm_type{$aa[1]}=1;
	}
}
close IN;
my @check=sort {$a<=>$b} keys %check;

#undefined code
my %point;
for(my $i=0;$i<@check-1;$i++){
	if($check{$check[$i+1]}[1]>$check{$check[$i]}[1]){
		$point{$i+1}[0]=$check{$check[$i+1]}[1];
		my $cut=$check{$check[$i+1]}[1];
		my %tmp;
		for(my $j=0;$j<=$i;$j++){
			if($cut>$check{$check[$j]}[1]){$tmp{$check{$check[$j]}[1]}=$j;}
		}
		my $lest=(sort {$a<=>$b} keys %tmp)[0];
		$point{$i+1}[1]=$tmp{$lest};
		%tmp=();
		for(my $j=$i+1;$j<@check;$j++){
			if($lest>=$check{$check[$j]}[1]){
				$point{$i+1}[2]=$j;
				$point{$i+1}[3]=0;
				if($lest>=$check{$check[$j]}[1]){$point{$i+1}[3]=1;}
				last;
			}
		}
	}
}
#undefined code

#get cm sort
my @cm;
if($map_sort>0){@cm=sort {$a<=>$b} keys %cm_type;}
else{@cm=sort {$b<=>$a} keys %cm_type;}
my $j=0;my $pre;
for(my $i=0;$i<@check;$i++){
	unless($pre){$pre=$check{$check[$i]}[1];}
	if($pre != $check{$check[$i]}[1]){$check{$check[$i]}[1]=$cm[$j];$j++;}
	else{$check{$check[$i]}[1]=$cm[$j];}
	if($check{$check[$i]}[0]=~/block/){
		my $block=$check{$check[$i]}[0];
		foreach my $ctg(@{$block{$block}}){
			my $part1=$map1{$ctg}[1];
			my @part1=split /\s+/,$part1;
			$part1[1]=$check{$check[$i]}[1];
			my $out1=join "\t",@part1;
			print "$out1\t$map2{$ctg}[1]\t$block\n";
		}
	}else{
		my $out=join "\t",@{$check{$check[$i]}};
		print "$out\n";
	}
	$pre=$check{$check[$i]}[1];
}



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


