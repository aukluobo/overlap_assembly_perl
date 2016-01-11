use strict;use warnings;

#huangguodong@genomics.cn 20141014

die "perl $0 <s2.o> <oo> <raw_map1> <raw_map2> <map_sort:1/-1>" if(@ARGV==0);

my ($result,$oo,$map1,$map2,$map_sort)=@ARGV;

my %map1;my %map2;my %map1_ctg_sort;my %map2_ctg_sort;my %cm1;my %cm2;
read_map($map1,\%map1,\%map1_ctg_sort,\%cm1);
read_map($map2,\%map2,\%map2_ctg_sort,\%cm2);

my %cm1_sort;my $c_n=1;
my @key_cm1;
if($map_sort>0){@key_cm1=sort {$a<=>$b} keys %cm1;}else{@key_cm1=sort {$a<=>$b} keys %cm1;}
foreach(@key_cm1){$cm1_sort{$c_n}=$_;$cm1{$_}=$c_n;$c_n++;}
print STDERR "$c_n\n";

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
my %all;my $num=1;my %cm_type;my %ex;my %at;
my $pre;my @cm_type;
while(<IN>){
	my @aa=split;
	$all{$num}=[@aa];
	if(($map1{$aa[0]}[1] and $aa[1]==(split /\s+/,$map1{$aa[0]}[1])[1]) or $aa[0]=~/block/){
		push @cm_type,$aa[1];
		$ex{$num}=1;
		unless(defined $pre){$pre=$aa[1];push @{$cm_type{$aa[1]}},$num;$at{$num}=$aa[1];}
		if($aa[1] != $pre){push @{$cm_type{$aa[1]}},$num;$at{$num}=$aa[1];}
		#print "map1\t$_";
		$pre=$aa[1];
	}
	#else{print "map2\t$_";}
	$num++;
}
close IN;

foreach my $key(keys %cm1){
	print STDERR "$key\t$cm1{$key}\n";
}

my %new_cm_type;
foreach my $cmm(keys %cm_type){
	if(@{$cm_type{$cmm}}>1){
		#need to reset the first one
		my @aa=@{$cm_type{$cmm}};
		my $dv=@aa;
		my $sort=$cm1{$cmm}-1;
		my $pre_cmm=$cm1_sort{$sort};
		my $block=$cmm-$pre_cmm;
		#print STDERR "$cmm\t$cm1{$cmm}\t$sort\t$pre_cmm\n";
		$block/=$dv;
		for(my $i=0;$i<@aa-1;$i++){
			my $cmmm=$cmm-$block;
			$new_cm_type{$cmmm}=[$aa[$i]];
			$all{$aa[$i]}[1]=$cmmm;
			$at{$aa[$i]}=$cmmm;
		}
		$new_cm_type{$cmm}=[$aa[-1]];
	}else{
		$new_cm_type{$cmm}=[$cm_type{$cmm}[0]];
	}
}

my @old=sort {$a<=>$b} keys %at;
my @new;
if($map_sort>0){@new=sort {$a<=>$b} keys %new_cm_type;}else{@new=sort {$b<=>$a} keys %new_cm_type;}
if(@old != @new){die "error in cm sort\t$#old\t$#new\n";}
my %refer;
for(my $i=0;$i<@old;$i++){
	$refer{$at{$old[$i]}}=$new[$i];
	#print "$at{$old[$i]} refer $new[$i]\n";
}

#foreach my $nn(sort {$a<=>$b} keys %all){
#	my $out=join "\t",@{$all{$nn}};
#	print "$out\n";
#}

foreach my $nn(sort {$a<=>$b} keys %all){
	if($ex{$nn}){
		#reset cM
		$all{$nn}[1]=$refer{$all{$nn}[1]};
		if($all{$nn}[0]=~/block/){
			my $block=$all{$nn}[0];
			foreach my $ctg(@{$block{$block}}){
				my $part1=$map1{$ctg}[1];
				my @part1=split /\s+/,$part1;
				$part1[1]=$all{$nn}[1];
				my $out1=join "\t",@part1;
				print "$out1\t$map2{$ctg}[1]\t$block\n";
			}
		}else{
			my $out=join "\t",@{$all{$nn}};
			print "$out\n";
		}
	}else{
		my $out=join "\t",@{$all{$nn}};
		print "$out\n";
	}
}





sub read_map{
	my ($map,$hash,$hash1,$cm)=@_;
	my $num=1;
	open IN,$map or die;
	while(<IN>){
		my @aa=split;chomp;
		$hash->{$aa[0]}=[$num,$_];
		$hash1->{$num}=$aa[0];
		for(my $i=1;$i<@aa;$i+=2){
			$cm->{$aa[$i]}=1;
		}
		$num++;
	}
	close IN;
}
