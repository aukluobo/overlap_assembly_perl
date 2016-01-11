use strict;use warnings;

#huangguodong@genomics.cn 20141014

die "perl $0 <merge_map> <raw_map1> <raw_map2>" if(@ARGV==0);

my ($map,$map1,$map2)=@ARGV;

my %map1;my %map2;my %map1_ctg_sort;my %map2_ctg_sort;my %cm1;my %cm2;
read_map($map1,\%map1,\%map1_ctg_sort,\%cm1);
read_map($map2,\%map2,\%map2_ctg_sort,\%cm2);

my %cm1_sort;my %cm2_sort;my $c_n=1;
my @key_cm1=sort {$b<=>$a} keys %cm1;
foreach(@key_cm1){$cm1_sort{$c_n}=$_;$cm1{$_}=$c_n;$c_n++;}
$c_n=1;
my @key_cm2=sort {$a<=>$b} keys %cm2;
foreach(@key_cm2){$cm2_sort{$c_n}=$_;$cm2{$_}=$c_n;$c_n++;}

open IN,$map or die;
my %cm_sort;my $cm_num=1;
my %all;my $num=1;
while(<IN>){
	my @aa=split;my $rec=$_;chomp $rec;
	unless($cm_sort{$aa[1]}){$cm_sort{$aa[1]}=$cm_num;$cm_num++;}
	my $len=@aa;my $big=0;
	for(my $i=3;$i<@aa;$i+=2){
		if($aa[$i]>$big){$big=$aa[$i];}
	}
	$all{$aa[1]}{$aa[0]}=[$num,$len,$big,$rec];
	$num++;
}
close IN;

foreach my $cM(sort {$b<=>$a} keys %all){
	my @ctg=sort {$all{$cM}{$a}[0]<=>$all{$cM}{$b}[0]} keys $all{$cM};
	if(@ctg==1){next;}
	my %len;my %big;
	for(my $i=0;$i<@ctg;$i++){
		if($all{$cM}{$ctg[$i]}[2]){$len{$ctg[$i]}=1;$big{$all{$cM}{$ctg[$i]}[2]}=1;}
	}
	if(keys %big == 1){
		my $change=0;my $st_co;
		for(my $i=0;$i<@ctg;$i++){
			if($len{$ctg[$i]}){
				#delete the aft-fix
				my @aa=split /\s+/,$all{$cM}{$ctg[$i]}[3];
				my $recrec="$aa[0]\t$aa[1]\t$aa[2]";
				$all{$cM}{$ctg[$i]}[3]=$recrec;
			}
		}
	}else{
		my %ins;my %inb;my $pre_big;
		for(my $i=0;$i<@ctg-1;$i++){
			my $first=$ctg[$i];my $second=$ctg[$i+1];
			unless($pre_big){
				if($all{$cM}{$first}[2]>0){
					$pre_big=$all{$cM}{$first}[2];
				}else{
					next;
				}
			}
			if($all{$cM}{$first}[2]>0){unless($pre_big){$pre_big=$all{$cM}{$first}[2];}}
			if($all{$cM}{$first}[1]<=$all{$cM}{$second}[1]){
				#print STDERR "$all{$cM}{$first}[3]\n$all{$cM}{$second}[3]\n";
				if($all{$cM}{$first}[2]<$all{$cM}{$second}[2]){
					if($all{$cM}{$second}[2]>0 && $all{$cM}{$second}[2]>$pre_big){
						#
						$ins{$i+1}=1;
						$inb{$i+1}=1;
						$pre_big=$all{$cM}{$second}[2];
					}
				}
			}else{
				#got change
				$ins{$i+1}=1;
			}
		}
		my $len=0;
		for(my $i=0;$i<@ctg;$i++){
			my $now_len=$all{$cM}{$ctg[$i]}[1];
			my @aa=split /\s+/,$all{$cM}{$ctg[$i]}[3];
			my $append="\t$aa[1]\t$aa[2]";
			if($ins{$i}){
				unless($len){$len=$all{$cM}{$ctg[$i-1]}[1];}
				my $pre_ctg=$ctg[$i-1];
				#if($inb{$i}){$len+=2;}
				if($inb{$i} or ($map2{$pre_ctg}[2] && $map2{$ctg[$i]}[2] && ($map2{$pre_ctg}[2] ne $map2{$ctg[$i]}[2]))){
					$len+=2;
				}
			}
			if($len && $now_len>$len){
				#
				print STDERR "trace $ctg[$i]\n";
				$len=$now_len;
			}
			while($now_len<$len){
					print STDERR "$ctg[$i]\t$now_len\t$len\n";
					$all{$cM}{$ctg[$i]}[3].=$append;
					$now_len+=2;
			}
			
		}

	}
}

foreach my $cM(sort {$b<=>$a} keys %all){
	foreach my $ctg(sort {$all{$cM}{$a}[0]<=>$all{$cM}{$b}[0]} keys $all{$cM}){
		print "$all{$cM}{$ctg}[3]\n";
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
		my $out;my @out;
		for(my $i=1;$i<@aa;$i+=2){
			$cm->{$aa[$i]}=1;
			push @out,$aa[$i];
		}
		$out=join "\t",@out;
		push @{$hash->{$aa[0]}},$out;
		$num++;
	}
	close IN;
}