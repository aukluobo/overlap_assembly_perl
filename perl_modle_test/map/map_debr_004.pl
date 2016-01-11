use strict;use warnings;

#huangguodong@genomics.cn 20141014

die "perl $0 <merge_map> <raw_map1> <raw_map2>" if(@ARGV==0);

my ($map,$map1,$map2)=@ARGV;

my %map1;my %map2;my %map1_ctg_sort;my %map2_ctg_sort;my %cm1;my %cm2;
read_map($map1,\%map1,\%map1_ctg_sort,\%cm1);
read_map($map2,\%map2,\%map2_ctg_sort,\%cm2);

my %cm1_sort;my %cm2_sort;my $c_n=1;
my @key_cm1=sort {$a<=>$b} keys %cm1;
foreach(@key_cm1){$cm1_sort{$c_n}=$_;$cm1{$_}=$c_n;}
$c_n=1;
my @key_cm2=sort {$a<=>$b} keys %cm2;
foreach(@key_cm2){$cm2_sort{$c_n}=$_;$cm2{$_}=$c_n;}

my %map;my %block;my $blo=1;my %exblock;my %map_ctg_sort;
#find the insert block by map2 and recalculate the cM
open IN,$map or die;
my $num=1;
while(<IN>){
	my @aa=split;chomp;my $rec=$_;
	my $ctg=$aa[0];
	$map{$ctg}=[$num,$rec];
	$map_ctg_sort{$num}=$ctg;
	if($map2{$ctg}[1] && $aa[1]==(split /\s+/,$map2{$ctg}[1])[1]){
		#rec show map2
		$block{$blo}{$num}=$rec;
		$exblock{$num}=$blo;
	}else{
		$blo++;
	}
	$num++;
}
close IN;
my %block_new_cm;
foreach my $bb(sort {$a<=>$b} keys %block){
	my @num=sort {$a<=>$b} keys $block{$bb};

	my $pre=$num[0]-1;my $aft=$num[-1]+1;
	my $pre_cm;my $aft_cm;my $pre_ctg;my $aft_ctg;
	if($pre<1){
		#no pre
		$aft_cm=$map{$map_ctg_sort{$aft}}[1];$aft_cm=(split /\s+/,$aft_cm)[1];
		$pre_cm=$aft_cm;
		$pre_ctg=0;$aft_ctg=$map_ctg_sort{$aft};
	}else{
		$pre_cm=$map{$map_ctg_sort{$pre}}[1];$pre_cm=(split /\s+/,$pre_cm)[1];
		$aft_cm=$map{$map_ctg_sort{$aft}}[1];$aft_cm=(split /\s+/,$aft_cm)[1];
		$pre_ctg=$map_ctg_sort{$pre};$aft_ctg=$map_ctg_sort{$aft};
	}
	#check if map2 have different block. if pre_cm == aft_cm, then the aft_cm need to redefined by next cm in map1.
	my %block2;my %type;my $ao=1;my %type_sort;
	foreach my $nn(@num){
		my @aa=split /\s+/,$block{$bb}{$nn};
		my $type;my $type_out;my $keep=0;
		for(my $i=1;$i<@aa;$i+=2){
			if($aa[$i]=~/a/){next;}
			$block2{$aa[$i]}=1;
			$type.="$aa[$i]\t";
			$type_out.="$aa[$i]\t$aa[$i+1]\t";$keep++;
		}
		push @{$type{$type}{$nn}},$aa[0];
		unless($type_sort{$type}){$type_sort{$type}=[$ao,$type_out,$keep];$ao++;}
	}
	if($pre_ctg && $map2{$pre_ctg}[1]){
		my @tmp_cm=split /\s+/,$map2{$pre_ctg}[1];
	}
	if($aft_ctg && $map2{$aft_ctg}[1]){
		my @tmp_cm=split /\s+/,$map2{$aft_ctg}[1];
	}
	my $block2=keys %block2;
	my $cha=$aft_cm-$pre_cm;
	my $dv=$cha/($block2+1);
	#print "$cha\t$pre_cm\t$aft_cm\n";
	if($cha==0){
		my $ty_cu=1;my $ext=1;my $last_ext;
		foreach my $ty(sort {$type_sort{$a}[0]<=>$type_sort{$b}[0]} keys %type_sort){
			unless($last_ext){$last_ext=$type_sort{$ty}[2];}
			foreach my $nn(sort {$a<=>$b} keys $type{$ty}){
				foreach my $pctg(@{$type{$ty}{$nn}}){
					my $out=1;
					$block_new_cm{$pctg}="$pre_cm\t1/0\t$type_sort{$ty}[1]\t";
					$ext=$type_sort{$ty}[2];
					my $tmp="";
					while($out<$ty_cu){$tmp.="$type_sort{$ty}[1]\t";$ext+=$type_sort{$ty}[2];$out++;}
					if($ext>$last_ext){
						if($ext-$last_ext>1){
							#delete some
							my $over=$ext-$last_ext-1;
							$ext-=$over;
							$over*=2;
							my @tmp=split /\s+/,$tmp;
							my $new_tmp;my $tmp_c=0;
							for(my $i=0;$i<@tmp-$over;$i++){
								$new_tmp.="$tmp[$i]\t";$tmp_c++;
							}
							$tmp=$new_tmp;
						}
					}elsif($ext<$last_ext){
						#supply more
						my $over=$last_ext-$ext+1;
						while($over>0){$tmp.="$type_sort{$ty}[1]\t";$over--;}
					}
					$block_new_cm{$pctg}.="$tmp\n";
					##$block_new_cm{$pctg}="$pre_cm\t1\t$block{$bb}{$nn}\n";
					#print "$pctg\t$block_new_cm{$pctg}\n";
					##print "$pctg\t$pre_cm\t$map_ctg_sort{$pre}\t$aft_cm\t$new_cm\n";
					$last_ext=$ext;
				}
			}
			$ty_cu++;
		}
	}else{
		my $new_cm=$pre_cm;
		foreach my $ty(sort {$type_sort{$a}[0]<=>$type_sort{$b}[0]} keys %type_sort){
			$new_cm+=$dv;
			#if map1 eq map2 in pre then no change
			foreach my $nn(sort {$a<=>$b} keys $type{$ty}){
				foreach my $pctg(@{$type{$ty}{$nn}}){
					$block_new_cm{$pctg}="$new_cm\t1/1\n";
					#print "$pctg\t$pre_cm\t$map_ctg_sort{$pre}\t$aft_cm\t$new_cm\n";
				}
			}
		}
	}
}


foreach my $ctg(sort {$map{$a}[0]<=>$map{$b}[0]} keys %map){
	if($exblock{$map{$ctg}[0]}){
		#my $bl=$exblock{$map{$ctg}[0]};
		#print "$ctg\t$block_new_cm{$ctg}\t1\tnew\t$map{$ctg}[1]\n";
		print "$ctg\t$block_new_cm{$ctg}";
	}else{
		#print "$map{$ctg}[1]\n";
		my @aa=split /\s+/,$map{$ctg}[1];
		for(my $i=0;$i<@aa;$i++){
			if($i>0 && ($aa[$i]=~/CTG/)){next;}
			if($aa[$i]=~/add/ or $aa[$i]=~/block/){last;}
			print "$aa[$i]\t";
		}
		print "\n";
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
		$num++;
		for(my $i=1;$i<@aa;$i+=2){
			$cm->{$aa[$i]}=1;
		}
	}
	close IN;
}