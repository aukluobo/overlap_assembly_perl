use strict;
use warnings;

#huangguodong@genomics.cn 2014.09.23

#this program is used to add the longest scaffold in the unused scafflod list based on the agp
#two main conditions:
#1.in the begining  and the ending 
#2.in the main gap between the two ends of one BAC
#the un-overlap gap will be treated as un-anchor scaffold if the genectic map can anchor more.
#two type conditions need to treat:
#1:same BAC id surround N gap
#2:neibor BAC id 1-->1 2-->1(30k) 2-->2;  note: 1-->2 means the block is overlap.

die "perl $0 <agp> <unused_scaffold_list.lst.len> <mtp_rel> <data_dir> <output_file>\ndepend on perl 5.16 or newer" if(@ARGV==0);

my ($agp,$unused_scf,$rel,$data_dir,$out)=@ARGV;
my %rel;my %bacSort;my %BACinCTG;my %last_sort;
read_rel($rel,\%rel,$data_dir,\%bacSort,\%BACinCTG);
my %unused_scf_len;my %bacTolen;my %exCTG;my %add_longest;
read_scf_len($unused_scf,\%unused_scf_len,\%bacTolen,\%BACinCTG,\%exCTG,\%add_longest);

my %block;my $num_c=1;my %bac12;
open IN,$agp or die;
open OUT,">$out.new.agp" or die; 

my %block_d;my $sort_d=1;my %row_sort;

while(<IN>){
	my @aa=split;
	my $re=$_;
	if($exCTG{$aa[0]}){
		if($re=~/bac/){
			#calculate the exist BAC scaffold length
			my $scf=$aa[5];my $id=$scf;$id=~s/_\d+$//;
			if($aa[-2]=~/(d\d+)_/){
				my $ascf=$';my $dd=$1;
				$ascf=~s/_\d+$//;
				$block_d{$ascf}{$dd}{$sort_d}=1;
			}
		}
		$block{$aa[0]}{$num_c}=[@aa];
		$row_sort{$sort_d}=$num_c;
		$num_c++;
		if(/scaffold/){$sort_d++;}
	}else{
		print OUT $re;
	}
	if($aa[4] eq "W"){
		my $id=$aa[5];$id=~s/_\d+$//;
		$bac12{$id}{$aa[-1]}=1;
	}
}
close IN;

my %head_tail;
foreach my $ctg(keys %block){
	my ($head,$tail)=(sort {$a<=>$b} keys $block{$ctg})[0,-1];
	my $head_bac=$block{$ctg}{$head}[5];my $tail_bac=$block{$ctg}{$tail}[5];
	$head_bac=~s/_\d+$//;$tail_bac=~s/_\d+$//;
	my $head_rel_sort=$bacSort{$head_bac};my $tail_rel_sort=$bacSort{$tail_bac};
	if($head_rel_sort != 0){$head_tail{$ctg}{"head"}=[1,$head_rel_sort];}
	my $last_sort=(sort {$a<=>$b} keys $rel{$ctg})[-1];
	if($tail_rel_sort != $last_sort){$head_tail{$ctg}{"tail"}=[1,$tail_rel_sort];print STDERR "$ctg\t$tail_rel_sort\n";}
}


#find the insert point
my %get_point;
foreach my $ctg(keys %block){
	my %add_point;
	foreach my $row_sort(sort {$a<=>$b} keys $block{$ctg}){
		my $record=$block{$ctg}{$row_sort}[6];
		if($record=~/scaffold/){
			my $before=$row_sort-1;my $after=$row_sort+1;
			my @aa=@{$block{$ctg}{$before}};my @bb=@{$block{$ctg}{$after}};
			my $ida=$aa[5];$ida=~s/_\d+$//;
			my $idb=$bb[5];$idb=~s/_\d+$//;
			if($ida eq $idb){
				if(($aa[-1]==2 && $aa[-1]>$bb[-1]) or (!$bac12{$ida}{"2"} && $aa[-1]>$bb[-1])){
					$add_point{$ida}{$row_sort}=1;
				}
			}
		}
	}
	#detete lone 1
	my %lone_1;
	my ($st,$ed)=(sort {$a<=>$b} keys $block{$ctg})[0,-1];
	foreach my $row_sort(sort {$a<=>$b} keys $block{$ctg}){
		my $record=$block{$ctg}{$row_sort}[5];
		if($record=~/bac/ && ($row_sort>$st && $row_sort<$ed)){
			my $pre=$row_sort-1;my $aft=$row_sort+1;
			if($block{$ctg}{$pre}[6]=~/scaffold/ && $block{$ctg}{$aft}[6]=~/scaffold/){
				my @aa=@{$block{$ctg}{$row_sort}};
				if($aa[-1]==1){
					my $ida=$aa[5];$ida=~s/_\d+$//;
					$lone_1{$ida}{$row_sort}=1;
				}
			}
		}
	}
	#d1d2
	foreach my $id_scf(keys %block_d){
		if($block_d{$id_scf}{"d1"} && $block_d{$id_scf}{"d2"}){
			my @d1=keys $block_d{$id_scf}{"d1"};my @d2=keys $block_d{$id_scf}{"d2"};
			foreach my $d1(@d1){
				foreach my $d2(@d2){
					if($d2-$d1==1){
						my $t_row=$row_sort{$d1};
						undef $add_point{$id_scf};delete $add_point{$id_scf};
						$add_point{$id_scf}{$t_row}=1;
						print "d1\td2\t$t_row\n";
						last;
					}
				}
			}
		}
	}
	#

	foreach my $bac_id(keys %add_point){
		if(keys $add_point{$bac_id} > 1){
			my $lone1=(sort {$a<=>$b} keys $lone_1{$bac_id})[0];
			my $keep;my @keep;
			foreach my $snt(sort {$a<=>$b} keys $add_point{$bac_id}){
				if($snt<$lone1){push @keep,$snt;}
			}
			if(@keep>0){
				$keep=$keep[-1];
				$get_point{$bac_id}{$keep}=1;
			}else{
				print STDERR "$bac_id no insert point\n";
			}
		}else{
			my $keep=(keys $add_point{$bac_id})[0];
			$get_point{$bac_id}{$keep}=1;
		}
	}
}
my %ex_point;
foreach my $bac_id(keys %get_point){
	print "$bac_id\thave get_point\n";
	$ex_point{$bac_id}=1;
}

my %haveadd;
foreach my $ctg(keys %block){
	my $st=1;my $ed=1;my $num=1;
	my ($first,$last)=(sort {$a<=>$b} keys $block{$ctg})[0,-1];
	foreach my $sort(sort {$a<=>$b} keys $block{$ctg}){
		if($first==$sort){
			#need to add but need to check the 1/2
			my @tmp=@{$block{$ctg}{$sort}};my $tmp_bac=$tmp[5];$tmp_bac=~s/_\d+$//;
			if($tmp[-1]==1){
				if($add_longest{$tmp_bac}){
					my $insert_scf=$add_longest{$tmp_bac}[1];
					my $insert_scf_len=$add_longest{$tmp_bac}[0];
					$haveadd{$tmp_bac}=1;
					$ed=$st+$insert_scf_len-1;
					print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddfirst\n";
					$st=$ed+1;$num++;
					$ed=$st+100-1;
					print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
					$st=$ed+1;$num++;
				}
			}
		}
		#if($last==$sort){}
		if($block{$ctg}{$sort}[4] eq "N"){			
			print "check\t$sort\t$block{$ctg}{$sort}[5]\n";
			my $pre=$sort-1;my $aft=$sort+1;
			my @pre=@{$block{$ctg}{$pre}};
			my @aft=@{$block{$ctg}{$aft}};
			my @now=@{$block{$ctg}{$sort}};
			my $pre_bac=$pre[5];$pre_bac=~s/_\d+$//;
			my $aft_bac=$aft[5];$aft_bac=~s/_\d+$//;
			#check condition
			if($bacSort{$aft_bac}-$bacSort{$pre_bac}<0){
					print "<0\n";
					my @aa=@{$block{$ctg}{$sort}};
					my $len=$aa[2]-$aa[1]+1;
					$ed=$st+$len-1;
					print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
					$st=$ed+1;$num++;
			}elsif($bacSort{$aft_bac}-$bacSort{$pre_bac}==0){
				#same BAC
				print "=0\t$pre_bac\t$aft_bac\t$aft[-1]\t$pre[-1]\n";
				#if(($aa[-1]==2 && $aa[-1]>$bb[-1]) or (!$check_123{$ida}{"2"} && $aa[-1]>$bb[-1])){
				#if($aft[-1]<$pre[-1]){
				if($get_point{$pre_bac}{$sort}){
					if(($pre[-1]==2 && $aft[-1]<$pre[-1]) or (!$bac12{$pre_bac}{"2"} && $aft[-1]<$pre[-1])){
						#condition 1
						print "condition 1\t$aft_bac\n";
						if($add_longest{$pre_bac}){
							my $now_gap_len=$now[5];
							my $new_pre_gap_len=100;my $new_aft_gap_len=100;
							if($now_gap_len != 100){
								#known gap
								my $sp=int($now_gap_len/2);
								if($sp>100){
									$new_pre_gap_len=$sp;$new_aft_gap_len=$sp;
								}
							}
							$ed=$st+$new_pre_gap_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_pre_gap_len\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
							my $insert_scf=$add_longest{$pre_bac}[1];
							my $insert_scf_len=$add_longest{$pre_bac}[0];
							$haveadd{$pre_bac}=1;
							$ed=$st+$insert_scf_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddc1\n";
							$st=$ed+1;$num++;
							$ed=$st+$new_aft_gap_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_aft_gap_len\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
						}else{
							print "=0=";
							my @aa=@{$block{$ctg}{$sort}};
							my $len=$aa[2]-$aa[1]+1;
							$ed=$st+$len-1;
							print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
							$st=$ed+1;$num++;
						}
					}else{
						print "=0==\n";
						my @aa=@{$block{$ctg}{$sort}};
						my $len=$aa[2]-$aa[1]+1;
						$ed=$st+$len-1;
						print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
						$st=$ed+1;$num++;
					}
				}else{
					my @aa=@{$block{$ctg}{$sort}};
					my $len=$aa[2]-$aa[1]+1;
					$ed=$st+$len-1;
					print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
					$st=$ed+1;$num++;
				}
			}elsif($bacSort{$aft_bac}-$bacSort{$pre_bac}==1){
					#neibor BAC . three condition
					my $flag=1;
					#1-->1
					if($pre[-1] == 1 && $aft[-1] == 1){
						#insert aft longest. aft must don't have 2
						if($bac12{$aft_bac}{"2"}){
							my @aa=@{$block{$ctg}{$sort}};
							my $len=$aa[2]-$aa[1]+1;
							$ed=$st+$len-1;
							print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
							$st=$ed+1;$num++;
						}else{
							if(!$ex_point{$aft_bac} && $add_longest{$aft_bac}){
								my $now_gap_len=$now[5];
								my $new_pre_gap_len=100;my $new_aft_gap_len=100;
								if($now_gap_len != 100){
									#known gap
									my $sp=int($now_gap_len/2);
									if($sp>100){
										$new_pre_gap_len=$sp;$new_aft_gap_len=$sp;
									}
								}
								$ed=$st+$new_pre_gap_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_pre_gap_len\tscaffold\tno\tna\n";
								$st=$ed+1;$num++;
								my $insert_scf=$add_longest{$aft_bac}[1];
								my $insert_scf_len=$add_longest{$aft_bac}[0];
								$haveadd{$aft_bac}=1;
								$ed=$st+$insert_scf_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd11\n";
								$st=$ed+1;$num++;
								$ed=$st+$new_aft_gap_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_aft_gap_len\tscaffold\tno\tna\n";
								$st=$ed+1;$num++;
							}else{
								my @aa=@{$block{$ctg}{$sort}};
								my $len=$aa[2]-$aa[1]+1;
								$ed=$st+$len-1;
								print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
								$st=$ed+1;$num++;
							}
						}
						$flag=0;
					}
					#2-->2
					if($pre[-1] == 2 && $aft[-1] == 2){
						#insert pre longest. pre must don't have 1
						if($bac12{$pre_bac}{"1"}){
							my @aa=@{$block{$ctg}{$sort}};
							my $len=$aa[2]-$aa[1]+1;
							$ed=$st+$len-1;
							print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
							$st=$ed+1;$num++;
						}else{
							if(!$ex_point{$pre_bac} && $add_longest{$pre_bac}){
								my $now_gap_len=$now[5];
								my $new_pre_gap_len=100;my $new_aft_gap_len=100;
								if($now_gap_len != 100){
									#known gap
									my $sp=int($now_gap_len/2);
									if($sp>100){
										$new_pre_gap_len=$sp;$new_aft_gap_len=$sp;
									}
								}
								$ed=$st+$new_pre_gap_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_pre_gap_len\tscaffold\tno\tna\n";
								$st=$ed+1;$num++;
								my $insert_scf=$add_longest{$pre_bac}[1];
								my $insert_scf_len=$add_longest{$pre_bac}[0];
								$haveadd{$pre_bac}=1;
								$ed=$st+$insert_scf_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd22\n";
								$st=$ed+1;$num++;
								$ed=$st+$new_aft_gap_len-1;
								print OUT "$ctg\t$st\t$ed\t$num\tN\t$new_aft_gap_len\tscaffold\tno\tna\n";
								$st=$ed+1;$num++;
							}else{
								my @aa=@{$block{$ctg}{$sort}};
								my $len=$aa[2]-$aa[1]+1;
								$ed=$st+$len-1;
								print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
								$st=$ed+1;$num++;
							}
						}
						$flag=0;
					}
					#2-->1
					if($pre[-1] == 2 && $aft[-1] == 1){
						#two BAC no overlap.change 10K gap. add both. the pre don't have 1 and the aft don't have 2
						if(!$ex_point{$pre_bac} && $add_longest{$pre_bac}){
							#add 100 gap
							$ed=$st+100-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
							#add pre
							my $insert_scf=$add_longest{$pre_bac}[1];
							my $insert_scf_len=$add_longest{$pre_bac}[0];
							$haveadd{$pre_bac}=1;
							$ed=$st+$insert_scf_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd21\n";
							$st=$ed+1;$num++;
						}
						#10k gap
						$ed=$st+10000-1;
						print OUT "$ctg\t$st\t$ed\t$num\tN\t10000\tscaffold\tno\tna\n";
						$st=$ed+1;$num++;
						#add aft
						if(!$ex_point{$aft_bac} && $add_longest{$aft_bac}){
							my $insert_scf=$add_longest{$aft_bac}[1];
							my $insert_scf_len=$add_longest{$aft_bac}[0];
							$haveadd{$aft_bac}=1;
							$ed=$st+$insert_scf_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd21\n";
							$st=$ed+1;$num++;
							#add 100 gap
							$ed=$st+100-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
						}
						$flag=0;
					}
					if($flag){
						my @aa=@{$block{$ctg}{$sort}};
						my $len=$aa[2]-$aa[1]+1;
						$ed=$st+$len-1;
						print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
						$st=$ed+1;$num++;
					}
			}else{
					#BAC is missing. if the BAC is exists, we insert longest.
					my $miss_bac_num=$bacSort{$aft_bac}-$bacSort{$pre_bac}-1;
					my $cu=1;my $flag=0;
					while($cu<=$miss_bac_num){
						#
						my $bacSort=$bacSort{$pre_bac}+$cu;
						my $bac_miss=$rel{$ctg}{$bacSort}[0];
						if(!$ex_point{$bac_miss} && $add_longest{$bac_miss}){
							#add longest
							$flag=1;
							$ed=$st+100-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
							my $insert_scf=$add_longest{$bac_miss}[1];
							my $insert_scf_len=$add_longest{$bac_miss}[0];
							$haveadd{$bac_miss}=1;
							$ed=$st+$insert_scf_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddmiss\n";
							$st=$ed+1;$num++;
						}
						$cu++;
					}
					my $bacmiss=$aft_bac;
					if($aft[-1]==1){
						if(!$ex_point{$bacmiss} && $add_longest{$bacmiss}){
							$flag=1;
							$ed=$st+100-1;
							print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
							$st=$ed+1;$num++;
							my $insert_scf=$add_longest{$bacmiss}[1];
							my $insert_scf_len=$add_longest{$bacmiss}[0];
							$haveadd{$bacmiss}=1;
							$ed=$st+$insert_scf_len-1;
							print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddmissaft\n";
							$st=$ed+1;$num++;
						}
					}
					if($flag){
						$ed=$st+100-1;
						print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
						$st=$ed+1;$num++;
					}else{
						my @aa=@{$block{$ctg}{$sort}};
						my $len=$aa[2]-$aa[1]+1;
						$ed=$st+$len-1;
						print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\n";
						$st=$ed+1;$num++;
					}
			}
		}else{
			my @aa=@{$block{$ctg}{$sort}};
			my $len=$aa[2]-$aa[1]+1;
			$ed=$st+$len-1;
			#print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\t$aa[10]\n";
			$aa[1]=$st;$aa[2]=$ed;my $out=join "\t",@aa;
			print OUT "$out\n";
			$st=$ed+1;$num++;
			if($last==$sort){
				#need to add last
				my @tmp=@{$block{$ctg}{$sort}};my $tmp_bac=$tmp[5];$tmp_bac=~s/_\d+$//;
				if($tmp[-1]==2){
					if(!$ex_point{$tmp_bac} && $add_longest{$tmp_bac}){
						my $insert_scf=$add_longest{$tmp_bac}[1];
						my $insert_scf_len=$add_longest{$tmp_bac}[0];
						$haveadd{$tmp_bac}=1;
						$ed=$st+100-1;
						print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
						$st=$ed+1;$num++;
						$ed=$st+$insert_scf_len-1;
						print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddlast\n";
						$st=$ed+1;$num++;	
					}
				}
			}
		}
	}
}
close IN;
close OUT;
open IN,"$out.new.agp" or die;
open OUT,">$out.new.agp.ht" or die; 
my %newblock;my $newnum=1;
while(<IN>){
	my @aa=split;
	if($head_tail{$aa[0]}){
		$newblock{$aa[0]}{$newnum}=[@aa];
		$newnum++;
	}else{
		print OUT $_;
	}
}
close IN;

#the second first and the second last also need to add in the scaffold
#%haveadd
foreach my $ctg(keys %newblock){
	my $st=1;my $ed=1;my $num=1;
	if($head_tail{$ctg}{"head"}[0]){
		#add head
		my $head_sort=$head_tail{$ctg}{"head"}[1];
		my $cu=0;my $flag=0;
		while($cu<$head_sort){
			my $need_bac=$rel{$ctg}{$cu}[0];
			if($add_longest{$need_bac}){
				my $insert_scf=$add_longest{$need_bac}[1];
				my $insert_scf_len=$add_longest{$need_bac}[0];
				$ed=$st+$insert_scf_len-1;
				print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd_head\n";
				$st=$ed+1;$num++;
				$ed=$st+100-1;
				print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
				$st=$ed+1;$num++;
			}
			$cu++;
		}
	}
	foreach my $sort(sort {$a<=>$b} keys $newblock{$ctg}){
		my @aa=@{$newblock{$ctg}{$sort}};
		my $len=$aa[2]-$aa[1]+1;
		$ed=$st+$len-1;
		#print OUT "$aa[0]\t$st\t$ed\t$num\t$aa[4]\t$aa[5]\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\n";
		$aa[1]=$st;$aa[2]=$ed;$aa[3]=$num;
		my $out=join "\t",@aa;
		print OUT "$out\n";
		$st=$ed+1;$num++;
	}
	if($head_tail{$ctg}{"tail"}[0]){
		#add tail
		print STDERR "$ctg\taddtail\n";
		my $tail_sort=$head_tail{$ctg}{"tail"}[1];
		my $cu=$tail_sort+1;my $flag=0;
		my $last_sort=(sort {$a<=>$b} keys $rel{$ctg})[-1];
		while($cu<=$last_sort){
			my $need_bac=$rel{$ctg}{$cu}[0];
			if($add_longest{$need_bac}){
				print STDERR "$need_bac add tail\n";
				$ed=$st+100-1;
				print OUT "$ctg\t$st\t$ed\t$num\tN\t100\tscaffold\tno\tna\n";
				$st=$ed+1;$num++;
				my $insert_scf=$add_longest{$need_bac}[1];
				my $insert_scf_len=$add_longest{$need_bac}[0];
				$ed=$st+$insert_scf_len-1;
				print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\tadd_tail\n";
				$st=$ed+1;$num++;
			}
			$cu++;
		}
	}
}

close OUT;


sub read_scf_len{
 	my ($len,$hash,$tolen,$BACinCTG,$exCTG,$add_longest)=@_;
 	open IN,$len or die;
 	while(<IN>){
   		my @aa=split;
    	$hash->{$aa[0]}=$aa[1];
    	my $id=$aa[0];$id=~s/_\d+$//;
    	#$aa[0]=~s/_\d+$//;
    	$tolen->{$id}+=$aa[1];
    	if($add_longest->{$id}){
    		if($aa[1]>$add_longest->{$id}[0]){
    			$add_longest->{$id}[0]=$aa[1];$add_longest->{$id}[1]=$aa[0];
    		}
    	}
    	else{
    		$add_longest->{$id}[0]=$aa[1];$add_longest->{$id}[1]=$aa[0];
    	}
    	if($BACinCTG->{$id}){
    		my $ctg=$BACinCTG->{$id};
    		if($exCTG->{$ctg}){next;}
    		$exCTG->{$ctg}=1;
    	}
 	}
	close IN;
}

sub read_rel{
  	my ($rel,$hash,$dir,$bacSort,$BACinCTG)=@_;
  	open IN,$rel or die;
  	while(<IN>){
    	if(/bac/){
     	my @aa=split;
      		if($aa[-1]=~/(CTG\d+)_(bac.*)_(\d+)$/){
		        my $ctg=$1;my $bac=$2;my $sort=$3;
				my $ex=0;
				if(-s "$dir/$ctg/$aa[-1].fa"){$ex=1;}
				#bac_id exists_or_not expect_len overlap    $aa[1] is expect total length; $aa[3] is expect overlap length
				$hash->{$ctg}{$sort}=[$bac,$ex,$aa[1],$aa[3]];
				$bacSort->{$bac}=$sort;
				$BACinCTG->{$bac}=$ctg;
      		}
    	}
  	}
  	close IN;
}

