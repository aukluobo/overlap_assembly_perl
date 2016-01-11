use strict;
use warnings;

package Resolve_AGP_1_change;

#huangguodong@genomics.cn

sub re_agp{
	my ($self,$out,$all_scf,$bacTolen,$rel_ctg,$ctg,$agp)=@_;#all the scalar is hash reference
	#firstly, we need to resolve the repeat scaffold and make the order of overlap
	#make sure the scaffold was ordered by the CTG sort
	my %repeat;my %repeat_tmp;
	foreach my $sort(sort {$a<=>$b} keys %{$out}){
		foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$sort}){
			if($out->{$sort}{$ctg_sort}[0][0] && @{$out->{$sort}{$ctg_sort}[0]}>1){
				$repeat_tmp{$out->{$sort}{$ctg_sort}[0][0]}{$sort}++;
			}
			if($out->{$sort}{$ctg_sort}[1][0] && @{$out->{$sort}{$ctg_sort}[1]}>1){
				$repeat_tmp{$out->{$sort}{$ctg_sort}[1][0]}{$sort}++;
			}
		}
	}
	foreach my $key(keys %repeat_tmp){
		if((keys $repeat_tmp{$key})>1){$repeat{$key}=1;}
	}
	my @key_sort=sort {$a<=>$b} keys %{$out};
	my $key_sort=@key_sort;
	#my $pre_sort=shift @key_sort;
	#my %pre_scf;
	#my %order;my $order=1;
	my %ovl;my %revers_ovl;my %pre_ovl;my %result;my $cal_result=1;my %rep_data1;my %rep_data2;
	my %next;my $cal_next=1;
	my %result_pre2;my @result_pre2;
	my $cal_agp=1;my $agp_st=1;my $agp_ed=1;
	#the sort was start from 1 depend on the Paragraph module.
	print "start to output agp\n";
	open OUT,">$agp" or die "no agp\n";
	if(@key_sort==1){
		#only two BAC ouput all record;
		print STDERR "only one or two BACs please select the longest scaffold please\n";
		my $aft_s=$key_sort[0];
		foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$aft_s}){
			if($out->{$aft_s}{$ctg_sort}[0][0] && @{$out->{$aft_s}{$ctg_sort}[0]}>1){
				if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
					#ovl check st ed
					my $stat1=$out->{$aft_s}{$ctg_sort}[0][1];
					my $stat2=$out->{$aft_s}{$ctg_sort}[1][1];
					my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];
					my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];
					my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
					my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
					if($stat1 eq "st"){
						$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
						if($stat2 eq "ed"){
							$typ2="-";$st2=1;$ed2=$bb[2]-1;
						}else{
							$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
						}
					}else{
						$st1=1;$ed1=$aa[2]-1;$typ1="+";
						if($stat2 eq "st"){
							$typ2="+";
							$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
						}else{
							$typ2="-";
							$st2=1;$ed2=$bb[2]-1;
						}
					}
					if($ovl1>$ovl2){
						#$onetwo=1;
						$st1=1;$ed1=$all_scf->{$scf1};
						$agp_ed=$agp_st+$all_scf->{$scf1}-1;
						print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
						$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
						$agp_ed=$agp_st+$len2-1;
						print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
						$cal_agp++;$agp_st=$agp_ed+1;
						$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
						$agp_st=$agp_ed+1;$cal_agp++;
					}else{
						#$onetwo=2;
						$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
						$agp_ed=$agp_st+$len1-1;
						print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
						$cal_agp++;$agp_st=$agp_ed+1;
						$agp_ed=$agp_st+$all_scf->{$scf2}-1;
						print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
						$cal_agp++;$agp_st=$agp_ed+1;
						$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";	
						print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
						$agp_st=$agp_ed+1;$cal_agp++;
					}
				}else{
					#only 0
					my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
					my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
					if($stat1 eq "st"){$typ="-";}
					my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
					print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
					print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
					$cal_agp++;$agp_st=$agp_ed+1;
					$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
					print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
					$agp_st=$agp_ed+1;$cal_agp++;
				}
			}else{
				#check if only 1
				if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
					my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
					my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
					if($stat2 eq "ed"){$typ="-";}
					my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
					print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
					print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
					$cal_agp++;$agp_st=$agp_ed+1;
					$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
					print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
					$agp_st=$agp_ed+1;$cal_agp++;
				}
			}
		}
	}
	my $last=0;
	for(my $i=0;$i<$key_sort-1;$i++){
		print "work\t";
		my $first=$key_sort[$i];my $second=$key_sort[$i+1];
		print "$first\t$second\n";
		my $last_check=$i+1;
		if($last_check==$key_sort-1){$last=1;}
		if($second-$first==1){
			print "well\n";
			my $check_scf=$rel_ctg->{$first}[0];print "$check_scf is the repeat check scaffold\n";
			my @rep=grep /${check_scf}_/,(keys %repeat);
			if(@rep>0){
				#got repeat scaffold
				my %rep_scf;my %rep_scf2;
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$first}){
					#just check 1!
					if($out->{$first}{$ctg_sort}[1][0] && @{$out->{$first}{$ctg_sort}[1]}>1){
						#my $scaffold1=$out->{$first}{$ctg_sort}[0][0];
						my $scaffold1;
						my $scaffold2=$out->{$first}{$ctg_sort}[1][0];
						if($scaffold2=~/${check_scf}_/){
							#my @tmp=grep {$_ eq $scaffold2} @rep;
							#if(@tmp>0){$rep_scf{$scaffold2}=$scaffold1;}
							if($repeat{$scaffold2}){
								print "prepare data for $scaffold2\n";
								if($out->{$first}{$ctg_sort}[0][0] && @{$out->{$first}{$ctg_sort}[0]}>1){
									$scaffold1=$out->{$first}{$ctg_sort}[0][0];
									$rep_scf{$scaffold2}=$scaffold1;
									$rep_data1{$scaffold2}=[@{$out->{$first}{$ctg_sort}[1]}];
									$rep_data1{$scaffold1}=[@{$out->{$first}{$ctg_sort}[0]}];
								}else{
									$rep_scf{$scaffold2}=1;
									$rep_data1{$scaffold2}=[@{$out->{$first}{$ctg_sort}[1]}];
								}
							}else{
								if($out->{$first}{$ctg_sort}[0][0] && @{$out->{$first}{$ctg_sort}[0]}>1){
									my @aa=@{$out->{$first}{$ctg_sort}[0]};my $fout=join "\t",@aa;print "d1\t$fout\n";
									$result{"1"}{$cal_result}[0]=[@aa];
								}
								if($out->{$first}{$ctg_sort}[1][0] && @{$out->{$first}{$ctg_sort}[1]}>1){
									my @aa=@{$out->{$first}{$ctg_sort}[1]};my $fout=join "\t",@aa;print "d1\t$fout\n";
									$result{"1"}{$cal_result}[1]=[@aa];
								}
								$cal_result++;
							}
							if($scaffold1 && $scaffold2){$ovl{$scaffold1}=$scaffold2;$revers_ovl{$scaffold2}=$scaffold1;}
						}
					}else{
						#be careful with the con scaffold in 0 
						if($out->{$first}{$ctg_sort}[0][0] && @{$out->{$first}{$ctg_sort}[0]}>1){
							my @aa=@{$out->{$first}{$ctg_sort}[0]};my $fout=join "\t",@aa;
							$result{"1"}{$cal_result}[0]=[@aa];
							$cal_result++;
						}
					}
				}
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$second}){
					#both: element or ovl. but just need to check 0
					print "second have element\n";
					if($out->{$second}{$ctg_sort}[0][0] && @{$out->{$second}{$ctg_sort}[0]}>1){
						my $scaffold1=$out->{$second}{$ctg_sort}[0][0];
						print $scaffold1,"\n";
						if($scaffold1=~/${check_scf}_/){
							if($repeat{$scaffold1}){
								print "prepare data for $scaffold1\n";
								if($out->{$second}{$ctg_sort}[1][0] && @{$out->{$second}{$ctg_sort}[1]}>1){
									#ovl
									my $scaffold2=$out->{$second}{$ctg_sort}[1][0];
									$rep_scf2{$scaffold1}=$scaffold2;
									$rep_data2{$scaffold1}=[@{$out->{$second}{$ctg_sort}[0]}];
									$rep_data2{$scaffold2}=[@{$out->{$second}{$ctg_sort}[1]}];
									$ovl{$scaffold1}=$scaffold2;$revers_ovl{$scaffold2}=$scaffold1;
								}else{
									$rep_scf2{$scaffold1}=1;
									$rep_data2{$scaffold1}=[@{$out->{$second}{$ctg_sort}[0]}];
								}
							}else{
								my @aa=@{$out->{$second}{$ctg_sort}[0]};my $fout=join "\t",@aa;print "d2\t$fout\n";
								$result{"2"}{$cal_result}[0]=[@aa];
								if($out->{$second}{$ctg_sort}[1][0] && @{$out->{$second}{$ctg_sort}[1]}>1){
									@aa=@{$out->{$second}{$ctg_sort}[1]};$fout=join "\t",@aa;print "d2\t$fout\n";
									$result{"2"}{$cal_result}[1]=[@aa];
								}
								$result{"2"}{$cal_result}[2]=$second;
								print "update result 2 $cal_result\t$second\t$out->{$second}{$ctg_sort}[0][0]\n";
								$cal_result++;
							}
						}
					}else{
						#be careful with the con scaffold in 1
						if($out->{$second}{$ctg_sort}[1][0] && @{$out->{$second}{$ctg_sort}[1]}>1){
							my @aa=@{$out->{$second}{$ctg_sort}[1]};
							$result{"2"}{$cal_result}[0]=();
							$result{"2"}{$cal_result}[1]=[@aa];
							$result{"2"}{$cal_result}[2]=$second;
							my $fout=join "\t",@aa;print "d2\t$fout\n";
							print "update result 2 $cal_result\t$second\t$out->{$second}{$ctg_sort}[1][0]\n";
							$cal_result++;
						}
					}
				}
				#deal result 1
				print "result 1\n";
				if($result{"1"}){
					foreach my $key(sort {$a<=>$b} keys $result{"1"}){
						#may need check some pre result{2}
						if($result{"1"}{$key}[0][0]){
							if($repeat{$result{"1"}{$key}[0][0]}){
								next;
							}
						}
						if($result{"1"}{$key}[1][0]){
							if($repeat{$result{"1"}{$key}[1][0]}){next;}
						}
						#output the 0 single
						if($result{"1"}{$key}[0][0] && @{$result{"1"}{$key}[0]}>1){
							if($result{"1"}{$key}[1][0] && @{$result{"1"}{$key}[1]}>1){
								#ovl check st ed
								next;
							}else{
								#con
								 my $st1=1;my $ed1=$all_scf->{$result{"1"}{$key}[0][0]};my $typ1="+";
								 my $len1=$ed1-$st1+1;$agp_ed=$agp_st+$len1-1;
								 my $stat1=$result{"1"}{$key}[0][1];
								 if($stat1 eq "st"){$typ1="-";}
								 print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$result{1}{$key}[0][0]\t$st1\t$ed1\t$typ1\n";
								 print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$result{1}{$key}[0][0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
								 $cal_agp++;$agp_st=$agp_ed+1;
								 $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								 print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								 $agp_st=$agp_ed+1;$cal_agp++;
							}
						}else{
							next;
						}
					}
					foreach my $key(sort {$a<=>$b} keys $result{"1"}){
						#may need check some pre result{2}
						if($result{"1"}{$key}[0][0]){
							if($repeat{$result{"1"}{$key}[0][0]}){
								next;
							}
						}
						if($result{"1"}{$key}[1][0]){
							if($repeat{$result{"1"}{$key}[1][0]}){next;}
						}
						#output the ovl 
						if($result{"1"}{$key}[0][0] && @{$result{"1"}{$key}[0]}>1){
							if($result{"1"}{$key}[1][0] && @{$result{"1"}{$key}[1]}>1){
								#ovl check st ed
								my $stat1=$result{"1"}{$key}[0][1];
								my $stat2=$result{"1"}{$key}[1][1];
								my @aa=@{$result{"1"}{$key}[0]};
								my @bb=@{$result{"1"}{$key}[1]};
								my $scf1=$aa[0];my $scf2=$bb[0];
								my $ovl1=$bb[3]-$bb[2]+1;
								my $ovl2=$bb[3]-$bb[2]+1;
								my ($st1,$ed1,$st2,$ed2,$len1,$len2,$type1,$type2);
								if($stat1 eq "st"){
									$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$type1="-";
									if($stat2 eq "st"){
										$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};$type2="+";
									}else{
										$st2=1;$ed2=$bb[2]-1;$type2="-";
									}
								}else{
									$st1=1;$ed1=$aa[2]-1;$type1="+";
									if($stat2 eq "st"){
										$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};$type2="+";
									}else{
										$st2=1;$ed2=$bb[2]-1;$type2="-";
									}
								}
								if($ovl1>$ovl2){
									$st1=1;$ed1=$all_scf->{$aa[0]};
									$len1=$ed1-$st1+1;$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$type1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$type1\tnod\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$len2=$ed2-$st2+1;$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$type2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$type2\tnod\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									$st2=1;$ed2=$all_scf->{$bb[0]};
									$len1=$ed1-$st1+1;$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$type1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$type1\tnod\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$len2=$ed2-$st2+1;$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$type2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$type2\tnod\t2\n";
									#insert 100bp Ns
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;
									$cal_agp++;
								}
							}else{
								#con
								next;
							}
						}
					}
					foreach my $key(sort {$a<=>$b} keys $result{"1"}){
						#may need check some pre result{2}
						if($result{"1"}{$key}[0][0]){
							if($repeat{$result{"1"}{$key}[0][0]}){
								next;
							}
						}
						if($result{"1"}{$key}[1][0]){
							if($repeat{$result{"1"}{$key}[1][0]}){next;}
						}
						if($result{"1"}{$key}[0][0] && @{$result{"1"}{$key}[0]}>1){
							if($result{"1"}{$key}[1][0] && @{$result{"1"}{$key}[1]}>1){
								#ovl check st ed
								next;
							}else{
								#con
								next;
							}
						}else{
							#output 1 single
							if($result{"1"}{$key}[1][0] && @{$result{"1"}{$key}[1]}>1){
								my $st1=1;my $ed1=$all_scf->{$result{"1"}{$key}[1][0]};my $typ1="+";
								my $stat1=$result{"1"}{$key}[1][1];
								if($stat1 eq "ed"){$typ1="-";}
								my $len1=$ed1-$st1+1;$agp_ed=$agp_st+$len1-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$result{1}{$key}[1][0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$result{1}{$key}[1][0]\t$st1\t$ed1\t$typ1\tnod\t2\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
						}
					}
					$result{"1"}=();
				}
				#deal repeat region
				#rep_scf rep_scf2
				print "repeat\n";
				#calculate the length to decide which part is better to place
				my %decide;
				foreach my $scf_rep(@rep){
					print "check repeat\t$scf_rep\n";
					my @aa=@{$rep_data1{$scf_rep}};
					print "rep_data1 get\n";
					my @bb=@{$rep_data2{$scf_rep}};
					print "rep_data2 get\n";
					my $len1=$aa[3]-$aa[2]+1;my $len2=$bb[3]-$bb[2]+1;
					if($len1>$len2){$decide{"1"}{$scf_rep}=$len1;}else{$decide{"2"}{$scf_rep}=$len2;}
				}
				#output ovl0 belong to 0 but no ovl1
				#output single belong to 0
				if($decide{"1"}){
					foreach my $scf_rep(sort {$decide{"1"}{$b}<=>$decide{"1"}{$a}} keys $decide{"1"}){
						if($rep_scf{$scf_rep} eq "1"){
							if($rep_scf2{$scf_rep} eq "1"){
								#both con
								my @aa=@{$rep_data1{$scf_rep}};
								my @bb=@{$rep_data2{$scf_rep}};
								if($aa[1] eq $bb[1]){print "$scf_rep\tsame end in contain;\n";}
								my $st=1;my $ed=$all_scf->{$scf_rep};my $typ="+";
								my $stat1=$aa[1];
								if($stat1 eq "ed"){$typ="-";}
								my $len1=$ed-$st+1;$agp_ed=$agp_st+$len1-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\td1_$scf_rep\t2\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
							if(!($rep_scf2{$scf_rep} eq 1)){
								#second have ovl
								$pre_ovl{$rep_scf2{$scf_rep}}=1;
								my @aa=@{$rep_data1{$scf_rep}};
								my @bb=@{$rep_data2{$scf_rep}};
								if($aa[1] eq $bb[1]){print "$scf_rep\tsame end in contain and ovl;\n";}
								#just ouput the ovl (2)
								my $scf_ovl=$rep_scf2{$scf_rep};
								my @cc=@{$rep_data2{$scf_ovl}};
								my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
								my $stat1=$bb[1];my $stat2=$cc[1];
								my $ovl1=$bb[3]-$bb[2]+1;my $ovl2=$cc[3]-$cc[2]+1;
								if($stat1 eq "st"){
									$st1=$bb[3]+1;$ed1=$all_scf->{$scf_rep};
									if($stat2 eq "ed"){
										$typ1="-";$typ2="-";
										$st2=1;$ed2=$cc[2]-1;
									}else{
										$typ1="-";$typ2="+";
										$st2=$cc[3]+1;$ed2=$all_scf->{$scf_ovl};
									}
								}else{
									$st1=1;$ed1=$bb[2]-1;
									if($stat2 eq "st"){
										$typ1="+";$typ2="+";
										$st2=$cc[3]+1;$ed2=$all_scf->{$scf_ovl};
									}else{
										$typ1="+";$typ2="-";
										$st2=1;$ed2=$cc[2]-1;
									}
								}
								if($ovl1>$ovl2){
									#$onetwo=1;
									$st1=1;$ed1=$all_scf->{$scf_rep};
									$agp_ed=$agp_st+$all_scf->{$scf_rep}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
									$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf_ovl};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf_ovl}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
							}
						}
						if($rep_scf{$scf_rep} ne "1"){
							if($rep_scf2{$scf_rep} eq "1"){
	                            #second no ovl
	                            my $scf_ovl=$rep_scf{$scf_rep};
	                            my @bb=@{$rep_data1{$scf_rep}};
	                            my @aa=@{$rep_data1{$scf_ovl}};
	                            my @cc=@{$rep_data2{$scf_rep}};
	                            if($bb[1] eq $cc[1]){print "$scf_rep\tsame end in ovl and contain;\n";}
	                            my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
	                            my $stat1=$aa[1];my $stat2=$bb[1];
	                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
	                            if($stat1 eq "st"){
	                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf_ovl};
	                                if($stat2 eq "ed"){
	                                    $typ1="-";$typ2="-";
	                                    $st2=1;$ed2=$bb[2]-1;
	                                }else{
	                                    $typ1="-";$typ2="+";
	                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf_rep};
	                                }
	                            }else{
	                                $st1=1;$ed1=$aa[2]-1;
	                                if($stat2 eq "st"){
	                                    $typ1="+";$typ2="+";
	                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf_rep};
	                                }else{
	                                    $typ1="+";$typ2="-";
	                                    $st2=1;$ed2=$bb[2]-1;
	                                }
	                            }
								if($ovl1>$ovl2){
									#$onetwo=1;
									$st1=1;$ed1=$all_scf->{$scf_ovl};
									$agp_ed=$agp_st+$all_scf->{$scf_ovl}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
									$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf_rep};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf_rep}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
	                        }
	                        if($rep_scf2{$scf_rep} ne "1"){
	                            #second have ovl && a little complicated situation
								my $scf_ovl1=$rep_scf{$scf_rep};
								my $scf_ovl2=$rep_scf2{$scf_rep};
								my @aa=@{$rep_data1{$scf_ovl1}};
								my @bb=@{$rep_data1{$scf_rep}};
								my @cc=@{$rep_data2{$scf_rep}};
								my @dd=@{$rep_data2{$scf_ovl2}};
								my $stat1=$aa[1];my $stat3=$dd[1];
								if($bb[1] eq $cc[1]){
									print "$scf_rep head and end stat is same and can not link\n";
									#output all record. the following pm will deal with this situation
									#output aa bb
									my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
		                            my $stat1=$aa[1];my $stat2=$bb[1];
		                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
		                            my $scf1=$aa[0];my $scf2=$bb[0];
		                            if($stat1 eq "st"){
		                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
		                                if($stat2 eq "ed"){
		                                    $typ2="-";$st2=1;$ed2=$bb[2]-1;
		                                }else{
		                                     $typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
		                                }
		                            }else{
		                                $st1=1;$ed1=$aa[2]-1;$typ1="+";
		                                if($stat2 eq "st"){
		                                    $typ2="+";
		                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
		                                }else{
		                                    $typ2="-";
		                                    $st2=1;$ed2=$bb[2]-1;
		                                }
		                            }
		                            if($ovl1>$ovl2){
		                                #$onetwo=1;
		                                $st1=1;$ed1=$all_scf->{$scf1};
		                                $agp_ed=$agp_st+$all_scf->{$scf1}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
		                                $agp_ed=$agp_st+$len2-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }else{
		                                #$onetwo=2;
		                                $st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
		                                $agp_ed=$agp_st+$len1-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+$all_scf->{$scf2}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }
		                            #output cc dd
		                            #my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
		                            $stat1=$cc[1];$stat2=$dd[1];
		                            $scf1=$cc[0];$scf2=$dd[0];
		                            $ovl1=$cc[3]-$cc[2]+1;$ovl2=$dd[3]-$dd[2]+1;
		                            if($stat1 eq "st"){
		                                $st1=$cc[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
		                                if($stat2 eq "ed"){
		                                    $typ2="-";$st2=1;$ed2=$dd[2]-1;
		                                }else{
		                                     $typ2="+";$st2=$dd[3]+1;$ed2=$all_scf->{$scf2};
		                                }
		                            }else{
		                                $st1=1;$ed1=$cc[2]-1;$typ1="+";
		                                if($stat2 eq "st"){
		                                    $typ2="+";
		                                    $st2=$dd[3]+1;$ed2=$all_scf->{$scf2};
		                                }else{
		                                    $typ2="-";
		                                    $st2=1;$ed2=$dd[2]-1;
		                                }
		                            }
		                            if($ovl1>$ovl2){
		                                #$onetwo=1;
		                                $st1=1;$ed1=$all_scf->{$scf1};
		                                $agp_ed=$agp_st+$all_scf->{$scf1}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
		                                $agp_ed=$agp_st+$len2-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }else{
		                                #$onetwo=2;
		                                $st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
		                                $agp_ed=$agp_st+$len1-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\td1_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+$all_scf->{$scf2}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\td1_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }
								}else{
									my ($st1,$ed1,$st2,$ed2,$st3,$ed3,$len1,$len2,$len3,$typ1,$typ2,$typ3);
									my $test_mid=0;
									my %tmp_mid;
									$tmp_mid{"1"}=$bb[2];$tmp_mid{"2"}=$bb[3];
									$tmp_mid{"3"}=$cc[2];$tmp_mid{"4"}=$cc[3];
									my @key_mid=sort {$tmp_mid{$a}<=>$tmp_mid{$b}} keys %tmp_mid;
									my $key_mid=join "",@key_mid;
									if($key_mid=~/^13/ || $key_mid=~/^31/){
										#check BAC is redundant;
										print "$scf_rep is redundant. need to re-blast $scf_ovl1 and $scf_ovl2. discard BAC $scf_rep\n";
										#deal with the overlap
										#aa bb
										my ($st21,$ed21,$st22,$ed22,$typ21,$typ22);
										my $stat1=$aa[1];my $stat2=$bb[1];
			                            my $scf1=$aa[0];my $scf2=$bb[0];
			                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
			                            if($stat1 eq "st"){
			                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
			                                if($stat2 eq "ed"){
			                                    $typ21="-";$st21=1;$ed21=$bb[2]-1;
			                                }else{
			                                     $typ21="+";$st21=$bb[3]+1;$ed21=$all_scf->{$scf2};
			                                }
			                            }else{
			                                $st1=1;$ed1=$aa[2]-1;$typ1="+";
			                                if($stat2 eq "st"){
			                                    $typ21="+";
			                                    $st21=$bb[3]+1;$ed21=$all_scf->{$scf2};
			                                }else{
			                                    $typ21="-";
			                                    $st21=1;$ed21=$bb[2]-1;
			                                }
			                            }
			                            if($ovl1>$ovl2){
			                                #$onetwo=1;
			                                $st1=1;$ed1=$all_scf->{$scf1};
			                            }else{
			                                #$onetwo=2;
			                                $st21=1;$ed21=$all_scf->{$scf2};
			                            }
			                            #cc dd
			                            $stat1=$cc[1];$stat2=$dd[1];
			                            $scf1=$cc[0];$scf2=$dd[0];
			                            $ovl1=$cc[3]-$cc[2]+1;$ovl2=$dd[3]-$dd[2]+1;
			                            if($stat1 eq "st"){
			                                $st22=$cc[3]+1;$ed22=$all_scf->{$scf1};$typ22="-";
			                                if($stat2 eq "ed"){
			                                    $typ3="-";$st3=1;$ed3=$dd[2]-1;
			                                }else{
			                                     $typ3="+";$st3=$dd[3]+1;$ed3=$all_scf->{$scf2};
			                                }
			                            }else{
			                                $st22=1;$ed22=$cc[2]-1;$typ22="+";
			                                if($stat2 eq "st"){
			                                    $typ3="+";$st3=$dd[3]+1;$ed3=$all_scf->{$scf2};
			                                }else{
			                                    $typ3="-";$st3=1;$ed3=$dd[2]-1;
			                                }
			                            }
			                            if($ovl1>$ovl2){
			                                #$onetwo=1;
			                                $st22=1;$ed22=$all_scf->{$scf1};
			                            }else{
			                                #$onetwo=2;
			                                $st3=1;$ed3=$all_scf->{$scf2};
			                            }
			                            #end
			                            #check bb cc coordonation
			                            my $check_ovl=check($st21,$ed21,$st22,$ed22);
			                            if($check_ovl){
											if($st21<$st22){$st2=$st22;}else{$st2=$st21;}
											if($ed21<$ed22){$ed2=$ed21;}else{$ed2=$ed22;}
											$typ2=$typ21;
											if($typ21 ne $typ22){
												$typ3=~tr/-+/+-/;
											}
											#output
											$len1=$ed1-$st1+1;$len2=$ed2-$st2+1;$len3=$ed3-$st3+1;
											$agp_ed=$agp_st+$len1-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len2-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len3-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td1_$scf_rep\t3\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											$agp_st=$agp_ed+1;$cal_agp++;
			                            }else{
			                            	#need trip lenght and need to discard scf_rep
			                            	my $p_len1=$st22-$ed21+1;my $p_len2=$st21-$ed22+1;
											my $len_trip;
											if($p_len1>$p_len2){$len_trip=$p_len1;}else{$len_trip=$p_len2;}
											if($typ1 eq "+"){
												$ed1=$ed1-$len_trip;
											}else{
												$st1=$st1+$len_trip;
											}
											if($typ21 ne $typ22){
												$typ3=~tr/-+/+-/;
											}
											#output
											$len1=$ed1-$st1+1;$len3=$ed3-$st3+1;
											$agp_ed=$agp_st+$len1-1;
											print "decide 1 $st1\t$ed1\t$st3\t$ed3\t$len1\t$len3\n";
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											#$agp_ed=$agp_st+$len2-1;
											#print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
											#print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\t2\n";
											#$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len3-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td1_$scf_rep\t3\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											$agp_st=$agp_ed+1;$cal_agp++;
			                            }
									}else{
										if($bb[3]<$cc[3]){
											$st2=$bb[3];$ed2=$cc[2];$test_mid=1;
										}else{
											$st2=$cc[3];$ed2=$bb[2];$test_mid=0;
										}
										if($stat1 eq "st"){
											if($test_mid){
												$typ1="-";$typ2="+";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}else{
												$typ1="-";$typ2="-";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}
										}else{
											if($test_mid){
												$typ1="+";$typ2="+";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}else{
												$typ1="+";$typ2="-";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}
										}
										my $ovl1=$aa[3]-$aa[2]+1;my $ovl21=$bb[3]-$bb[2]+1;
										my $ovl3=$dd[3]-$dd[2]+1;my $ovl22=$cc[3]-$cc[2]+1;
										if($ovl1>$ovl21){
											#keep 0
											$st1=1;$ed1=$all_scf->{$scf_ovl1};
										}else{
											#keek 11
											$st2=1;
											if($stat1 eq "st"){
												$st1=$aa[3]+1;$ed1=$all_scf->{$scf_ovl1};
											}else{
												$st1=1;$ed1=$aa[2]-1;
											}
										}
										if($ovl3>$ovl22){
											#keep 3
											$st3=1;$ed3=$all_scf->{$scf_ovl2};
										}else{
											#keep 21
											$ed2=$all_scf->{$scf_rep};
											if($stat3 eq "st"){
												$st3=$dd[3]+1;$ed3=$all_scf->{$scf_ovl2};
											}else{
												$st3=1;$ed3=$dd[2]-1;
											}
										}
										$len1=$ed1-$st1+1;$len2=$ed2-$st2+1;$len3=$ed3-$st3+1;
										$agp_ed=$agp_st+$len1-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td1_$scf_rep\t1\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+$len2-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td1_$scf_rep\t2\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+$len3-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td1_$scf_rep\t3\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
									}
								}
	                        }
						}
					}
				}
				#output single belong to 1
				if($decide{"2"}){
					print "decide 2\n";
					foreach my $scf_rep(sort {$decide{"2"}{$a}<=>$decide{"2"}{$b}} keys $decide{"2"}){
						if($rep_scf{$scf_rep} eq "1"){
							if($rep_scf2{$scf_rep} eq "1"){
								#both con
								my @aa=@{$rep_data1{$scf_rep}};
								my @bb=@{$rep_data2{$scf_rep}};
								if($aa[1] eq $bb[1]){print "$scf_rep\tsame end in contain;\n";}
								my $st=1;my $ed=$all_scf->{$scf_rep};my $typ="+";
								my $stat1=$aa[1];
								if($stat1 eq "ed"){$typ="-";}
								my $len1=$ed-$st+1;$agp_ed=$agp_st+$len1-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\td2_$scf_rep\t1\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
							if(!($rep_scf2{$scf_rep} eq 1)){
								#second have ovl
								$pre_ovl{$rep_scf2{$scf_rep}}=1;
								my @aa=@{$rep_data1{$scf_rep}};
								my @bb=@{$rep_data2{$scf_rep}};
								if($aa[1] eq $bb[1]){print "$scf_rep\tsame end in contain and ovl;\n";}
								#just ouput the ovl (2)
								my $scf_ovl=$rep_scf2{$scf_rep};
								my @cc=@{$rep_data2{$scf_ovl}};
								my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
								my $stat1=$bb[1];my $stat2=$cc[1];
								my $ovl1=$bb[3]-$bb[2]+1;my $ovl2=$cc[3]-$cc[2]+1;
								if($stat1 eq "st"){
									$st1=$bb[3]+1;$ed1=$all_scf->{$scf_rep};
									if($stat2 eq "ed"){
										$typ1="-";$typ2="-";
										$st2=1;$ed2=$cc[2]-1;
									}else{
										$typ1="-";$typ2="+";
										$st2=$cc[3]+1;$ed2=$all_scf->{$scf_ovl};
									}
								}else{
									$st1=1;$ed1=$bb[2]-1;
									if($stat2 eq "st"){
										$typ1="+";$typ2="+";
										$st2=$cc[3]+1;$ed2=$all_scf->{$scf_ovl};
									}else{
										$typ1="+";$typ2="-";
										$st2=1;$ed2=$cc[2]-1;
									}
								}
								if($ovl1>$ovl2){
									#$onetwo=1;
									$st1=1;$ed1=$all_scf->{$scf_rep};
									$agp_ed=$agp_st+$all_scf->{$scf_rep}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
									$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf_ovl};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf_ovl}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
							}
						}
						if($rep_scf{$scf_rep} ne "1"){
							if($rep_scf2{$scf_rep} eq "1"){
	                            #second no ovl
	                            my $scf_ovl=$rep_scf{$scf_rep};
	                            my @bb=@{$rep_data1{$scf_rep}};
	                            my @aa=@{$rep_data1{$scf_ovl}};
	                            my @cc=@{$rep_data2{$scf_rep}};
	                            if($bb[1] eq $cc[1]){print "$scf_rep\tsame end in ovl and contain;\n";}
	                            my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
	                            my $stat1=$aa[1];my $stat2=$bb[1];
	                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
	                            if($stat1 eq "st"){
	                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf_ovl};
	                                if($stat2 eq "ed"){
	                                    $typ1="-";$typ2="-";
	                                    $st2=1;$ed2=$bb[2]-1;
	                                }else{
	                                    $typ1="-";$typ2="+";
	                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf_rep};
	                                }
	                            }else{
	                                $st1=1;$ed1=$aa[2]-1;
	                                if($stat2 eq "st"){
	                                    $typ1="+";$typ2="+";
	                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf_rep};
	                                }else{
	                                    $typ1="+";$typ2="-";
	                                    $st2=1;$ed2=$bb[2]-1;
	                                }
	                            }
								if($ovl1>$ovl2){
									#$onetwo=1;
									$st1=1;$ed1=$all_scf->{$scf_ovl};
									$agp_ed=$agp_st+$all_scf->{$scf_ovl}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
									$agp_ed=$agp_st+$len2-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf_rep};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf_rep}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
	                        }
	                        if($rep_scf2{$scf_rep} ne "1"){
	                            #second have ovl && a little complicated situation
								my $scf_ovl1=$rep_scf{$scf_rep};
								my $scf_ovl2=$rep_scf2{$scf_rep};
								my @aa=@{$rep_data1{$scf_ovl1}};
								my @bb=@{$rep_data1{$scf_rep}};
								my @cc=@{$rep_data2{$scf_rep}};
								my @dd=@{$rep_data2{$scf_ovl2}};
								my $stat1=$aa[1];my $stat3=$dd[1];
								if($bb[1] eq $cc[1]){
									print "$scf_rep head and end stat is same and can not link\n";
									#output all record. the following pm will deal with this situation
									#output aa bb
									my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
		                            my $stat1=$aa[1];my $stat2=$bb[1];
		                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
		                            my $scf1=$aa[0];my $scf2=$bb[0];
		                            if($stat1 eq "st"){
		                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
		                                if($stat2 eq "ed"){
		                                    $typ2="-";$st2=1;$ed2=$bb[2]-1;
		                                }else{
		                                     $typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
		                                }
		                            }else{
		                                $st1=1;$ed1=$aa[2]-1;$typ1="+";
		                                if($stat2 eq "st"){
		                                    $typ2="+";
		                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
		                                }else{
		                                    $typ2="-";
		                                    $st2=1;$ed2=$bb[2]-1;
		                                }
		                            }
		                            if($ovl1>$ovl2){
		                                #$onetwo=1;
		                                $st1=1;$ed1=$all_scf->{$scf1};
		                                $agp_ed=$agp_st+$all_scf->{$scf1}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
		                                $agp_ed=$agp_st+$len2-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }else{
		                                #$onetwo=2;
		                                $st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
		                                $agp_ed=$agp_st+$len1-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+$all_scf->{$scf2}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }
		                            #output cc dd
		                            #my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
		                            $stat1=$cc[1];$stat2=$dd[1];
		                            $scf1=$cc[0];$scf2=$dd[0];
		                            $ovl1=$cc[3]-$cc[2]+1;$ovl2=$dd[3]-$dd[2]+1;
		                            if($stat1 eq "st"){
		                                $st1=$cc[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
		                                if($stat2 eq "ed"){
		                                    $typ2="-";$st2=1;$ed2=$dd[2]-1;
		                                }else{
		                                     $typ2="+";$st2=$dd[3]+1;$ed2=$all_scf->{$scf2};
		                                }
		                            }else{
		                                $st1=1;$ed1=$cc[2]-1;$typ1="+";
		                                if($stat2 eq "st"){
		                                    $typ2="+";
		                                    $st2=$dd[3]+1;$ed2=$all_scf->{$scf2};
		                                }else{
		                                    $typ2="-";
		                                    $st2=1;$ed2=$dd[2]-1;
		                                }
		                            }
		                            if($ovl1>$ovl2){
		                                #$onetwo=1;
		                                $st1=1;$ed1=$all_scf->{$scf1};
		                                $agp_ed=$agp_st+$all_scf->{$scf1}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
		                                $agp_ed=$agp_st+$len2-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }else{
		                                #$onetwo=2;
		                                $st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
		                                $agp_ed=$agp_st+$len1-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$cc[0]\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+$all_scf->{$scf2}-1;
		                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$dd[0]\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
		                                $cal_agp++;$agp_st=$agp_ed+1;
		                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
		                            }
								}else{
									my ($st1,$ed1,$st2,$ed2,$st3,$ed3,$len1,$len2,$len3,$typ1,$typ2,$typ3);
									my $test_mid=0;
									my %tmp_mid;
									$tmp_mid{"1"}=$bb[2];$tmp_mid{"2"}=$bb[3];
									$tmp_mid{"3"}=$cc[2];$tmp_mid{"4"}=$cc[3];
									my @key_mid=sort {$tmp_mid{$a}<=>$tmp_mid{$b}} keys %tmp_mid;
									my $key_mid=join "",@key_mid;
									if($key_mid=~/^13/ || $key_mid=~/^31/){
										#check BAC is redundant;
										print "$scf_rep is redundant. need to re-blast $scf_ovl1 and $scf_ovl2. discard BAC $scf_rep\n";
										#deal with the overlap
										#aa bb
										my ($st21,$ed21,$st22,$ed22,$typ21,$typ22);
										my $stat1=$aa[1];my $stat2=$bb[1];
			                            my $scf1=$aa[0];my $scf2=$bb[0];
			                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
			                            if($stat1 eq "st"){
			                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
			                                if($stat2 eq "ed"){
			                                    $typ21="-";$st21=1;$ed21=$bb[2]-1;
			                                }else{
			                                     $typ21="+";$st21=$bb[3]+1;$ed21=$all_scf->{$scf2};
			                                }
			                            }else{
			                                $st1=1;$ed1=$aa[2]-1;$typ1="+";
			                                if($stat2 eq "st"){
			                                    $typ21="+";
			                                    $st21=$bb[3]+1;$ed21=$all_scf->{$scf2};
			                                }else{
			                                    $typ21="-";
			                                    $st21=1;$ed21=$bb[2]-1;
			                                }
			                            }
			                            if($ovl1>$ovl2){
			                                #$onetwo=1;
			                                $st1=1;$ed1=$all_scf->{$scf1};
			                            }else{
			                                #$onetwo=2;
			                                $st21=1;$ed21=$all_scf->{$scf2};
			                            }
			                            #cc dd
			                            $stat1=$cc[1];$stat2=$dd[1];
			                            $scf1=$cc[0];$scf2=$dd[0];
			                            $ovl1=$cc[3]-$cc[2]+1;$ovl2=$dd[3]-$dd[2]+1;
			                            if($stat1 eq "st"){
			                                $st22=$cc[3]+1;$ed22=$all_scf->{$scf1};$typ22="-";
			                                if($stat2 eq "ed"){
			                                    $typ3="-";$st3=1;$ed3=$dd[2]-1;
			                                }else{
			                                     $typ3="+";$st3=$dd[3]+1;$ed3=$all_scf->{$scf2};
			                                }
			                            }else{
			                                $st22=1;$ed22=$cc[2]-1;$typ22="+";
			                                if($stat2 eq "st"){
			                                    $typ3="+";$st3=$dd[3]+1;$ed3=$all_scf->{$scf2};
			                                }else{
			                                    $typ3="-";$st3=1;$ed3=$dd[2]-1;
			                                }
			                            }
			                            if($ovl1>$ovl2){
			                                #$onetwo=1;
			                                $st22=1;$ed22=$all_scf->{$scf1};
			                            }else{
			                                #$onetwo=2;
			                                $st3=1;$ed3=$all_scf->{$scf2};
			                            }
			                            #end
			                            #check bb cc coordonation
			                            my $check_ovl=check($st21,$ed21,$st22,$ed22);
			                            if($check_ovl){
											if($st21<$st22){$st2=$st22;}else{$st2=$st21;}
											if($ed21<$ed22){$ed2=$ed21;}else{$ed2=$ed22;}
											$typ2=$typ21;
											if($typ21 ne $typ22){
												$typ3=~tr/-+/+-/;
											}
											#output
											$len1=$ed1-$st1+1;$len2=$ed2-$st2+1;$len3=$ed3-$st3+1;
											$agp_ed=$agp_st+$len1-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len2-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len3-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td2_$scf_rep\t3\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											$agp_st=$agp_ed+1;$cal_agp++;
			                            }else{
			                            	#need trip lenght and need to discard scf_rep
			                            	my $p_len1=$st22-$ed21+1;my $p_len2=$st21-$ed22+1;
											my $len_trip;
											if($p_len1>$p_len2){$len_trip=$p_len1;}else{$len_trip=$p_len2;}
											if($typ1 eq "+"){
												$ed1=$ed1-$len_trip;
											}else{
												$st1=$st1+$len_trip;
											}
											if($typ21 ne $typ22){
												$typ3=~tr/-+/+-/;
											}
											#output
											$len1=$ed1-$st1+1;
											$len3=$ed3-$st3+1;
											$agp_ed=$agp_st+$len1-1;
											print "decide 2 $st1\t$ed1\t$st3\t$ed3\t$len1\t$len3\n";
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											#$agp_ed=$agp_st+$len2-1;
											#print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
											#print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\t2\n";
											#$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+$len3-1;
											print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td2_$scf_rep\t3\n";
											$cal_agp++;$agp_st=$agp_ed+1;
											$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
											$agp_st=$agp_ed+1;$cal_agp++;
			                            }
									}else{
										if($bb[3]<$cc[3]){
											$st2=$bb[3];$ed2=$cc[2];$test_mid=1;
										}else{
											$st2=$cc[3];$ed2=$bb[2];$test_mid=0;
										}
										if($stat1 eq "st"){
											if($test_mid){
												$typ1="-";$typ2="+";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}else{
												$typ1="-";$typ2="-";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}
										}else{
											if($test_mid){
												$typ1="+";$typ2="+";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}else{
												$typ1="+";$typ2="-";
												if($stat3 eq "st"){$typ3="+";}else{$typ3="-";}
											}
										}
										my $ovl1=$aa[3]-$aa[2]+1;my $ovl21=$bb[3]-$bb[2]+1;
										my $ovl3=$dd[3]-$dd[2]+1;my $ovl22=$cc[3]-$cc[2]+1;
										if($ovl1>$ovl21){
											#keep 0
											$st1=1;$ed1=$all_scf->{$scf_ovl1};
										}else{
											#keek 11
											$st2=1;
											if($stat1 eq "st"){
												$st1=$aa[3]+1;$ed1=$all_scf->{$scf_ovl1};
											}else{
												$st1=1;$ed1=$aa[2]-1;
											}
										}
										if($ovl3>$ovl22){
											#keep 3
											$st3=1;$ed3=$all_scf->{$scf_ovl2};
										}else{
											#keep 21
											$ed2=$all_scf->{$scf_rep};
											if($stat3 eq "st"){
												$st3=$dd[3]+1;$ed3=$all_scf->{$scf_ovl2};
											}else{
												$st3=1;$ed3=$dd[2]-1;
											}
										}
										$len1=$ed1-$st1+1;$len2=$ed2-$st2+1;$len3=$ed3-$st3+1;
										$agp_ed=$agp_st+$len1-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl1\t$st1\t$ed1\t$typ1\td2_$scf_rep\t1\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+$len2-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_rep\t$st2\t$ed2\t$typ2\td2_$scf_rep\t2\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+$len3-1;
										print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$scf_ovl2\t$st3\t$ed3\t$typ3\td2_$scf_rep\t3\n";
										$cal_agp++;$agp_st=$agp_ed+1;
										$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
										$agp_st=$agp_ed+1;$cal_agp++;
									}
								}
	                        }
						}
					}
				}
			}else{
				#no repeat scaffold
				#output the non-repeat sort;
				print "nonrepeat\n";
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$first}){
					#check if pre 0 contain repeat
					if($out->{$first}{$ctg_sort}[0][0] && @{$out->{$first}{$ctg_sort}[0]}>1){
						if($repeat{$out->{$first}{$ctg_sort}[0][0]}){
							print "$first\t$out->{$first}{$ctg_sort}[0][0]\tpre 0 have repeat\n";next;
						}
					}
					if($out->{$first}{$ctg_sort}[1][0] && @{$out->{$first}{$ctg_sort}[1]}>1){
						if($repeat{$out->{$first}{$ctg_sort}[1][0]}){print "$first\t$out->{$first}{$ctg_sort}[1][0]\tpre 1 have repeat\n";next;}
					}
					#output other
					if($out->{$first}{$ctg_sort}[0][0] && @{$out->{$first}{$ctg_sort}[0]}>1){
						if($out->{$first}{$ctg_sort}[1][0] && @{$out->{$first}{$ctg_sort}[1]}>1){
							#ovl relation
							my $scf1=$out->{$first}{$ctg_sort}[0][0];my $scf2=$out->{$first}{$ctg_sort}[1][0];
                            my @aa=@{$out->{$first}{$ctg_sort}[0]};
                            my @bb=@{$out->{$first}{$ctg_sort}[1]};
                            my ($st1,$ed1,$typ1,$st2,$ed2,$typ2,$len1,$len2);
                            my $stat1=$aa[1];my $stat2=$bb[1];
                            my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
                            if($stat1 eq "st"){
                                $st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
                                if($stat2 eq "ed"){
                                    $typ2="-";$st2=1;$ed2=$bb[2]-1;
                                }else{
                                     $typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                                }
                            }else{
                                $st1=1;$ed1=$aa[2]-1;$typ1="+";
                                if($stat2 eq "st"){
                                    $typ2="+";
                                    $st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                                }else{
                                    $typ2="-";
                                    $st2=1;$ed2=$bb[2]-1;
                                }
                            }
                            if($ovl1>$ovl2){
                                #$onetwo=1;
                                $st1=1;$ed1=$all_scf->{$scf1};
                                $agp_ed=$agp_st+$all_scf->{$scf1}-1;
                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                                $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
                                $agp_ed=$agp_st+$len2-1;
                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                                $cal_agp++;$agp_st=$agp_ed+1;
                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                            }else{
                                #$onetwo=2;
                                $st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
                                $agp_ed=$agp_st+$len1-1;
                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                                $cal_agp++;$agp_st=$agp_ed+1;
                                $agp_ed=$agp_st+$all_scf->{$scf2}-1;
                                print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                                $cal_agp++;$agp_st=$agp_ed+1;
                                $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                            }
						}else{
							#only 0
							my $scf1=$out->{$first}{$ctg_sort}[0][0];
							my @aa=@{$out->{$first}{$ctg_sort}[0]};
							my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
							my $stat1=$out->{$first}{$ctg_sort}[0][1];
							if($stat1 eq "st"){$typ="-";}
							my $len1=$ed-$st+1;$agp_ed=$agp_st+$len1-1;
							print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
							$cal_agp++;$agp_st=$agp_ed+1;
							$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
						}
					}else{
						#only 1
						if($out->{$first}{$ctg_sort}[1][0] && @{$out->{$first}{$ctg_sort}[1]}>1){
							my $scf2=$out->{$first}{$ctg_sort}[1][0];
							my @bb=@{$out->{$first}{$ctg_sort}[1]};
							my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
							my $stat2=$out->{$first}{$ctg_sort}[1][1];
							if($stat2 eq "ed"){$typ="-";}
							my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
							print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
							$cal_agp++;$agp_st=$agp_ed+1;
							$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
						}
					}
				}
				#make second ready to output if second-first>1 in the next round
				
			}
			#check the if in last situation and output the left record
			if($last){
				#output the left recond.
				#they are only the second sort and the second BAC to be output
				#just result{"2"} if in repeat condition and whole record in non-repeat condition. 
				#the non-neighbor will be just the whole record
				print "last check result 2\n";
				my @need_key;
				if($result{"2"}){
					#first check the sort is correct or not
					my @key=sort {$a<=>$b} keys $result{"2"};
					foreach my $need_key(@key){
						if($result{"2"}{$need_key}[2] == $second){
							push @need_key,$need_key;
						}
					}
				}
				if(@need_key>0){
					print "in result 2\n";
					foreach my $key(sort {$a<=>$b} @need_key){
						if($result{"2"}{$key}[0][0] && @{$result{"2"}{$key}[0]}>1){
							if($result{"2"}{$key}[1][0] && @{$result{"2"}{$key}[1]}>1){
						    	#ovl check st ed
					        	my $stat1=$result{"2"}{$key}[0][1];
					        	my $stat2=$result{"2"}{$key}[1][1];
								my @aa=@{$result{"2"}{$key}[0]};my $scf1=$aa[0];
								my @bb=@{$result{"2"}{$key}[1]};my $scf2=$bb[0];
								my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
					        	my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
								if($stat1 eq "st"){
									$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
							    	if($stat2 eq "ed"){
									    $typ2="-";$st2=1;$ed2=$bb[2]-1;
							   		}else{
							   	    	$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
							    	}
								}else{
									$st1=1;$ed1=$aa[2]-1;$typ1="+";
							    	if($stat2 eq "st"){
							    		$typ2="+";
							    		$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
							    	}else{
							   	 		$typ2="-";
							    		$st2=1;$ed2=$bb[2]-1;
							    	}
								}
								if($ovl1>$ovl2){
								    #$onetwo=1;
								    $st1=1;$ed1=$all_scf->{$scf1};
								    $agp_ed=$agp_st+$all_scf->{$scf1}-1;
								    print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
								    $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
								    $agp_ed=$agp_st+$len2-1;
								    print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
								    $cal_agp++;$agp_st=$agp_ed+1;
								    $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf2}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
							}else{
								#only 0
								my @aa=@{$result{"2"}{$key}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
								my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
								if($stat1 eq "st"){$typ="-";}
								my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
						}else{
							#check if only 1
							if($result{"2"}{$key}[1][0] && @{$result{"2"}{$key}[1]}>1){
								my @bb=@{$result{"2"}{$key}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
								my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
								if($stat2 eq "ed"){$typ="-";}
								my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
						}
					}	
				}else{
					#in non-repeat case. output the whole record
					print "in whole record\n";
					my $aft_s=$second;
					foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$aft_s}){
						if($out->{$aft_s}{$ctg_sort}[0][0] && @{$out->{$aft_s}{$ctg_sort}[0]}>1){
							if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
								#ovl check st ed
								my $stat1=$out->{$aft_s}{$ctg_sort}[0][1];
	                         	my $stat2=$out->{$aft_s}{$ctg_sort}[1][1];
	                         	my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];
	                         	my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];
	                         	my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
	                         	my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
	                         	if($stat1 eq "st"){
	                            	$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
	                            	if($stat2 eq "ed"){
	                                	$typ2="-";$st2=1;$ed2=$bb[2]-1;
	                             	}else{
	                                	$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
	                             	}
	                         	}else{
	                             	$st1=1;$ed1=$aa[2]-1;$typ1="+";
	                             	if($stat2 eq "st"){
	                                 	$typ2="+";
	                                 	$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
	                             	}else{
	                                 	$typ2="-";
	                                 	$st2=1;$ed2=$bb[2]-1;
	                             	}
	                         	}
	                         	if($ovl1>$ovl2){
	                             	#$onetwo=1;
	                             	$st1=1;$ed1=$all_scf->{$scf1};
	                             	$agp_ed=$agp_st+$all_scf->{$scf1}-1;
	                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
	                             	$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
	                             	$agp_ed=$agp_st+$len2-1;
	                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
	                             	$cal_agp++;$agp_st=$agp_ed+1;
	                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
	                         	}else{
	                             	#$onetwo=2;
	                            	$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
	                             	$agp_ed=$agp_st+$len1-1;
	                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
	                             	$cal_agp++;$agp_st=$agp_ed+1;
	                             	$agp_ed=$agp_st+$all_scf->{$scf2}-1;
	                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
	                             	$cal_agp++;$agp_st=$agp_ed+1;
	                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";	
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
	                         	}
	                     	}else{
	                         	#only 0
	                         	my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
	                         	my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
	                         	if($stat1 eq "st"){$typ="-";}
	                         	my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
	                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
	                         	$cal_agp++;$agp_st=$agp_ed+1;
	                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
	                     	}
	                 	}else{
	                     	#check if only 1
	                     	if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
	                        	my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
	                         	my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
	                         	if($stat2 eq "ed"){$typ="-";}
	                         	my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
	                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
	                         	$cal_agp++;$agp_st=$agp_ed+1;
	                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
	                     	}
                 		}
					}
				}
				last;
			}
		}else{
			#not neighbor BAC relation. need to output all first non-repeat record. 
			#but if the second is the last sort them we will need to
			#check the sort and output all the last second sort
			#one special case: the first is the start point;
			my $pre_s=$first;my $aft_s=$second;
			my $bac1=$rel_ctg->{$pre_s}[0];
			my $bac2=$rel_ctg->{$aft_s}[0];
			print "$bac1\t$bac2\tnot neighbor result\n";
			#first ==1; output whole record
			if($first==$key_sort[0]){
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$pre_s}){
					if($out->{$pre_s}{$ctg_sort}[0][0] && @{$out->{$pre_s}{$ctg_sort}[0]}>1){
						if($out->{$pre_s}{$ctg_sort}[1][0] && @{$out->{$pre_s}{$ctg_sort}[1]}>1){
							#ovl check st ed
							my $stat1=$out->{$pre_s}{$ctg_sort}[0][1];
                         	my $stat2=$out->{$pre_s}{$ctg_sort}[1][1];
                         	my @aa=@{$out->{$pre_s}{$ctg_sort}[0]};my $scf1=$aa[0];
                         	my @bb=@{$out->{$pre_s}{$ctg_sort}[1]};my $scf2=$bb[0];
                         	my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
                         	my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
                         	if($stat1 eq "st"){
                            	$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
                            	if($stat2 eq "ed"){
                                	$typ2="-";$st2=1;$ed2=$bb[2]-1;
                             	}else{
                                	$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}
                         	}else{
                             	$st1=1;$ed1=$aa[2]-1;$typ1="+";
                             	if($stat2 eq "st"){
                                 	$typ2="+";
                                 	$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}else{
                                 	$typ2="-";
                                 	$st2=1;$ed2=$bb[2]-1;
                             	}
                         	}
                         	if($ovl1>$ovl2){
                             	#$onetwo=1;
                             	$st1=1;$ed1=$all_scf->{$scf1};
                             	$agp_ed=$agp_st+$all_scf->{$scf1}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
                             	$agp_ed=$agp_st+$len2-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}else{
                             	#$onetwo=2;
                            	$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
                             	$agp_ed=$agp_st+$len1-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+$all_scf->{$scf2}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";	
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}
                     	}else{
                         	#only 0
                         	my @aa=@{$out->{$pre_s}{$ctg_sort}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
                         	my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
                         	if($stat1 eq "st"){$typ="-";}
                         	my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}else{
                     	#check if only 1
                     	if($out->{$pre_s}{$ctg_sort}[1][0] && @{$out->{$pre_s}{$ctg_sort}[1]}>1){
                        	my @bb=@{$out->{$pre_s}{$ctg_sort}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
                         	my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
                         	if($stat2 eq "ed"){$typ="-";}
                         	my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}
				}
			}
			#print the end of the  pre. we have two condition:pre have repeat and pre didn't have repeat.
			#pre have repeat . output the result{"2"}
			print "check result 2\n";
			my @need_key;
			if($result{"2"}){
				#first check the sort is correct or not
				my @key=sort {$a<=>$b} keys $result{"2"};
				foreach my $need_key(@key){
					if($result{"2"}{$need_key}[2] == $first){
						push @need_key,$need_key;
						print "$need_key\t$first\n";
					}
				}
				
				if(@need_key>0){
					print "have result2\t$first\n";
					foreach my $key(sort {$a<=>$b} @need_key){
						print "key\t$key\n";
						if($result{"2"}{$key}[0][0] && @{$result{"2"}{$key}[0]}>1){
							if($result{"2"}{$key}[1][0] && @{$result{"2"}{$key}[1]}>1){
							    #ovl check st ed
						        my $stat1=$result{"2"}{$key}[0][1];
						        my $stat2=$result{"2"}{$key}[1][1];
								my @aa=@{$result{"2"}{$key}[0]};my $scf1=$aa[0];
								my @bb=@{$result{"2"}{$key}[1]};my $scf2=$bb[0];
								print "$scf1\t$scf2\tresult 2\t$first\n";
								my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
						        my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
								if($stat1 eq "st"){
									$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
								    if($stat2 eq "ed"){
									    $typ2="-";$st2=1;$ed2=$bb[2]-1;
								    }else{
								        $typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
								    }
								}else{
									$st1=1;$ed1=$aa[2]-1;$typ1="+";
								    if($stat2 eq "st"){
								    	$typ2="+";
								    	$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
								    }else{
								    	$typ2="-";
								    	$st2=1;$ed2=$bb[2]-1;
								    }
								}
								if($ovl1>$ovl2){
								    #$onetwo=1;
								    $st1=1;$ed1=$all_scf->{$scf1};
								    $agp_ed=$agp_st+$all_scf->{$scf1}-1;
								    print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
								    $cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
								    $agp_ed=$agp_st+$len2-1;
								    print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
								    $cal_agp++;$agp_st=$agp_ed+1;
								    $agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}else{
									#$onetwo=2;
									$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
									$agp_ed=$agp_st+$len1-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\t1\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+$all_scf->{$scf2}-1;
									print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
									$cal_agp++;$agp_st=$agp_ed+1;
									$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
									$agp_st=$agp_ed+1;$cal_agp++;
								}
							}else{
								#only 0
								my @aa=@{$result{"2"}{$key}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
								my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
								if($stat1 eq "st"){$typ="-";}
								my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
						}else{
							#check if only 1
							if($result{"2"}{$key}[1][0] && @{$result{"2"}{$key}[1]}>1){
								my @bb=@{$result{"2"}{$key}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
								my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
								if($stat2 eq "ed"){$typ="-";}
								my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
								print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
								$cal_agp++;$agp_st=$agp_ed+1;
								$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
							}
						}
					}
				}
				$result{"2"}=();
			}
			#pre didn't have repeat. output the whole record
			if(@need_key<1 && $result_pre2[0] && $result_pre2[1] == $first){
				#$pre_s is the sort
				print "whole record in nonrepeat result_pre2\t$first\t$result_pre2[1]\n";
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$pre_s}){
					if($out->{$pre_s}{$ctg_sort}[0][0] && @{$out->{$pre_s}{$ctg_sort}[0]}>1){
						if($out->{$pre_s}{$ctg_sort}[1][0] && @{$out->{$pre_s}{$ctg_sort}[1]}>1){
							#ovl check st ed
							my $stat1=$out->{$pre_s}{$ctg_sort}[0][1];
                         	my $stat2=$out->{$pre_s}{$ctg_sort}[1][1];
                         	my @aa=@{$out->{$pre_s}{$ctg_sort}[0]};my $scf1=$aa[0];
                         	my @bb=@{$out->{$pre_s}{$ctg_sort}[1]};my $scf2=$bb[0];
                         	my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
                         	my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
                         	if($stat1 eq "st"){
                            	$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
                            	if($stat2 eq "ed"){
                                	$typ2="-";$st2=1;$ed2=$bb[2]-1;
                             	}else{
                                	$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}
                         	}else{
                             	$st1=1;$ed1=$aa[2]-1;$typ1="+";
                             	if($stat2 eq "st"){
                                 	$typ2="+";
                                 	$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}else{
                                 	$typ2="-";
                                 	$st2=1;$ed2=$bb[2]-1;
                             	}
                         	}
                         	if($ovl1>$ovl2){
                             	#$onetwo=1;
                             	$st1=1;$ed1=$all_scf->{$scf1};
                             	$agp_ed=$agp_st+$all_scf->{$scf1}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
                             	$agp_ed=$agp_st+$len2-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}else{
                             	#$onetwo=2;
                            	$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
                             	$agp_ed=$agp_st+$len1-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+$all_scf->{$scf2}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";	
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}
                     	}else{
                         	#only 0
                         	my @aa=@{$out->{$pre_s}{$ctg_sort}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
                         	my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
                         	if($stat1 eq "st"){$typ="-";}
                         	my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}else{
                     	#check if only 1
                     	if($out->{$pre_s}{$ctg_sort}[1][0] && @{$out->{$pre_s}{$ctg_sort}[1]}>1){
                        	my @bb=@{$out->{$pre_s}{$ctg_sort}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
                         	my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
                         	if($stat2 eq "ed"){$typ="-";}
                         	my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}
				}
				@result_pre2=();
			}

			#output 30Kbp ns to represent the miss BAC and can be insert bac if the BACs is exists but not mapped
			$agp_ed=$agp_st+30000-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tn\t30000\tscaffold\tno\tna\n";
			print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tn\t30000\tscaffold\tno\tna\n";
			$agp_st=$agp_ed+1;$cal_agp++;

			if($last){
				#output all second sort record
				foreach my $ctg_sort(sort {$a<=>$b} keys $out->{$aft_s}){
					if($out->{$aft_s}{$ctg_sort}[0][0] && @{$out->{$aft_s}{$ctg_sort}[0]}>1){
						if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
							#ovl check st ed
							my $stat1=$out->{$aft_s}{$ctg_sort}[0][1];
                         	my $stat2=$out->{$aft_s}{$ctg_sort}[1][1];
                         	my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];
                         	my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];
                         	my ($st1,$ed1,$st2,$ed2,$len1,$len2,$typ1,$typ2);
                         	my $ovl1=$aa[3]-$aa[2]+1;my $ovl2=$bb[3]-$bb[2]+1;
                         	if($stat1 eq "st"){
                            	$st1=$aa[3]+1;$ed1=$all_scf->{$scf1};$typ1="-";
                            	if($stat2 eq "ed"){
                                	$typ2="-";$st2=1;$ed2=$bb[2]-1;
                             	}else{
                                	$typ2="+";$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}
                         	}else{
                             	$st1=1;$ed1=$aa[2]-1;$typ1="+";
                             	if($stat2 eq "st"){
                                 	$typ2="+";
                                 	$st2=$bb[3]+1;$ed2=$all_scf->{$scf2};
                             	}else{
                                 	$typ2="-";
                                 	$st2=1;$ed2=$bb[2]-1;
                             	}
                         	}
                         	if($ovl1>$ovl2){
                             	#$onetwo=1;
                             	$st1=1;$ed1=$all_scf->{$scf1};
                             	$agp_ed=$agp_st+$all_scf->{$scf1}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;$len2=$ed2-$st2+1;
                             	$agp_ed=$agp_st+$len2-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}else{
                             	#$onetwo=2;
                            	$st2=1;$ed2=$all_scf->{$scf2};$len1=$ed1-$st1+1;
                             	$agp_ed=$agp_st+$len1-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st1\t$ed1\t$typ1\tnod\t1\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+$all_scf->{$scf2}-1;
                             	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\n";
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st2\t$ed2\t$typ2\tnod\t2\n";
                             	$cal_agp++;$agp_st=$agp_ed+1;
                             	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";	
								print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
								$agp_st=$agp_ed+1;$cal_agp++;
                         	}
                     	}else{
                         	#only 0
                         	my @aa=@{$out->{$aft_s}{$ctg_sort}[0]};my $scf1=$aa[0];my $stat1=$aa[1];
                         	my $st=1;my $ed=$all_scf->{$scf1};my $typ="+";
                         	if($stat1 eq "st"){$typ="-";}
                         	my $len=$ed-$st+1;$agp_ed=$agp_st+$len-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$aa[0]\t$st\t$ed\t$typ\tnod\t1\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}else{
                     	#check if only 1
                     	if($out->{$aft_s}{$ctg_sort}[1][0] && @{$out->{$aft_s}{$ctg_sort}[1]}>1){
                        	my @bb=@{$out->{$aft_s}{$ctg_sort}[1]};my $scf2=$bb[0];my $stat2=$bb[1];
                         	my $st=1;my $ed=$all_scf->{$scf2};my $typ="+";
                         	if($stat2 eq "ed"){$typ="-";}
                         	my $len2=$ed-$st+1;$agp_ed=$agp_st+$len2-1;
                         	print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tW\t$bb[0]\t$st\t$ed\t$typ\tnod\t2\n";
                         	$cal_agp++;$agp_st=$agp_ed+1;
                         	$agp_ed=$agp_st+100-1;print "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							print OUT "$ctg\t$agp_st\t$agp_ed\t$cal_agp\tN\t100\tscaffold\tno\tna\n";
							$agp_st=$agp_ed+1;$cal_agp++;
                     	}
                 	}
				}
				print "last check sort\n";
				last;
			}
		}
		#%result_pre2;
		@result_pre2=(1,$second);
	}
	#$result{"2"}=();
	close OUT;
}

sub check{
	my ($st1,$ed1,$st2,$ed2)=@_;
	my %kk;
	$kk{"1"}=$st1;$kk{"2"}=$ed1;
	$kk{"3"}=$st2;$kk{"4"}=$ed2;
	my @key=sort {$kk{$a}<=>$kk{$b}} keys %kk;
	#foreach(@key){print "$kk{$_}\t";}print "\n";
	my $key=join "",@key;
	if($key=~/^12/ or $key=~/^34/ or $key=~/^43/ or $key=~/^21/){
		return 0;
	}else{
		return 1;
	}
}

1
