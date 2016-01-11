package Filter_Overlap_BAC_2;
use strict;use warnings;

#huangguodong@genomics.cn : <m8> <scaffold_len> 

#my %len;
#open IN,$ARGV[1] or die;
#while(<IN>){
#  my @aa=split;
#  $len{$aa[0]}=$aa[1];
#}
#close IN;

#filter($ARGV[0],\%len);

sub filter{
	my ($self,$m8,$hash,$del,$scf_len,$bacSort)=@_;
#	my ($m8)=@_;
	open IN,$m8 or die;
	my %tmp;
	my %check_reverse;
	my $num=1;
	while(<IN>){
	    my @aa=split;
		if($aa[11]<100){next;}
		if($aa[0] eq $aa[1]){next;}
		my $p_id1=$aa[0];if($del->{$p_id1}){next;}
		$p_id1=~s/_\d+$//;
		my $p_id2=$aa[1];if($del->{$p_id2}){next;}
		$p_id2=~s/_\d+$//;
		my $check_id="$p_id2#$p_id1";
		if($check_reverse{$check_id}){next;}
		my $type="+";
		if($aa[8]>$aa[9]){($aa[8],$aa[9])=($aa[9],$aa[8]);$type="-";}
		if($aa[12]){$type=$aa[12];}
		$tmp{$aa[0]}{$aa[1]}{$num}=[$aa[6],$aa[7],$aa[8],$aa[9],$type,$aa[2],$aa[3],$aa[4],$aa[5],$aa[10],$aa[11]];
		$num++;
		$check_reverse{"$p_id1#$p_id2"}=1;
	}
	close IN;
	$num=1;
	foreach my $key(keys %tmp){
	    foreach my $key1(keys $tmp{$key}){
		  	#check contain
		  	print STDERR "check\t$key\t$key1\n";
			my @pre=();my $pre=0;
		  	my %m_tmp;
		  	foreach my $key2(sort {$tmp{$key}{$key1}{$a}[0]<=>$tmp{$key}{$key1}{$b}[0]} keys $tmp{$key}{$key1}){
			    my $st1=$tmp{$key}{$key1}{$key2}[0];my $ed1=$tmp{$key}{$key1}{$key2}[1];
			    my $st2=$tmp{$key}{$key1}{$key2}[2];my $ed2=$tmp{$key}{$key1}{$key2}[3];
				my $type=$tmp{$key}{$key1}{$key2}[4];
				if($ed1-$st1<300 or $ed2-$st2<300){next;}
				unless(@pre){@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;next;}
				my %b_tmp;
				$b_tmp{"1"}=$pre[2];$b_tmp{"2"}=$pre[3];
				$b_tmp{"3"}=$st2;$b_tmp{"4"}=$ed2;
				my @k_b=sort {$b_tmp{$a}<=>$b_tmp{$b}} keys %b_tmp;my $k_b=join "",@k_b;
				my %a_tmp;
				$a_tmp{"1"}=$pre[0];$a_tmp{"2"}=$pre[1];
				$a_tmp{"3"}=$st1;$a_tmp{"4"}=$ed1;
				my @k_a=sort {$a_tmp{$a}<=>$a_tmp{$b}} keys %a_tmp;my $k_a=join "",@k_a;
				my $c_a=check($k_a);my $c_b=check($k_b);
				if($c_a==2){
					#contain
				  	if($c_b==0){
				    	my $len1=$st2-$pre[3]+1;my $len2=$pre[2]-$ed2+1;
						my $len;if($len1>0){$len=$len1;}else{$len=$len2;}
						if($len>500){
							if($pre[0] != $st1){
							  	if($k_a=~/^1/){
						        	next;
						      	}
						      	if($k_a=~/^3/){
						        	@pre=($st1,$ed1,$st2,$ed2);$pre=$key2;next;
						      	}
					    	}else{
					    		if($pre[1]<$ed1){
					      			@pre=($st1,$ed1,$st2,$ed2);$pre=$key2;next;
					      		}else{
					      			next;
					      		}
					    	}
						}
				  	}
				  	my $n_st1=$a_tmp{$k_a[0]};my $n_ed1=$a_tmp{$k_a[-1]};
				  	my $n_st2=$b_tmp{$k_b[0]};my $n_ed2=$b_tmp{$k_b[-1]};
				  	$type="$pre[4]$type";
				  	@pre=($n_st1,$n_ed1,$n_st2,$n_ed2,$type);$pre=$key2;
				  	next;
				}
				if($c_a==1){
					#neighbor
				  	if($c_b==0){
					    my $len1=$st2-$pre[3]+1;my $len2=$pre[2]-$ed2+1;
						my $len;
						if($len1>0){$len=$len1;}else{$len=$len2;}
						#print "$len1\t$len2\t$len\n";
						if($len>500){
							@{$m_tmp{$pre}}=@pre;
							#$m_tmp{$key2}=[$st1,$ed1,$st2,$ed2];
							#@pre=();
							@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;
							next;
						}
				  	}
					my $n_st1=$a_tmp{$k_a[0]};my $n_ed1=$a_tmp{$k_a[-1]};
					my $n_st2=$b_tmp{$k_b[0]};my $n_ed2=$b_tmp{$k_b[-1]};
					$type="$pre[4]$type";
					@pre=($n_st1,$n_ed1,$n_st2,$n_ed2,$type);$pre=$key2;
				}
				if($c_a==0){
					#overlap
				  	#if($c_b==2){
				  	#	if($pre[2] != $st2){
				    #		if($k_b=~/^1/){next;}
					#		if($k_b=~/^3/){@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;next;}
					#	}else{
					#		if($pre[3]<$ed2){
					#			@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;next;
					#		}else{
					#			next;
					#		}
					#	}
				  	#}
				  	my $len1=$pre[0]-$ed1+1;my $len2=$st1-$pre[1]+1;
				  	my $len;if($len1>0){$len=$len1;}else{$len=$len2;}
				  	if($len>500){
						@{$m_tmp{$pre}}=@pre;
						##$m_tmp{$key2}=[$st1,$ed1,$st2,$ed2];
						##@pre=();
						@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;
						next;
				  	}
				  	my $n_st1=$a_tmp{$k_a[0]};my $n_ed1=$a_tmp{$k_a[-1]};
				  	if($c_b==1 or $c_b==2){
					    #my $n_st1=$k_a[0];my $n_ed1=$k_a[-1];
						my $n_st2=$b_tmp{$k_b[0]};my $n_ed2=$b_tmp{$k_b[-1]};
						$type="$pre[4]$type";
						@pre=($n_st1,$n_ed1,$n_st2,$n_ed2,$type);$pre=$key2;
						next;
				  	}
				  	if($c_b==0){
					    $len1=$st2-$pre[3]+1;$len2=$pre[2]-$ed2+1;
						if($len1>0){$len=$len1;}else{$len=$len2;}
						if($len>500){
						  	@{$m_tmp{$pre}}=@pre;
						  	#$m_tmp{$key2}=[$st1,$ed1,$st2,$ed2];
						  	#@pre=();
						  	@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;
						 	next;
						}
						my $n_st2=$b_tmp{$k_b[0]};my $n_ed2=$b_tmp{$k_b[-1]};
						#$pre[2]=$n_st2;$pre[3]=$n_ed2;
						$type="$pre[4]$type";
						@pre=($n_st1,$n_ed1,$n_st2,$n_ed2,$type);$pre=$key2;
				  	}
				}#if
		    }#foreach
		  	if(@pre){@{$m_tmp{$pre}}=@pre;}
		  	#foreach my $ass(sort {$a<=>$b} keys %m_tmp){
		    #	my $out=join "\t",@{$m_tmp{$ass}};
			#	#print "$key\t$key1\t$ass\t$out\n";
			#	#print "$key1\t$key\t90\t100\t0\t0\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t0.0\t500\n";
			#	$hash->{$key}{$key1}{$num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
			#	print STDERR "final\t$key\t$key1\t$num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
			#	$num++;
			#	#
		  	#}
		  	#next;
		  	@pre=();
		  	#find the only best record for the next analysis
		  	my @last=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} keys %m_tmp;
		  	my @last2=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} keys %m_tmp;
		  	if(@last>1){
		  		#check end
		  		my $cut1=20;my $cut2=500;
    			my $bac1=$key;my $bac1_id=$bac1;$bac1_id=~s/_\d+$//;
    			my $bac2=$key1;my $bac2_id=$bac2;$bac2_id=~s/_\d+$//;
    			my %block;my %block_rev;
    			block(\@last,\@last2,\%m_tmp,\%block,\%block_rev);
    			my %mm_tmp;my @longest1=(0,0);my @longest2=(0,0);
    			#try to find the block that can merge which not just by neighbor
    			foreach my $rec(@last){
    				my $t_len1=$m_tmp{$rec}[1]-$m_tmp{$rec}[0]+1;
    				my $t_len2=$m_tmp{$rec}[3]-$m_tmp{$rec}[2]+1;
    				if($t_len1>$longest1[0]){@longest1=($t_len1,$rec);}
    				if($t_len2>$longest2[0]){@longest2=($t_len2,$rec);}
    			}
    			if($longest1[1] == $longest2[1]){
    				#best is the same block and try to extend by the block information
    				my $bloaa=$block_rev{$longest1[1]}[0];
    				my @aa=keys $block{$bloaa};
    				my $blobb=$block_rev{$longest2[1]}[1];
    				my @bb=keys $block{$blobb};
    				my @keep=@aa;
    				print STDERR "$longest1[1] == $longest2[1]\n";
    				if(@keep>1){
    					#mutiple block
    					my @keep_tmp=check_block(\@aa,\@bb,$bloaa,$blobb,$longest1[1],$longest2[1],\%block_rev);
    					foreach(@keep_tmp){print STDERR $_,"\n";}
    					my $check_block_tmp=1;
    					my $control=@keep;
    					my @last_tmp=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} @keep_tmp;
    					my @last2_tmp=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} @keep_tmp;
    					while(@keep_tmp<$control){
    						#rebuld block information
    						$control=@keep_tmp;
    						my %block_tmp;my %block_rev_tmp;
    						block(\@last_tmp,\@last2_tmp,\%m_tmp,\%block_tmp,\%block_rev_tmp);
    						my $bloaa_tmp=$block_rev_tmp{$longest1[1]}[0];
    						my $blobb_tmp=$block_rev_tmp{$longest2[1]}[1];
    						my @aa_tmp=keys $block_tmp{$bloaa_tmp};my @bb_tmp=keys $block_tmp{$blobb_tmp};
    						@keep_tmp=check_block(\@aa_tmp,\@bb_tmp,$bloaa,$blobb,$longest1[1],$longest2[1],\%block_rev_tmp);
    					}
    					@keep=@keep_tmp;
    				}
    				my $c_num=1;
			  		foreach my $ass(sort {$a<=>$b} @keep){
			  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
			  			print STDERR "filter2deal\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
			  			$c_num++;
			  		}
    			}else{
    				#have two best block
    				#first check if can link. if not choose the best.
    				print STDERR "$key\t$key1\ttwo best\n";
    				#check the unmap length between two best
    				my ($first,$last)=sort {$a<=>$b} ($longest1[1],$longest2[1]);

    				my $len1=$m_tmp{$first}[1]-$m_tmp{$first}[0]+1;
    				my $len2=$m_tmp{$last}[3]-$m_tmp{$last}[2]+1;
    				if($len1<=2000 or $len2 <= 2000){print STDERR "too short to trust\n";}

    				#check st ed
    				my $cut1=20;my $cut2=500;
    				my $bac1=$key;my $bac1_id=$bac1;$bac1_id=~s/_\d+$//;
    				my $bac2=$key1;my $bac2_id=$bac2;$bac2_id=~s/_\d+$//;
    				my $bac1_st;my $bac1_ed;my $bac2_st;my $bac2_ed;
    				foreach my $rect(@last){
    					if($m_tmp{$rect}[0]<=$cut1){$bac1_st=$rect;}
    					if($scf_len->{$bac1}-$m_tmp{$rect}[0]<=$cut1){$bac1_ed=$rect;}
    					if($m_tmp{$rect}[2]<=$cut1){$bac2_st=$rect;}
    					if($scf_len->{$bac2}-$m_tmp{$rect}[3]<=$cut1){$bac2_ed=$rect;}
    				}

    				my $bloaa=$block_rev{$longest1[1]}[0];
    				my @aa=keys $block{$bloaa};
    				my $blobb=$block_rev{$longest2[1]}[1];
    				my @bb=keys $block{$blobb};

    				if(($bac1_st or $bac1_ed) && ($bac2_st or $bac2_ed)){
    					#get st or ed. may be get overlap relation
    					print STDERR "$key\t$key1\t$scf_len->{$key}\t$scf_len->{$key1}\tmay get ovelap relation should check hand\n";
    					#foreach my $ass(sort {$a<=>$b} keys %m_tmp){
    					#	print STDERR "filter2dealtwoall\t$key\t$key1\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
    					#}
    					if($bac1_st){
    						if($bac1_ed){
    							#st and end. should be contain.
    							if($block_rev{$bac1_st}[0] == $block_rev{$bac1_ed}[0]){
    								#can connect together
    								my %tp;
    								foreach my $tp(keys %block_rev){
    									$tp{$block_rev{$tp}[1]}=1;
    								}
    								if(keys %tp == 1){
    									my $c_num=1;
								  		foreach my $ass(sort {$a<=>$b} keys %m_tmp){
								  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
								  			print STDERR "filter2dealtwocontain\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
								  			$c_num++;
								  		}
    								}
    							}
    						}else{
    							#only st and 
    							if($bac2_st && $bac2_ed){
    								#bac2 in contain
    								if($block_rev{$bac2_st}[1] == $block_rev{$bac2_ed}[1]){
    									my %tp;
	    								foreach my $tp(keys %block_rev){
	    									$tp{$block_rev{$tp}[0]}=1;
	    								}
	    								if(keys %tp == 1){
	    									my $c_num=1;
									  		foreach my $ass(sort {$a<=>$b} keys %m_tmp){
									  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
									  			print STDERR "filter2dealtwocontain\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
									  			$c_num++;
									  		}
	    								}
    								}
    							}elsif($bac2_st or $bac2_ed){
    								#only bac2 st: bac1 st --- bac2 st or ed
    								if($bac2_ed){$bac2_st=$bac2_ed;}
    								if($bac1_st==$bac2_st){
    									#same block record.try extend block;
    									my $bloaa=$block_rev{$bac1_st}[0];
					    				my @aa=keys $block{$bloaa};
					    				my $blobb=$block_rev{$bac2_st}[1];
					    				my @bb=keys $block{$blobb};
					    				my @keep=@aa;
					    				print STDERR "$bac1_st == $bac2_st\n";
					    				if(@keep>1){
					    					#mutiple block
					    					my @keep_tmp=check_block(\@aa,\@bb,$bloaa,$blobb,$bac1_st,$bac2_st,\%block_rev);
					    					foreach(@keep_tmp){print STDERR $_,"\n";}
					    					my $check_block_tmp=1;
					    					my $control=@keep;
					    					my @last_tmp=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} @keep_tmp;
					    					my @last2_tmp=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} @keep_tmp;
					    					while(@keep_tmp<$control){
					    						#rebuld block information
					    						$control=@keep_tmp;
					    						my %block_tmp;my %block_rev_tmp;
					    						block(\@last_tmp,\@last2_tmp,\%m_tmp,\%block_tmp,\%block_rev_tmp);
					    						my $bloaa_tmp=$block_rev_tmp{$bac1_st}[0];
					    						my $blobb_tmp=$block_rev_tmp{$bac2_st}[1];
					    						my @aa_tmp=keys $block_tmp{$bloaa_tmp};my @bb_tmp=keys $block_tmp{$blobb_tmp};
					    						@keep_tmp=check_block(\@aa_tmp,\@bb_tmp,$bloaa,$blobb,$bac1_st,$bac2_st,\%block_rev_tmp);
					    					}
					    					@keep=@keep_tmp;
					    				}
					    				my $c_num=1;
								  		foreach my $ass(sort {$a<=>$b} @keep){
								  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
								  			print STDERR "filter2deal\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
								  			$c_num++;
								  		}
    								}else{
    									#check if can link
    									if($block_rev{$bac1_st}[1] == $block_rev{$bac2_st}[1] && $block_rev{$bac1_st}[0] == $block_rev{$bac2_st}[0]){
    										#can link; output the block
    										my $longest1=$bac1_st;my $longest2=$bac1_st;
    										my $bloaa=$block_rev{$longest1}[0];
						    				my @aa=keys $block{$bloaa};
						    				my $blobb=$block_rev{$longest2}[1];
						    				my @bb=keys $block{$blobb};
						    				my @keep=@aa;
						    				print STDERR "$longest1 == $longest2\n";
						    				if(@keep>1){
						    					#mutiple block
						    					my @keep_tmp=check_block(\@aa,\@bb,$bloaa,$blobb,$longest1,$longest2,\%block_rev);
						    					foreach(@keep_tmp){print STDERR $_,"\n";}
						    					my $check_block_tmp=1;
						    					my $control=@keep;
						    					my @last_tmp=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} @keep_tmp;
						    					my @last2_tmp=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} @keep_tmp;
						    					while(@keep_tmp<$control){
						    						#rebuld block information
						    						$control=@keep_tmp;
						    						my %block_tmp;my %block_rev_tmp;
						    						block(\@last_tmp,\@last2_tmp,\%m_tmp,\%block_tmp,\%block_rev_tmp);
						    						my $bloaa_tmp=$block_rev_tmp{$longest1}[0];
						    						my $blobb_tmp=$block_rev_tmp{$longest2}[1];
						    						my @aa_tmp=keys $block_tmp{$bloaa_tmp};my @bb_tmp=keys $block_tmp{$blobb_tmp};
						    						@keep_tmp=check_block(\@aa_tmp,\@bb_tmp,$bloaa,$blobb,$longest1,$longest2,\%block_rev_tmp);
						    					}
						    					@keep=@keep_tmp;
						    				}
						    				my $flag=0;
						    				foreach my $rect_aa(@keep){
						    					if($rect_aa==$bac2_st){$flag=1;}
						    				}
						    				if($flag){
						    					my $c_num=1;
										  		foreach my $ass(sort {$a<=>$b} @keep){
										  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
										  			print STDERR "filter2deal\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
										  			$c_num++;
										  		}
						    				}
    									}
    								}
    							}
    						}
    					}else{
    						if($bac1_ed){
    							#only bac1-ed
    							if($bac2_st && $bac2_ed){
    								#bac2 in contain
    								if($block_rev{$bac2_st}[1] == $block_rev{$bac2_ed}[1]){
    									my %tp;
	    								foreach my $tp(keys %block_rev){
	    									$tp{$block_rev{$tp}[0]}=1;
	    								}
	    								if(keys %tp == 1){
	    									my $c_num=1;
									  		foreach my $ass(sort {$a<=>$b} keys %m_tmp){
									  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
									  			print STDERR "filter2dealtwocontain\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
									  			$c_num++;
									  		}
	    								}
    								}
    							}elsif($bac2_st or $bac2_ed){
    								#only bac2 st: bac1 st --- bac2 st or ed
    								if($bac2_ed){$bac2_st=$bac2_ed;}
    								if($bac1_ed==$bac2_st){
    									#same block record.try extend block;
    									my $bloaa=$block_rev{$bac1_ed}[0];
					    				my @aa=keys $block{$bloaa};
					    				my $blobb=$block_rev{$bac2_st}[1];
					    				my @bb=keys $block{$blobb};
					    				my @keep=@aa;
					    				print STDERR "$bac1_ed == $bac2_st\n";
					    				if(@keep>1){
					    					#mutiple block
					    					my @keep_tmp=check_block(\@aa,\@bb,$bloaa,$blobb,$bac1_ed,$bac2_st,\%block_rev);
					    					foreach(@keep_tmp){print STDERR $_,"\n";}
					    					my $check_block_tmp=1;
					    					my $control=@keep;
					    					my @last_tmp=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} @keep_tmp;
					    					my @last2_tmp=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} @keep_tmp;
					    					while(@keep_tmp<$control){
					    						#rebuld block information
					    						$control=@keep_tmp;
					    						my %block_tmp;my %block_rev_tmp;
					    						block(\@last_tmp,\@last2_tmp,\%m_tmp,\%block_tmp,\%block_rev_tmp);
					    						my $bloaa_tmp=$block_rev_tmp{$bac1_ed}[0];
					    						my $blobb_tmp=$block_rev_tmp{$bac2_st}[1];
					    						my @aa_tmp=keys $block_tmp{$bloaa_tmp};my @bb_tmp=keys $block_tmp{$blobb_tmp};
					    						@keep_tmp=check_block(\@aa_tmp,\@bb_tmp,$bloaa,$blobb,$bac1_ed,$bac2_st,\%block_rev_tmp);
					    					}
					    					@keep=@keep_tmp;
					    				}
					    				my $c_num=1;
								  		foreach my $ass(sort {$a<=>$b} @keep){
								  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
								  			print STDERR "filter2deal\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
								  			$c_num++;
								  		}
    								}else{
    									#check if can link
    									if($block_rev{$bac1_ed}[1] == $block_rev{$bac2_st}[1] && $block_rev{$bac1_ed}[0] == $block_rev{$bac2_st}[0]){
    										#can link; output the block
    										my $longest1=$bac1_ed;my $longest2=$bac1_ed;
    										my $bloaa=$block_rev{$longest1}[0];
						    				my @aa=keys $block{$bloaa};
						    				my $blobb=$block_rev{$longest2}[1];
						    				my @bb=keys $block{$blobb};
						    				my @keep=@aa;
						    				print STDERR "$longest1 == $longest2\n";
						    				if(@keep>1){
						    					#mutiple block
						    					my @keep_tmp=check_block(\@aa,\@bb,$bloaa,$blobb,$longest1,$longest2,\%block_rev);
						    					foreach(@keep_tmp){print STDERR $_,"\n";}
						    					my $check_block_tmp=1;
						    					my $control=@keep;
						    					my @last_tmp=sort {$m_tmp{$a}[0]<=>$m_tmp{$b}[0]} @keep_tmp;
						    					my @last2_tmp=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} @keep_tmp;
						    					while(@keep_tmp<$control){
						    						#rebuld block information
						    						$control=@keep_tmp;
						    						my %block_tmp;my %block_rev_tmp;
						    						block(\@last_tmp,\@last2_tmp,\%m_tmp,\%block_tmp,\%block_rev_tmp);
						    						my $bloaa_tmp=$block_rev_tmp{$longest1}[0];
						    						my $blobb_tmp=$block_rev_tmp{$longest2}[1];
						    						my @aa_tmp=keys $block_tmp{$bloaa_tmp};my @bb_tmp=keys $block_tmp{$blobb_tmp};
						    						@keep_tmp=check_block(\@aa_tmp,\@bb_tmp,$bloaa,$blobb,$longest1,$longest2,\%block_rev_tmp);
						    					}
						    					@keep=@keep_tmp;
						    				}
						    				my $flag=0;
						    				foreach my $rect_aa(@keep){
						    					if($rect_aa==$bac2_st){$flag=1;}
						    				}
						    				if($flag){
						    					my $c_num=1;
										  		foreach my $ass(sort {$a<=>$b} @keep){
										  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
										  			print STDERR "filter2deal\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
										  			$c_num++;
										  		}
						    				}
    									}
    								}
    							}
    						}
    					}
    				}else{
    					print STDERR "$key\t$key1\tno overlap relation detected\n";
    				}

    				#
    				my $ed_aa_1=$m_tmp{$first}[1];my $st_aa_2=$m_tmp{$last}[0];
    				my @rect_aa;my $flag=1;
    				foreach my $rect(@last){
    					if($m_tmp{$rect}[1] eq $ed_aa_1){$flag=2;}
    					if($m_tmp{$rect}[0] eq $st_aa_2){$flag=3;}
    					if($flag==1){next;}
    					if($flag==3){last;}
    					#
    				}

    			}
		  	}else{
		  		my $c_num=1;
		  		foreach my $ass(sort {$a<=>$b} keys %m_tmp){
		  			$hash->{$key}{$key1}{$c_num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
		  			print STDERR "filter2lone\t$key\t$key1\t$c_num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
		  			$c_num++;
		  		}
		  	}	
		}#foreach
	}#foreach
}


sub check{
  my ($string)=@_;
#  print "$string\n";
  if($string=~/^1..2$/ or $string=~/^3..4$/){
    #contain
	return 2;
  }
  if($string=~/^12..$/ or $string=~/^34..$/){
  	#neibor
	return 0;
  }
  if($string=~/^1.2.$/ or $string=~/^.1.2$/){
  	#overlap
	return 1;
  }
}

sub block{
	my ($last,$last2,$m_tmp,$blockc,$block_rev)=@_;
	my ($st1,$ed1,$st2,$ed2);my $block=1;my $pre_rec;
    	foreach my $rec(@{$last}){
    		unless($st1){
 				$st1=$m_tmp->{$rec}[0];$ed1=$m_tmp->{$rec}[1];
    			#$st2=$m_tmp{$rec}[2];$ed2=$m_tmp{$rec}[3];
    			$pre_rec=$rec;
    			$blockc->{$block}{$rec}=1;push @{$block_rev->{$rec}},$block;
    			next;
    		}
    		if($m_tmp->{$rec}[0]-$ed1>1000){
    			$block++;
    		}
    		$blockc->{$block}{$rec}=1;push @{$block_rev->{$rec}},$block;
    		$st1=$m_tmp->{$rec}[0];$ed1=$m_tmp->{$rec}[1];
    	}
    	$block++;
    	foreach my $rec(@{$last2}){
    		unless($st2){
    			$st2=$m_tmp->{$rec}[2];$ed2=$m_tmp->{$rec}[3];
    			#$st2=$m_tmp{$rec}[2];$ed2=$m_tmp{$rec}[3];
    			$pre_rec=$rec;$blockc->{$block}{$rec}=1;push @{$block_rev->{$rec}},$block;
    			next;
    		}
    		if($m_tmp->{$rec}[2]-$ed2>1000){
    			$block++;
    		}
    		$blockc->{$block}{$rec}=1;push @{$block_rev->{$rec}},$block;
    		$st2=$m_tmp->{$rec}[2];$ed2=$m_tmp->{$rec}[3];
    	}
}

sub check_block{
	my ($aa,$bb,$bloaa,$blobb,$long1,$long2,$block_rev)=@_;
	my @aa=@{$aa};my @bb=@{$bb};
	my @aa_pre=grep {$_ < $long1} @aa;
    my @aa_aft=grep {$_ > $long1} @aa;
    my @keep_aa;
    if(@aa_pre>0){
	    foreach my $rect(sort {$b<=>$a} @aa_pre){
	    	my $blo_bb_tmp=$block_rev->{$rect}[1];
	    	if($blo_bb_tmp==$blobb){push @keep_aa,$rect;}
	    	else{last;}
	    }
    }
    if(@aa_aft>0){
    	foreach my $rect(sort {$a<=>$b} @aa_aft){
    		my $blo_bb_tmp=$block_rev->{$rect}[1];
    		if($blo_bb_tmp==$blobb){push @keep_aa,$rect;}
    		else{last;}
    	}
    }
    push @keep_aa,$long1;
    #deal bb
    my @bb_pre=grep {$_ < $long2} @bb;
    my @bb_aft=grep {$_ > $long2} @bb;
    my @keep_bb;
    if(@bb_pre>0){
    	foreach my $rect(sort {$b<=>$a} @bb_pre){
			my $blo_aa_tmp=$block_rev->{$rect}[0];
			if($blo_aa_tmp==$bloaa){push @keep_bb,$rect;}
    		else{last;}
    	}
    }
	if(@bb_aft>0){
		foreach my $rect(sort {$a<=>$b} @bb_aft){
			my $blo_aa_tmp=$block_rev->{$rect}[0];
			if($blo_aa_tmp==$bloaa){push @keep_bb,$rect;}
			else{last;}
		}
	}
	push @keep_bb,$long2;
	#get keep_aa and keep_bb overlap
	my %tmp;
	foreach(@keep_aa){$tmp{$_}++;}foreach(@keep_bb){$tmp{$_}++;}
	my @keep;
	foreach(keys %tmp){if($tmp{$_}==2){push @keep,$_;}}
	return (@keep);
}



1;
