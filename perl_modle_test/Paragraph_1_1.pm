package Paragraph_1_1;

# huangguodong@genomics.cn

use strict;
use warnings;

sub para{
    my ($self,$hash,$ex_ovl_len,$out,$sort,$bac1_id,$bac2_id)=@_;
    #filter overlap which in contain
    my %use_contain;
    my $flag_con=0;
	my %ovl_blo;my %block;my $block=1;
	my $flag=0;
	my $bac1=0;my $bac2=0;
	my %scf_len;
    if($hash->{"contain"}){
		$flag++;
		$flag_con=1;
		foreach my $num(keys $hash->{"contain"}){
	    	my $tag=$hash->{"contain"}{$num}[0];
			my ($st1,$ed1,$st2,$ed2)=($hash->{"contain"}{$num}[2],$hash->{"contain"}{$num}[3],$hash->{"contain"}{$num}[4],$hash->{"contain"}{$num}[5]);
			my ($b_len1,$b_len2)=($hash->{"contain"}{$num}[6],$hash->{"contain"}{$num}[7]);
	    	my @tag=split /\|/,$tag;
			unless($bac1){if($tag[0]=~/bac_\d+_\d+/){$bac1=$&;}}
			unless($bac2){if($tag[1]=~/bac_\d+_\d+/){$bac2=$&;}}
			my $cu=0;
			my $check_C=0;
			foreach my $id(@tag){
				if($id=~/^C/){
					$check_C++;
				}
			}
			if($check_C==2){
				my $id1=(split /\s+/,$tag[0])[1];
				my $id2=(split /\s+/,$tag[1])[1];
				if($b_len1>$b_len2){
					my $type="pre";
					$block{$id1}{$id2}=[$type,$st1,$ed1,$st2,$ed2,$b_len1,$b_len2];
					$scf_len{$id1}=$b_len1;$scf_len{$id2}=$b_len2;
				}else{
					my $type="aft";
					$block{$id2}{$id1}=[$type,$st2,$ed2,$st1,$ed1,$b_len2,$b_len1];
					$scf_len{$id1}=$b_len1;$scf_len{$id2}=$b_len2;
				}
			}else{
				foreach my $id(@tag){
					if($id=~/^C/){
						my @id=split /\s+/,$id;
						my $id_a=$id[1];
						my $tmp_cu;my $type;
						if($cu){$tmp_cu=$cu-1;$type="pre";}else{$tmp_cu=$cu+1;$type="aft";($st1,$ed1,$st2,$ed2,$b_len1,$b_len2)=($st2,$ed2,$st1,$ed1,$b_len2,$b_len1);}
						$block{$tag[$tmp_cu]}{$id_a}=[$type,$st1,$ed1,$st2,$ed2,$b_len1,$b_len2];
						$scf_len{$tag[$tmp_cu]}=$b_len1;$scf_len{$id_a}=$b_len2;
					}
					$cu++;
				}
			}
	    	foreach my $id(@tag){
				if($id=~/^C/){
		    		my @id=split /\s+/,$id;
		    		$use_contain{$id[1]}++;;
		    		next;
				}
				if($use_contain{$id}){next;}else{$use_contain{$id}=0;}
	    	}
		}
    }
    if($hash->{"overlap"}){
		foreach my $num(keys $hash->{"overlap"}){
	    	#filter overlap which in contain
	    	my $tag=$hash->{"overlap"}{$num}[0];
			my ($st1,$ed1,$st2,$ed2)=($hash->{"overlap"}{$num}[2],$hash->{"overlap"}{$num}[3],$hash->{"overlap"}{$num}[4],$hash->{"overlap"}{$num}[5]);
			my ($b_len1,$b_len2)=($hash->{"overlap"}{$num}[6],$hash->{"overlap"}{$num}[7]);
			my $len1=$ed1-$st1;my $len2=$ed2-$st2;
			my @tag=split /\|/,$tag;
			unless($bac1){if($tag[0]=~/bac_\d+_\d+/){$bac1=$&;}}
			unless($bac2){if($tag[1]=~/bac_\d+_\d+/){$bac2=$&;}}
			my $dv_len;
			if($len1>$len2){$dv_len=$len2;}else{$dv_len=$len1;}
			if($use_contain{$tag[0]} || $use_contain{$tag[1]}){
				#print STDERR "$tag[0]\t$tag[1] overlap is not real\n";
				print "$tag[0]\t$tag[1] overlap is not real\n";
			}elsif($len1-$len2>500 && ($len1-$len2)/$dv_len*100>20){
				#print STDERR "$tag[0]\t$tag[1] overlap length distribute too large\n";
				print "$tag[0]\t$tag[1] overlap length distribute too large\n";
			}elsif($len1<=200 || $len2<=200){
				#print STDERR "$tag[0]\t$tag[1] overlap length is too small to trust\n";
				print "$tag[0]\t$tag[1] overlap length is too small to trust\n";
			}else{
				#$ovl_blo{"overlap"}=$hash->{"overlap"};
				#make overlap block ready
				print "$tag[0]\t$tag[1]\tovl_block\n";
				print "$st1,$ed1,$st2,$ed2,$b_len1,$b_len2,$tag[2],$tag[3]\n";
				$ovl_blo{$tag[0]}{$tag[1]}=["$tag[2]|$tag[3]",$st1,$ed1,$st2,$ed2,$b_len1,$b_len2,$tag[2],$tag[3]];
				#$ovl_blo{$tag[1]}{$tag[0]}=["$tag[3]|$tag[2]",$st2,$ed2,$st1,$ed1,$b_len2,$b_len1,$tag[3],$tag[2]];
				$scf_len{$tag[0]}=$b_len1;$scf_len{$tag[1]}=$b_len2;
				$flag++;
			}
		}
    }
	#build and sort block
	if($flag){
		#check contain block
		foreach my $key1(keys %block){
			my $big_pre_ed=0;
			foreach my $key2(sort {$block{$key1}{$a}[1]<=>$block{$key1}{$b}[1]} keys $block{$key1}){
				unless($big_pre_ed){$big_pre_ed=$block{$key1}{$key2}[2];next;}
				if($block{$key1}{$key2}[2] > $big_pre_ed){
					if($block{$key1}{$key2}[1]>$big_pre_ed && $block{$key1}{$key2}[1]-$big_pre_ed>=5000){
						print "find block biggest than 5Kbp\t$key1\t$key2\t$block{$key1}{$key2}[1]\t$block{$key1}{$key2}[2]\n";
					}
					$big_pre_ed=$block{$key1}{$key2}[2];
				}
			}
		}
		#find block
		my @bac_key=keys %block;
		my @bac_key1=grep /$bac1/,@bac_key;
		my @bac_key2=grep /$bac2/,@bac_key;
		#find block relation depend on overlap
		my @ovl_keys=keys %ovl_blo;
		my %rel_ovl;my %rel_ovl_reverse;my %rel_add;
		if(@ovl_keys != 0){
			#get overlap and check relation with block
			foreach my $key(@ovl_keys){
				my @key1=keys $ovl_blo{$key};
				my $key1=(keys $ovl_blo{$key})[0];
				my @check_key;
				if(@key1>1){
					#should check the repeat scaffold is link or not
					print "$key is special scaffold hand check please\tcheck if link next BAC\n";
					if(@key1>2){
						#multi-link. get the right link.
						print "fuck you . why you are so complex!\n";
						my @tag_sel;
						my %chose;
						foreach my $key_tmp(@key1){
							push @tag_sel,$ovl_blo{$key}{$key_tmp}[7];
							$chose{$key_tmp}=$ovl_blo{$key}{$key_tmp}[6];
						}
						my @stt=grep /st/,@tag_sel;my @edd=grep /ed/,@tag_sel;
						if (@stt>0 and @edd>0){
							my $st_len;my $ed_len;my $st_sel;my $ed_sel;
							foreach my $key_tmp(@key1){
								if ($ovl_blo{$key}{$key_tmp}[7] eq 'st'){
									if(!$st_len or $st_len < $ovl_blo{$key}{$key_tmp}[6]){
										$st_len=$ovl_blo{$key}{$key_tmp}[6];
										$st_sel=$key_tmp;
									}
								}
								if ($ovl_blo{$key}{$key_tmp}[7] eq 'ed'){
									if(!$ed_len or $ed_len < $ovl_blo{$key}{$key_tmp}[6]){
										$ed_len=$ovl_blo{$key}{$key_tmp}[6];
										$ed_sel=$key_tmp;
									}
								}
							}
							print "$key\tcomplex select $st_sel $ed_sel\n";
							#chose for the st
							push @check_key,$st_sel;
							#chose for the ed
							push @check_key,$ed_sel;
						}else{
							$key1=(sort {$chose{$b}<=>$chose{$a}} keys %chose)[0];
							push @check_key,$key1;
						}
					}elsif(@key1==2){
						#check if it it st ed
						print "two end  check happy\n";
						my $tag1=$ovl_blo{$key}{$key1[0]}[7];
						my $tag2=$ovl_blo{$key}{$key1[1]}[7];
						if($tag1 eq $tag2){
							#not correct link. choose the long scaffold or the long overlap
							if($ovl_blo{$key}{$key1[1]}[6] > $ovl_blo{$key}{$key1[0]}[6]){
								print "$key1\tis redundunt\n";
								$key1=(keys $ovl_blo{$key})[1];
								push @check_key,$key1;
							}else{
								print "keep $key1\n";
							}
						}else{
							push @check_key,@key1;
						}
					}else{
						print "error\n";
					}
				}else{
					push @check_key,@key1;
				}
				#my $key1=(keys $ovl_blo{$key})[0];
				if (@check_key>1){
					#foreach my $check_key(@check_key){
						print "$key\ttwo end need check \n";
						my $check_key=$check_key[0];my $add_key=$check_key[1];
						if($block{$key}){
							if($block{$check_key}){
								#bac1 and bac2 have relation
								#$rel_ovl{$check_key}=$key;
								$rel_ovl{$key}=$check_key;
								$rel_add{$key}=$add_key;
								my @bac_key1_tmp=grep {!($_ eq $key)} @bac_key1;
								@bac_key1=();
								if(@bac_key1_tmp != 0){@bac_key1=@bac_key1_tmp;}
								my @bac_key2_tmp=grep {!($_ eq $check_key)} @bac_key2;
								@bac_key2_tmp=grep {!($_ eq $add_key)} @bac_key2;
								@bac_key2=();
								if(@bac_key2_tmp != 0){@bac_key2=@bac_key2_tmp;}
							}else{
								#$rel_ovl{$check_key}=$key;
								$rel_ovl{$key}=$check_key;
								$rel_add{$key}=$add_key;
							}
						}elsif($block{$check_key}){
							##$rel_ovl{$check_key}=$key;
							#$rel_ovl{$check_key}=$key;
							$rel_ovl{$key}=$check_key;
							$rel_add{$key}=$add_key;
						}else{
							#$rel_ovl{$check_key}=$key;
							$rel_ovl{$key}=$check_key;
							$rel_add{$key}=$add_key;
						}
						$rel_ovl_reverse{$check_key}=$key;
						$rel_ovl_reverse{$add_key}=$key;
					#}
				}else{
					if($block{$key}){
						if($block{$key1}){
							#bac1 and bac2 have relation
							$rel_ovl{$key}=$key1;
							my @bac_key1_tmp=grep {!($_ eq $key)} @bac_key1;
							@bac_key1=();
							if(@bac_key1_tmp != 0){@bac_key1=@bac_key1_tmp;}
							my @bac_key2_tmp=grep {!($_ eq $key1)} @bac_key2;
							@bac_key2=();
							if(@bac_key2_tmp != 0){@bac_key2=@bac_key2_tmp;}
						}else{
							$rel_ovl{$key}=$key1;
						}
					}elsif($block{$key1}){
						#$rel_ovl{$key1}=$key;
						$rel_ovl{$key}=$key1;
					}else{
						$rel_ovl{$key}=$key1;
					}
					$rel_ovl_reverse{$key1}=$key;
				}
			}
		}
		#combine block and check the overlap length to make sure the overlap is correct
		my %st_ed;my $cal_overlap_len=0;
		#my %scf_len;
		my @all_contain_block=(@bac_key1,@bac_key2);
		if(@all_contain_block){
		foreach my $key1(@all_contain_block){
			print "checking non overlap block $key1\n";
			my $pre_st=0;my $pre_ed;
			my @block_cal;
			foreach my $key2(sort {$block{$key1}{$a}[1]<=>$block{$key1}{$b}[1]} keys $block{$key1}){
				unless($pre_st){
					$pre_st=$block{$key1}{$key2}[1];$pre_ed=$block{$key1}{$key2}[2];
					#$scf_len{$key1}=$block{$key1}{$key2}[-2];$scf_len{$key2}=$block{$key1}{$key2}[-1];
					next;
				}
				#$scf_len{$key1}=$block{$key1}{$key2}[-2];$scf_len{$key2}=$block{$key1}{$key2}[-1];
				if($block{$key1}{$key2}[1] > $pre_ed){push @block_cal,$pre_st,$pre_ed;$pre_st=$block{$key1}{$key2}[1];$pre_ed=$block{$key1}{$key2}[2];}
				elsif($block{$key1}{$key2}[1] <= $pre_ed){
					if($block{$key1}{$key2}[2] <= $pre_ed){print "$key1\t$key2\tblock contain by other self\n";next;}
					else{
						$pre_ed=$block{$key1}{$key2}[2];
					}
				}
			}
			if($pre_st){push @block_cal,$pre_st,$pre_ed;}
			#if($rel_ovl{$key1}){
			#	#some problem
			#	#this part should only deal with the overlap scaffold was contained by other non-overlap scaffold
			#	#don't deal with the scaffold contain other contained scaffold
			#	#the overlap relation should be first priority
			#	#this part is wrong 20150817
			#	my $tmp_st;my $tmp_ed;
			#	print "$key1\t$rel_ovl{$key1}\tcontain_overlap_rel be \n";
			#	if($ovl_blo{$key1}{$rel_ovl{$key1}}){$tmp_st=$ovl_blo{$key1}{$rel_ovl{$key1}}[1];$tmp_ed=$ovl_blo{$key1}{$rel_ovl{$key1}}[2];}
			#	elsif($ovl_blo{$rel_ovl{$key1}}{$key1}){$tmp_st=$ovl_blo{$rel_ovl{$key1}}{$key1}[1];$tmp_ed=$ovl_blo{$rel_ovl{$key1}}{$key1}[2];}
			#	else{print "overlap test fail\n";}
			#	for(my $i=0;$i<@block_cal;$i+=2){
			#		my $b_st=$block_cal[$i];my $b_ed=$block_cal[$i+1];
			#		my $check=check_ovl_stat($tmp_st,$tmp_ed,$b_st,$b_ed);
			#		if($check){
			#			if($tmp_st<=$b_st){
			#				$tmp_ed=$b_st;
			#			}elsif($tmp_ed>=$b_ed){
			#				$tmp_st=$b_ed;
			#			}
			#		}
			#		if($tmp_st>=$tmp_ed){last;}
			#	}
			#	if($tmp_st<$tmp_ed){push @block_cal,$tmp_st,$tmp_ed;}
			#}
			if($rel_ovl{$key1}){
				#next;
				;
			}
			my @block_cal_sort=sort {$a<=>$b} @block_cal;
			my $scaf_st=$block_cal_sort[0];my $scaf_ed=$block_cal_sort[-1];
			my $to_ed_len=$scf_len{$key1}-$scaf_ed;
			#check which end belong to the overlap(use the shorter end)
			if($to_ed_len > $scaf_st){
				#start in overlap
				$st_ed{$key1}=["st",$scaf_st,$to_ed_len,$scaf_st,$scaf_ed];print "$key1\tst_ed\tst\n";
				$cal_overlap_len+=$scaf_ed;
			}elsif($to_ed_len < $scaf_st){
				#end in overlap
				$st_ed{$key1}=["ed",$to_ed_len,$scaf_st,$scaf_st,$scaf_ed];print "$key1\tst_ed\ted\n";
				$cal_overlap_len+=$scf_len{$key1}-$scaf_st;
			}else{
				#both is possible
				$st_ed{$key1}=["bo",$scaf_st,$to_ed_len,$scaf_st,$scaf_ed];print "$key1\tst_ed\tbo\n";
				$cal_overlap_len+=$scaf_ed;
			}
		}
		}
		#add the overlap relation to the overlap %ovl_blo
		my %st_ed_ovl;
		my @key_rel_ovl=keys %rel_ovl;
		if(@key_rel_ovl != 0){
			foreach my $key(keys %rel_ovl){
				my $key2=$rel_ovl{$key};
				#if($block{$key}){next;}
				if($rel_add{$key2}){
					print "reverse result; discard $key belong to $key2\n";
				}
				my @tmp_block_check=grep {$_ eq $key} @all_contain_block;
				if(@tmp_block_check){
					print "have contain ele $key\n";
					#may need to discard the record in contain condiction
				#	next;
					#20150817 checking if the scaffold was contained or contain
					print "$key\tcontain block discard\n";
					#$st_ed{$key}=();
					undef $st_ed{$key};delete  $st_ed{$key};
					my @tmp_bac1=grep {!($_ eq $key)} @bac_key1;
					@bac_key1=@tmp_bac1;
				}
				@tmp_block_check=grep {$_ eq $key2} @all_contain_block;
				if(@tmp_block_check){
					print "$key2\tcontain block discard\n";
					#$st_ed{$key2}=();
					undef $st_ed{$key2};delete  $st_ed{$key2};
					my @tmp_bac2=grep {!($_ eq $key2)} @bac_key2;
					@bac_key2=@tmp_bac2;
				}
				print "checking overlap $key\n";
				if($ovl_blo{$key}{$key2}[5]<$ovl_blo{$key}{$key2}[6]){
					#key in overlap
					$cal_overlap_len+=$ovl_blo{$key}{$key2}[5];
				}elsif($ovl_blo{$key}{$key2}[5]>=$ovl_blo{$key}{$key2}[6]){
					#key2 in overlap
					$cal_overlap_len+=$ovl_blo{$key}{$key2}[6];
				}
				#$scf_len{$key}=$ovl_blo{$key}{$key2}[-2];$scf_len{$key2}=$ovl_blo{$key}{$key2}[-1];
				my $key_type=$ovl_blo{$key}{$key2}[7];
				if($key_type eq "st"){
					my $key_type_len=$ovl_blo{$key}{$key2}[5]-$ovl_blo{$key}{$key2}[2];
					$st_ed_ovl{"1"}{$key}=$key_type_len;
				}elsif($key_type eq "ed"){
					my $key_type_len=$ovl_blo{$key}{$key2}[1];$st_ed_ovl{"1"}{$key}=$key_type_len;#may cause a bug
				}
				my $key2_type=$ovl_blo{$key}{$key2}[8];
				if($key2_type eq "st"){
					my $key2_type_len=$ovl_blo{$key}{$key2}[6]-$ovl_blo{$key}{$key2}[4];$st_ed_ovl{"2"}{$key2}=[$key2_type_len,$key];
				}elsif($key2_type eq "ed"){
					my $key2_type_len=$ovl_blo{$key}{$key2}[1];$st_ed_ovl{"2"}{$key2}=[$key2_type_len,$key];
				}
			}
		}
		print "$bac1\t$bac2\toverlap length\t$cal_overlap_len\t$ex_ovl_len\n";
		if($ex_ovl_len < $cal_overlap_len){print "expect overlap length is shorter than actual\n";}
		elsif($ex_ovl_len == $cal_overlap_len){print "expect overlap length is equal actual\n";}
		else{print "expect overlap length is longer than actual\n";}
		
		#find the two end and place the scaffold to make the overlap
		#%st_ed @bac_key1 @bac_key2 
		#make sure first end
		my %tmp_me;
		foreach my $scf_1(@bac_key1){
			if($st_ed{$scf_1}){$tmp_me{$scf_1}=$st_ed{$scf_1};}
		}
		my $first_end_scf=(sort {$tmp_me{$a}[2]<=>$tmp_me{$b}[2]} keys %tmp_me)[-1];
		my $first_end_scf_len=$tmp_me{$first_end_scf}[2];
		#check overlap end
		if($st_ed_ovl{"1"}){
			foreach my $key_st_ed_ovl_1(keys $st_ed_ovl{"1"}){
				if($st_ed_ovl{"1"}{$key_st_ed_ovl_1}>$first_end_scf_len){
					$first_end_scf=$key_st_ed_ovl_1;$first_end_scf_len=$st_ed_ovl{"1"}{$key_st_ed_ovl_1};
				}
			}
		}
		#make sure the second end
		my %tmp_ms;
		foreach my $scf_2(@bac_key2){
			if($st_ed{$scf_2}){$tmp_ms{$scf_2}=$st_ed{$scf_2};}
		}
		my $secon_end_scf=(sort {$tmp_ms{$a}[2]<=>$tmp_ms{$b}[2]} keys %tmp_ms)[-1];
		my $secon_end_scf_len;
		if($secon_end_scf){$secon_end_scf_len=$tmp_ms{$secon_end_scf}[2];}
		#check overlap end
		if($st_ed_ovl{"2"}){
		foreach my $key_st_ed_ovl_2(keys $st_ed_ovl{"2"}){
			if($st_ed_ovl{"2"}{$key_st_ed_ovl_2}[0]>$secon_end_scf_len){
				$secon_end_scf=$key_st_ed_ovl_2;$secon_end_scf_len=$st_ed_ovl{"2"}{$key_st_ed_ovl_2}[0];
			}
		}
		}
		#mkdir sure the middle block
		my $first1;my $first2;my $end1;my $end2;
		$first1=$first_end_scf;$end1=$secon_end_scf;
		if($st_ed_ovl{"1"}{$first_end_scf}){
			#first
			$first2=$rel_ovl{$first_end_scf};
			#if($rel_ovl_reverse{$first_end_scf}){
			#	($first1,$first2)=($first2,$first1);
			#}
		}
		if($rel_ovl{$first1}){
			$first2=$rel_ovl{$first1};
			#if($rel_ovl_reverse{$first1}){
			#	($first1,$first2)=($first2,$first1);
			#}
		}
		if($st_ed_ovl{"2"}{$secon_end_scf}){
			#second
			$end2=$st_ed_ovl{"2"}{$secon_end_scf}[1];
		}
		if($rel_ovl{$end1}){
			$end2=$rel_ovl{$end1};
			#if($rel_ovl_reverse{$end1}){
			#	($end1,$end2)=($end2,$end1);
			#}
		}
		my @bac_key1_sort;my @bac_key2_sort;
		my @bac_key1_filter=grep {!($_ eq $first1)} @bac_key1;
		my @bac_key2_filter=grep {!($_ eq $end1)} @bac_key2;
		my @keys_ovl_filter=keys $st_ed_ovl{"1"};#only deal with not contained overlap
		if(@bac_key1_filter != 0){
			#%scf_len  sort the scaffold based on the length of the scaffold
			#remember the overlap include in the contain block
			foreach my $bac_key1_filter_key(sort {$scf_len{$b}<=>$scf_len{$a}} @bac_key1_filter){
				push @bac_key1_sort,$bac_key1_filter_key;
			}
		}
		if(@bac_key2_filter != 0){
			foreach my $bac_key2_filter_key(sort {$scf_len{$b}<=>$scf_len{$a}} @bac_key2_filter){
				push @bac_key2_sort,$bac_key2_filter_key;
			}
		}
		#output sort; %out $sort
		my $cal_sort=1;
		if($first1){
			#$out->{$sort}{$cal_sort}[0]=[$first1];print "result\t$first1\tfirst1\n";
			my @first1;$first1[0]=$first1;
			my @first2;my @end2;
			if($st_ed{$first1}){
				push @first1,$st_ed{$first1}[0],$st_ed{$first1}[-2],$st_ed{$first1}[-1],"con1";
			}
			if($st_ed_ovl{"1"}{$first1}){
				push @first1,$ovl_blo{$first1}{$first2}[-2],$ovl_blo{$first1}{$first2}[1],$ovl_blo{$first1}{$first2}[2],"ovl";
				#push @first2,$first2,$ovl_blo{$first1}{$first2}[-1],$ovl_blo{$first1}{$first2}[3],$ovl_blo{$first1}{$first2}[4];
			}
			if($rel_ovl{$first1} && ($end1 eq $rel_ovl{$first1})){
				##$out->{$sort}{$cal_sort}[1]=[$end1];
				###print "$end1\tovl_first1\t";
				#push @first2,$first2,$ovl_blo{$first1}{$first2}[-1],$ovl_blo{$first1}{$first2}[3],$ovl_blo{$first1}{$first2}[4];
				$end1=0;
			}
			if(@first1){$out->{$sort}{$cal_sort}[0]=[@first1];}
			#if(@first2){$out->{$sort}{$cal_sort}[1]=[@first2];}
		}
		if($first1 && $first2){
			#push @{$out->{$sort}{$cal_sort}},$first2;print "result\t$first2\tfirst2\n";
			my @first2;
			push @first2,$first2,$ovl_blo{$first1}{$first2}[-1],$ovl_blo{$first1}{$first2}[3],$ovl_blo{$first1}{$first2}[4],"ovl";
			if(@first2){$out->{$sort}{$cal_sort}[1]=[@first2];}
			if($rel_add{$first1} or $rel_add{$first2}){
				$cal_sort++;
				my @add1;my @add2;
				my $tkey1;my $tkey2;
				if($rel_add{$first1}){$tkey1=$first1;$tkey2=$rel_add{$first1};}
				if($rel_add{$first2}){$tkey1=$first2;$tkey2=$rel_add{$first2};}
				my @tt1;my @tt2;
				push @tt1,$tkey1,$ovl_blo{$tkey1}{$tkey2}[-2],$ovl_blo{$tkey1}{$tkey2}[1],$ovl_blo{$tkey1}{$tkey2}[2],"ovllf";
				push @tt2,$tkey2,$ovl_blo{$tkey1}{$tkey2}[-1],$ovl_blo{$tkey1}{$tkey2}[3],$ovl_blo{$tkey1}{$tkey2}[4],"ovllf";
				if(@tt1){$out->{$sort}{$cal_sort}[0]=[@tt1];}
				if(@tt2){$out->{$sort}{$cal_sort}[1]=[@tt2];}
			}
		}
		$cal_sort++;
		foreach my $ele(@bac_key1_sort){
			#$out->{$sort}{$cal_sort}=[$ele];print "result\t$ele\tele1\n";
			#if($rel_ovl{$ele}){$out->{$sort}{$cal_sort}[1]=$rel_ovl{$ele};print "result\t$rel_ovl{$ele}\tele1_ol\n";}
			my @ele;
			push @ele,$ele,$st_ed{$ele}[0],$st_ed{$ele}[-2],$st_ed{$ele}[-1],"con2";
			my @ele_ovl;my $ele_ovl;
			if($rel_ovl{$ele}){
				$ele_ovl=$rel_ovl{$ele};
				@ele=();
				push @ele,$ele,$ovl_blo{$ele}{$ele_ovl}[-2],$ovl_blo{$ele}{$ele_ovl}[1],$ovl_blo{$ele}{$ele_ovl}[2],"ovl";
				push @ele_ovl,$ele_ovl,$ovl_blo{$ele}{$ele_ovl}[-1],$ovl_blo{$ele}{$ele_ovl}[3],$ovl_blo{$ele}{$ele_ovl}[4],"ovl";
			}
			if(@ele){$out->{$sort}{$cal_sort}[0]=[@ele];}
			if(@ele_ovl){$out->{$sort}{$cal_sort}[1]=[@ele_ovl];}
			$cal_sort++;
			if($rel_add{$ele} or $rel_add{$ele_ovl}){
				my @add1;my @add2;
				my $tkey1;my $tkey2;
				if($rel_add{$ele}){$tkey1=$ele;$tkey2=$rel_add{$ele};}
				if($rel_add{$ele_ovl}){$tkey1=$ele_ovl;$tkey2=$rel_add{$ele_ovl};}
				my @tt1;my @tt2;
				push @tt1,$tkey1,$ovl_blo{$tkey1}{$tkey2}[-2],$ovl_blo{$tkey1}{$tkey2}[1],$ovl_blo{$tkey1}{$tkey2}[2],"ovll1";
				push @tt2,$tkey2,$ovl_blo{$tkey1}{$tkey2}[-1],$ovl_blo{$tkey1}{$tkey2}[3],$ovl_blo{$tkey1}{$tkey2}[4],"ovll1";
				if(@tt1){$out->{$sort}{$cal_sort}[0]=[@tt1];}
				if(@tt2){$out->{$sort}{$cal_sort}[1]=[@tt2];}
				$cal_sort++;
			}
		}
		foreach my $ele(@keys_ovl_filter){
			if($ele eq $first1 || $ele eq $end2){next;}
			if($ele eq $end1 || $ele eq $first2){next;}
			#$out->{$sort}{$cal_sort}=[$ele,$rel_ovl{$ele}];print "result\t$ele\t$rel_ovl{$ele}\tovl\n";
			my @ele;my @ele_ovl;my $ele_ovl=$rel_ovl{$ele};
			push @ele,$ele,$ovl_blo{$ele}{$ele_ovl}[-2],$ovl_blo{$ele}{$ele_ovl}[1],$ovl_blo{$ele}{$ele_ovl}[2],"ovl";
			push @ele_ovl,$ele_ovl,$ovl_blo{$ele}{$ele_ovl}[-1],$ovl_blo{$ele}{$ele_ovl}[3],$ovl_blo{$ele}{$ele_ovl}[4],"ovl";
			if(@ele){$out->{$sort}{$cal_sort}[0]=[@ele];}
			if(@ele_ovl){$out->{$sort}{$cal_sort}[1]=[@ele_ovl];}
			$cal_sort++;
			if($rel_add{$ele} or $rel_add{$ele_ovl}){
				my @add1;my @add2;
				my $tkey1;my $tkey2;
				if($rel_add{$ele}){$tkey1=$ele;$tkey2=$rel_add{$ele};}
				if($rel_add{$ele_ovl}){$tkey1=$ele_ovl;$tkey2=$rel_add{$ele_ovl};}
				my @tt1;my @tt2;
				push @tt1,$tkey1,$ovl_blo{$tkey1}{$tkey2}[-2],$ovl_blo{$tkey1}{$tkey2}[1],$ovl_blo{$tkey1}{$tkey2}[2],"ovll2";
				push @tt2,$tkey2,$ovl_blo{$tkey1}{$tkey2}[-1],$ovl_blo{$tkey1}{$tkey2}[3],$ovl_blo{$tkey1}{$tkey2}[4],"ovll2";
				if(@tt1){$out->{$sort}{$cal_sort}[0]=[@tt1];}
				if(@tt2){$out->{$sort}{$cal_sort}[1]=[@tt2];}
				$cal_sort++;
			}
		}
		foreach my $ele(@bac_key2_sort){
			#$out->{$sort}{$cal_sort}=[$ele];print "result\t$ele\tele2\n";
			#if($rel_ovl{$ele}){$out->{$sort}{$cal_sort}[1]=$rel_ovl{$ele};print "result\t$rel_ovl{$ele}\tele2_ol\n";}
			my @ele;
			push @ele,$ele,$st_ed{$ele}[0],$st_ed{$ele}[-2],$st_ed{$ele}[-1],"con3";
			my @ele_ovl;my $ele_ovl;
			if($rel_ovl{$ele}){
				$ele_ovl=$rel_ovl{$ele};
				@ele=();
				push @ele,$ele,$ovl_blo{$ele_ovl}{$ele}[-1],$ovl_blo{$ele_ovl}{$ele}[3],$ovl_blo{$ele_ovl}{$ele}[4],"ovl";
				push @ele_ovl,$ele_ovl,$ovl_blo{$ele_ovl}{$ele}[-2],$ovl_blo{$ele_ovl}{$ele}[1],$ovl_blo{$ele_ovl}{$ele}[2],"ovl";
			}
			if(@ele){$out->{$sort}{$cal_sort}[1]=[@ele];}
			if(@ele_ovl){$out->{$sort}{$cal_sort}[0]=[@ele_ovl];}
			$cal_sort++;
			if($rel_add{$ele} or $rel_add{$ele_ovl}){
				my @add1;my @add2;
				my $tkey1;my $tkey2;
				if($rel_add{$ele}){$tkey1=$ele;$tkey2=$rel_add{$ele};}
				if($rel_add{$ele_ovl}){$tkey1=$ele_ovl;$tkey2=$rel_add{$ele_ovl};}
				my @tt1;my @tt2;
				push @tt1,$tkey1,$ovl_blo{$tkey1}{$tkey2}[-2],$ovl_blo{$tkey1}{$tkey2}[1],$ovl_blo{$tkey1}{$tkey2}[2],"ovll3";
				push @tt2,$tkey2,$ovl_blo{$tkey1}{$tkey2}[-1],$ovl_blo{$tkey1}{$tkey2}[3],$ovl_blo{$tkey1}{$tkey2}[4],"ovll3";
				if(@tt1){$out->{$sort}{$cal_sort}[0]=[@tt1];}
				if(@tt2){$out->{$sort}{$cal_sort}[1]=[@tt2];}
				$cal_sort++;
			}
		}
		if($end1){
			#$out->{$sort}{$cal_sort}=[$end1];print "result\t$end1\tend1\n";
			my @end1;my @end2;
			if($st_ed{$end1}){
				push @end1,$end1,$st_ed{$end1}[0],$st_ed{$end1}[-2],$st_ed{$end1}[-1],"con4";
			}
			if($st_ed_ovl{"2"}{$end1}){
				#$end2=$st_ed_ovl{"2"}{$secon_end_scf}[1];
				push @end1,$end1,$ovl_blo{$end2}{$end1}[-1],$ovl_blo{$end2}{$end1}[3],$ovl_blo{$end2}{$end1}[4],"ovl";
				push @end2,$end2,$ovl_blo{$end2}{$end1}[-2],$ovl_blo{$end2}{$end1}[1],$ovl_blo{$end2}{$end1}[2],"ovl";
			}
			if(@end1){$out->{$sort}{$cal_sort}[1]=[@end1];}
			if(@end2){$out->{$sort}{$cal_sort}[0]=[@end2];}
			if($rel_add{$end1} or $rel_add{$end2}){
				$cal_sort++;
				my @add1;my @add2;
				my $tkey1;my $tkey2;
				if($rel_add{$end1}){$tkey1=$end1;$tkey2=$rel_add{$end1};}
				if($rel_add{$end2}){$tkey1=$end2;$tkey2=$rel_add{$end2};}
				my @tt1;my @tt2;
				push @tt1,$tkey1,$ovl_blo{$tkey1}{$tkey2}[-2],$ovl_blo{$tkey1}{$tkey2}[1],$ovl_blo{$tkey1}{$tkey2}[2],"ovll4";
				push @tt2,$tkey2,$ovl_blo{$tkey1}{$tkey2}[-1],$ovl_blo{$tkey1}{$tkey2}[3],$ovl_blo{$tkey1}{$tkey2}[4],"ovll4";
				if(@tt1){$out->{$sort}{$cal_sort}[0]=[@tt1];}
				if(@tt2){$out->{$sort}{$cal_sort}[1]=[@tt2];}
			}
		}
		if($end1 && $end2){
			##push @{$out->{$sort}{$cal_sort}},$end2;print "result\t$end2\tend2\n";
			#push @end2,$end2,$ovl_blo{$end2}{$end1}[-2],$ovl_blo{$end2}{$end1}[1],$ovl_blo{$end2}{$end1}[2];
			#if(@end2){$out->{$sort}{$cal_sort}[1]=[@end2];}
		}
		$cal_sort++;
		#foreach my $sort_c(sort {$a<=>$b} keys %{$out}){
		#	foreach my $sort_c1(sort {$a<=>$b} keys $out->{$sort_c}){
		#		if($out->{$sort_c}{$sort_c1}[0] and $out->{$sort_c}{$sort_c1}[1] and $out->{$sort_c}{$sort_c1}[0][-1] eq 'ovl' and $out->{$sort_c}{$sort_c1}[1][-1] eq 'ovl'){
		#			if($out->{$sort_c}{$sort_c1}[0][0] =~/$bac2_id/ and $out->{$sort_c}{$sort_c1}[1][0] =~/$bac1_id/){
		#				#reverse record with rel_ovl relation
		#				($out->{$sort_c}{$sort_c1}[0],$out->{$sort_c}{$sort_c1}[1])=($out->{$sort_c}{$sort_c1}[1],$out->{$sort_c}{$sort_c1}[0]);
		#			}
		#		}
		#	}
		#}
		print "checking finish\n\n";
		return \%scf_len;
	}else{
		print STDERR "no record\n";
		$scf_len{"no"}="record";
		return \%scf_len;
	}
}

sub check_ovl_stat{
	my ($st1,$ed1,$st2,$ed2)=@_;
	my %tmp;
	$tmp{"1"}=$st1;$tmp{"2"}=$ed1;$tmp{"3"}=$st2;$tmp{"4"}=$ed2;
	my @tmp=sort {$tmp{$a}<=>$tmp{$b}} keys %tmp;my $tmp=join "",@tmp;
	if($tmp=~/^12..$/ or $tmp=~/^34..$/ or $tmp=~/^21..$/ or $tmp=~/^43..$/){return 0;}
	else{return 1;}
}






1
