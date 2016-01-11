package Filter_Overlap_BAC_1_forsimple;
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
	my ($self,$m8,$hash,$del)=@_;
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
		#$p_id1=~s/_\d+$//;
		my $p_id2=$aa[1];if($del->{$p_id2}){next;}
		#$p_id2=~s/_\d+$//;
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
			my @pre=();my $pre=0;
		  	my %m_tmp;
		  	foreach my $key2(sort {$tmp{$key}{$key1}{$a}[0]<=>$tmp{$key}{$key1}{$b}[0]} keys $tmp{$key}{$key1}){
			    my $st1=$tmp{$key}{$key1}{$key2}[0];my $ed1=$tmp{$key}{$key1}{$key2}[1];
			    my $st2=$tmp{$key}{$key1}{$key2}[2];my $ed2=$tmp{$key}{$key1}{$key2}[3];
				my $type=$tmp{$key}{$key1}{$key2}[4];
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
					#noverlap
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
					#neibor
					my $len1=$pre[0]-$ed1+1;my $len2=$st1-$pre[1]+1;
				  	my $len;if($len1>0){$len=$len1;}else{$len=$len2;}
				  	if($len>500){
						@{$m_tmp{$pre}}=@pre;
						##$m_tmp{$key2}=[$st1,$ed1,$st2,$ed2];
						##@pre=();
						@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;
						next;
				  	}
				  	#if($c_b==2){
				  	#	#contain
				  	#	if($pre[2] != $st2){
					#    	if($k_b=~/^1/){
					#    		#
					#    		next;
					#    	}
					#		if($k_b=~/^3/){@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;next;}
					#	}else{
					#		if($pre[3]<$ed2){
					#			@pre=($st1,$ed1,$st2,$ed2,$type);$pre=$key2;next;
					#		}else{
					#			next;
					#		}
					#	}
				  	#}
				  	#
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
		  	foreach my $ass(sort {$a<=>$b} keys %m_tmp){
		    	my $out=join "\t",@{$m_tmp{$ass}};
				#print "$key\t$key1\t$ass\t$out\n";
				#print "$key1\t$key\t90\t100\t0\t0\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t0.0\t500\n";
				$hash->{$key}{$key1}{$num}=[$m_tmp{$ass}[0],$m_tmp{$ass}[1],$m_tmp{$ass}[2],$m_tmp{$ass}[3],$m_tmp{$ass}[4]];
				print STDERR "filter1\t$key\t$key1\t$num\t$m_tmp{$ass}[0]\t$m_tmp{$ass}[1]\t$m_tmp{$ass}[2]\t$m_tmp{$ass}[3]\t$m_tmp{$ass}[4]\n";
				$num++;
				#
		  	}
		 	next;
		  	@pre=();
		 	#check target2 contain
		  	my %con;my $t_cou=0;my @last=sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} keys %m_tmp;
		  	if(@last>1){
		    	my $first=shift @last;my $last=pop @last;
				if(@last){
			  	my $fr_check=check_contain($first,\@last,\%m_tmp)
				}
		  		foreach my $ass(sort {$m_tmp{$a}[2]<=>$m_tmp{$b}[2]} keys %m_tmp){
		    	my $c_st=$m_tmp{$ass}[2];my $c_ed=$m_tmp{$ass}[3];
				unless(@pre){@pre=($c_st,$c_ed);$t_cou=1;$pre=$ass;next;}
				my %ch_tmp;
				$ch_tmp{"1"}=$pre[0];$ch_tmp{"2"}=$pre[1];
				$ch_tmp{"3"}=$c_st;$ch_tmp{"4"}=$c_ed;
				my @k_ch=sort {$ch_tmp{$a}<=>$ch_tmp{$b}} keys %ch_tmp;
				my $k_ch=join "",@k_ch;
				if($k_ch=~/^1..2$/ or $k_ch=~/^3..4$/){
			 		my $len1=$m_tmp{$ass}[0]-$m_tmp{$pre}[1]+1;my $len2=$m_tmp{$pre}[0]-$m_tmp{$ass}[1]+1;
				}
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









1;
