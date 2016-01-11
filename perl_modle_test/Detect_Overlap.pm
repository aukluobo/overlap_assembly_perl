package Detect_Overlap;

#huangguodong@genomics.cn 

sub detect{
    my ($self,$m8r,$scf_len,$bacSort,$dev)=@_;
    foreach my $bac1(keys %{$m8r}){
		foreach my $bac2(keys $m8r->{$bac1}){
		    my @key=sort {$m8r->{$bac1}{$bac2}{$a}[0]<=>$m8r->{$bac1}{$bac2}{$b}[0]} keys $m8r->{$bac1}{$bac2};
		    #my $bac1_len=$scf_len{$bac1};my $bac2_len=$scf_len{$bac2};
		    my $bacID1=$bac1;$bacID1=~s/_\d+$//;my $bacID1old=$bacID1;
		    my $bacID2=$bac2;$bacID2=~s/_\d+$//;
		    ($bacID1,$bacID2)=sort {$bacSort->{$a}<=>$bacSort->{$b}} ($bacID1,$bacID2);
		    my $bac1_p=(grep /$bacID1\_/,($bac1,$bac2))[0];my $bac2_p=(grep /$bacID2\_/,($bac1,$bac2))[0];
		    my $bac1_len=$scf_len->{$bac1_p};my $bac2_len=$scf_len->{$bac2_p};
		    if(@key==1){
				#my $bac1_len=$scf_len{$bac1};my $bac2_len=$scf_len{$bac2};
				my $st1=$m8r->{$bac1}{$bac2}{$key[0]}[0];
				my $ed1=$m8r->{$bac1}{$bac2}{$key[0]}[1];
				my $st2=$m8r->{$bac1}{$bac2}{$key[0]}[2];
				my $ed2=$m8r->{$bac1}{$bac2}{$key[0]}[3];
				my $type=$m8r->{$bac1}{$bac2}{$key[0]}[4];
				#if($bacID1old eq $bacID2){($st1,$ed1,$st2,$ed2)=($st2,$ed2,$st1,$ed1);}
				my $jia=$type=~s/\+/+/g;my $jian=$type=~s/\-/-/g;
				my $tp;
				if($jia>$jian){$tp="+";}else{$tp="-";}
				if($bacID1old eq $bacID2){($st1,$ed1,$st2,$ed2)=($st2,$ed2,$st1,$ed1);$tp=~tr/\+\-/\-\+/;}
				#contain one bac scaffold
				my $con=contain($st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len);
				if($con){
				    if($con==3){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"C $bac1_p|C $bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
				    elsif($con==2){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"$bac1_p|C $bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
				    elsif($con==1){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"C $bac1_p|$bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
				    #push @{$dev->{$bacID1}{$bacID2}{"contain"}},"$bac1_p|$bac2_p",$tp;
				    ##$dev->{"$bac1|$bac2"}{"contain"}=[1,$tp];
				    next;
				}
				#overlap
				#if the map length small than 200 then discard the record
				if($ed1-$st1+1<200 || $ed2-$st2+1<200){print "overlaptooshort: $bac1_p\t$bac2_p\t$st1\t$ed1\t$st2\t$ed2\t$bac1_len\t$bac2_len\t$tp\n";next;}
				my $ove=0;$ove=overlap($st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len,$tp);
				if($ove){
				    #$dev->{"$bac1|$bac2"}{"overlap"}=[1,$tp];
				    push @{$dev->{$bacID1}{$bacID2}{"overlap"}},"$bac1_p|$bac2_p|$ove",$tp;
				    push @{$dev->{$bacID1}{$bacID2}{"overlap"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    next;
				}
				#$dev->{"$bac1|$bac2"}{"discard"}=[1,$tp];
				push @{$dev->{$bacID1}{$bacID2}{"discard"}},"$bac1_p|$bac2_p",$tp;
				push @{$dev->{$bacID1}{$bacID2}{"discard"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				next;
		    }
		    else{
				my ($st1,$ed1,$st2,$ed2);
				my $type;
				my $Used1=0;my $Used2=0;my $ToLen1=0;my $ToLen2=0;
				my @st1;my @ed1;my @st2;my @ed2;
				foreach my $rec(@key){
				    push @st1,$m8r->{$bac1}{$bac2}{$rec}[0];push @ed1,$m8r->{$bac1}{$bac2}{$rec}[1];
				    push @st2,$m8r->{$bac1}{$bac2}{$rec}[2];push @ed2,$m8r->{$bac1}{$bac2}{$rec}[3];
				    unless($st1){
						$st1=$m8r->{$bac1}{$bac2}{$rec}[0];
						$ed1=$m8r->{$bac1}{$bac2}{$rec}[1];
						$st2=$m8r->{$bac1}{$bac2}{$rec}[2];
						$ed2=$m8r->{$bac1}{$bac2}{$rec}[3];
						$type.=$m8r->{$bac1}{$bac2}{$rec}[4];
						next;
				    }
				    if($m8r->{$bac1}{$bac2}{$rec}[0]<$st1){$st1=$m8r->{$bac1}{$bac2}{$rec}[0];}
				    if($m8r->{$bac1}{$bac2}{$rec}[1]>$ed1){$ed1=$m8r->{$bac1}{$bac2}{$rec}[1];}
				    if($m8r->{$bac1}{$bac2}{$rec}[2]<$st2){$st2=$m8r->{$bac1}{$bac2}{$rec}[2];}
				    if($m8r->{$bac1}{$bac2}{$rec}[3]>$ed2){$ed2=$m8r->{$bac1}{$bac2}{$rec}[3];}
				    $type.=$m8r->{$bac1}{$bac2}{$rec}[4];
				}
				my $ctrl_w=1;
				#while(){}
				my %ToLen1;my %ToLen2;
				if($bacID1old eq $bacID2){
				    foreach my $nlu(($st1..$ed1)){$ToLen2{$nlu}=1;}
				    foreach my $nlu(($st2..$ed2)){$ToLen1{$nlu}=1;}
				}
				else{
				    foreach my $nlu(($st1..$ed1)){$ToLen1{$nlu}=1;}
				    foreach my $nlu(($st2..$ed2)){$ToLen2{$nlu}=1;}
				}
				foreach my $rec(@key){
				    foreach my $nlu(($m8r->{$bac1}{$bac2}{$rec}[0]..$m8r->{$bac1}{$bac2}{$rec}[1])){
					if($ToLen1{$nlu}){undef $ToLen1{$nlu};delete $ToLen1{$nlu};}
				    }
				    foreach my $nlu(($m8r->{$bac1}{$bac2}{$rec}[2]..$m8r->{$bac1}{$bac2}{$rec}[3])){
					if($ToLen2{$nlu}){undef $ToLen2{$nlu};delete $ToLen2{$nlu};}
				    }
				}
				my @UnUsed1=keys %ToLen1;my @UnUsed2=keys %ToLen2;
				my $UnUsed1=@UnUsed1;my $UnUsed2=@UnUsed2;
				my $DvUn1=$UnUsed1/($ed1-$st1+1)*100;
				my $DvUn2=$UnUsed2/($ed2-$st2+1)*100;
				if($DvUn1>=50 or $DvUn2>=50){
				    #print STDERR "Unmap bases are too much for $bac1_p $bac2_p\n";next;}#need to fix
				    #discard the biggest and find the second biggest
				    print STDERR "Unmap bases:\t$bac1_p\t$bac2_p\t$UnUsed1\t$bac1_len\t$UnUsed2\t$bac2_len\n";
				}

				my $jia=$type=~s/\+/+/g;my $jian=$type=~s/\-/-/g;
				my $tp;
				if($jia>$jian){$tp="+";}else{$tp="-";}
				if($bacID1old eq $bacID2){($st1,$ed1,$st2,$ed2)=($st2,$ed2,$st1,$ed1);$tp=~tr/\+\-/\-\+/;}
				#contain one bac scaffold
				my $con=contain($st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len);
	            if($con){
				    if($con==3){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"C $bac1_p|C $bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
				    elsif($con==2){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"$bac1_p|C $bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
				    elsif($con==1){
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},"C $bac1_p|$bac2_p",$tp;
						push @{$dev->{$bacID1}{$bacID2}{"contain"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				    }
		            ##$dev->{"$bac1|$bac2"}{"contain"}=[1,$tp];
				    #push @{$dev->{$bacID1}{$bacID2}{"contain"}},"$bac1_p|$bac2_p",$tp;
				    next;
	            }
	            #overlap
	            my $ove=0;$ove=overlap($st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len,$tp);
	            if($ove){
	                #$dev->{"$bac1|$bac2"}{"overlap"}=[1,$tp];
			    	push @{$dev->{$bacID1}{$bacID2}{"overlap"}},"$bac1_p|$bac2_p|$ove",$tp;
			    	push @{$dev->{$bacID1}{$bacID2}{"overlap"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
			    	next;
	            }
	            #$dev->{"$bac1|$bac2"}{"discard"}=[1,$tp];
				push @{$dev->{$bacID1}{$bacID2}{"discard"}},"$bac1_p|$bac2_p",$tp;
				push @{$dev->{$bacID1}{$bacID2}{"discard"}},$st1,$ed1,$st2,$ed2,$bac1_len,$bac2_len;
				next;
		    }
		}
    }
}

sub overlap{
    my ($st1,$ed1,$st2,$ed2,$len1,$len2,$tp)=@_;
    my $cut1=20;my $cut2=500;
    my $slen=(sort {$a<=>$b} ($len1,$len2))[0];
    if($slen>10000){$cut1=100;}
    if($slen>5000){$cut1=50;}
    #if($len1){}
    if($st1<=$cut1){
		if($st2<=$cut1){
		    if($tp=~/\-/){return "st|st";}
		    if($tp=~/\+/){
				if($len2-$ed2<=$cut2 || $len1-$ed1<$cut2){return "st|st";}
		    }
		}
		elsif($len2-$ed2<=$cut1){
		    if($tp=~/\+/){return "st|ed";}
		    if($tp=~/\-/){
				if($len1-$ed1<$cut2){return "st|ed";}
		    }
		}
    }
    elsif($len1-$ed1<=$cut1){
		if($st2<=$cut1){
		    if($tp=~/\+/){return "ed|st";}
		    if($tp=~/\-/){
				if($len2-$ed2<$cut2 || $st1<=$cut2){return "ed|st";}
		    }
		}
		elsif($len2-$ed2<=$cut1){
		    if($tp=~/\+/){
				if($st1<$cut2 || $st2<$cut2){return "ed|ed";}
		    }
		    if($tp=~/\-/){return "ed|ed";}
		}
    }
}
sub contain{
    my ($st1,$ed1,$st2,$ed2,$len1,$len2)=@_;
    my $cut1=20;
    my $slen=(sort {$a<=>$b} ($len1,$len2))[0];
    if($slen>10000){$cut1=500;}
    if($slen>5000){$cut1=200;}
    my $type1=0;my $type2=0;
    my $con1=0;my $con2=0;
    if($st1<=$cut1){
	if($len1-$ed1<=$cut1){
	    #bac1 was contained;
	    $type1=1;$con1=1;
	}
    }
    if($st2<=$cut1 && $len2-$ed2<=$cut1){
	#bac2 was contained
	$type2=1;$con2=1;
    }
    if($type1 && $type2){return 3;}
    elsif($type1){return 1;}
    elsif($type2){return 2;}
    else{return 0;}
}

1



