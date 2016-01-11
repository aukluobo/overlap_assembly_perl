use strict;
use warnings;

#huangguodong@genomics.cn 20140928

die "perl $0 <map1> <map2> <map1 cM_sort:1/-1> <map2 cM_sort:1/-1> <split bock or not:1/0> <build:1/0>" if(@ARGV==0);

my ($map1,$map2,$m1_sort,$m2_sort,$split,$build)=@ARGV;
my (%map1,%map2);
read_map($map1,\%map1);read_map($map2,\%map2);

if($split){
	#foreach(keys %map1){print "check\t$_\t$map1{$_}[0]\t$map1{$_}[1]\n";}
	resort(\%map1,\%map2,$m1_sort,$m2_sort);
	#print "resort 2\n";
	resort(\%map2,\%map1,$m2_sort,$m1_sort);
	#delete uncorrect link
	del_link(\%map1,\%map2);
	#merge neibor link
	my %block;
	merge(\%map1,\%map2,\%block);
	foreach my $key(sort {$a<=>$b} keys %block){
		print "block$key\t";
		foreach(@{$block{$key}}){print "$_\t";}
		print "\n";
	}
	foreach(sort {$map1{$a}[0]<=>$map1{$b}[0]} keys %map1){
		my $out=join "\t",@{$map1{$_}};print "map1\t$_\t$out\n";
	}
	foreach(sort {$map2{$a}[0]<=>$map2{$b}[0]} keys %map2){
		my $out=join "\t",@{$map2{$_}};print "map2\t$_\t$out\n";
	}
	#find cM link information
	#cM_link(\%map1,\%map2);
	#split block to make the program faster
	my @split;
	@split=split_block(\%map1,\%map2);
	my $num=1;
	my $s_sta=0;my $s_stb=0;
	open AA,">split_aa_$num.lst" or die;
	open BB,">split_bb_$num.lst" or die;
	my @aa=sort {$map1{$a}[0]<=>$map1{$b}[0]} keys %map1;
	my @bb=sort {$map2{$a}[0]<=>$map2{$b}[0]} keys %map2;
	foreach my $s_ctg(@split){
		for(my $i=$s_sta;$i<@aa;$i++){
			my @out=@{$map1{$aa[$i]}};shift @out;my $out=join "\t",@out;
			print AA "$aa[$i]\t$out\n";
			if($aa[$i] eq $s_ctg){
				$s_sta=$i+1;
				close AA;
				last;
			}
		}
		for(my $i=$s_stb;$i<@bb;$i++){
			my @out=@{$map2{$bb[$i]}};shift @out;my $out=join "\t",@out;
			print BB "$bb[$i]\t$out\n";
			if($bb[$i] eq $s_ctg){
				$s_stb=$i+1;
				close BB;
				last;
			}
		}
		$num++;
		open AA,">split_aa_$num.lst" or die;
		open BB,">split_bb_$num.lst" or die;
	}
	#get last part
	for(my $i=$s_sta;$i<@aa;$i++){
		my @out=@{$map1{$aa[$i]}};shift @out;my $out=join "\t",@out;print AA "$aa[$i]\t$out\n";
	}
	for(my $i=$s_stb;$i<@bb;$i++){
		my @out=@{$map2{$bb[$i]}};shift @out;my $out=join "\t",@out;print BB "$bb[$i]\t$out\n";
	}
	close AA;close BB;
}else{
	my %map1_rv;my %map2_rv;
	my @key1=sort {$map1{$a}[0]<=>$map1{$b}[0]} keys %map1;
	my @key2=sort {$map2{$a}[0]<=>$map2{$b}[0]} keys %map2;
	my $end1=$key1[-1];my $end2=$key2[-1];
	#if($end1 eq $end2){
	#	pop @key2;
	#	undef $map2{$end2};delete $map2{$end2};
	#	$end2=$key2[-1];
	#}
	my $aa_st=$key1[0];my $bb_st=$key2[0];
	my @same;
	for(my $i=0;$i<@key1;$i++){
		my $id=$key1[$i];
		if($map2{$id}){push @same,$id;}
		my $next_index;
		if($key1[$i+1]){
			$next_index=$map1{$key1[$i+1]}[0];
		}else{
			$next_index=-1;
		}
		$map1_rv{$map1{$id}[0]}=[$id,$next_index];
		my $out=join "\t",@{$map1{$id}};
		print "map1\t$id\t$out\n";
	}
	for(my $i=0;$i<@key2;$i++){
		my $id=$key2[$i];
		my $next_index;
		if($key2[$i+1]){
			$next_index=$map2{$key2[$i+1]}[0];
		}else{
			$next_index=-1;
		}
		$map2_rv{$map2{$id}[0]}=[$id,$next_index];
		my $out=join "\t",@{$map2{$id}};
		print "map2\t$id\t$out\n";
	}

	my $path="start";my %used;my $shortest=@same;
	if($build){build(\%map1,\%map2,\%map1_rv,\%map2_rv,\%used,$end1,$end2,$aa_st,$bb_st,$path,$shortest);}
}


#read map
sub read_map{
	my ($map,$hash)=@_;
	open IN,$map or die;
	my $num=1;
	while(<IN>){
		my @aa=split;
		#$hash->{$aa[0]}=[$num];
		my $id=shift @aa;push @{$hash->{$id}},$num,@aa;
		$num++;
		#if($aa[3]){push @{$hash->{$aa[0]}},$aa[3];}
		#if($aa[5]){push @{$hash->{$aa[0]}},$aa[5];}
	}
	close IN;
}
#resort itself based on each other
sub resort{
	my ($aa,$bb,$aa_s,$bb_s)=@_;
	my @key=sort {$aa->{$a}[0]<=>$aa->{$b}[0]} keys %{$aa};
	my @key1;
	foreach my $key(@key){
		if($bb->{$key}){push @key1,$key;}
	}
	for(my $i=1;$i<@key1;$i++){
		my $pre=$i-1;$pre=$key1[$pre];my $aft=$key1[$i];
		my $aa_cou_pre=@{$aa->{$pre}};my $aa_cou_aft=@{$aa->{$aft}};
		if($aa->{$pre}[1] == $aa->{$aft}[1] && $aa_cou_pre == $aa_cou_aft){
			my $bb_cou_pre=@{$bb->{$pre}};my $bb_cou_aft=@{$bb->{$aft}};
			if($bb->{$pre}[1] != $bb->{$aft}[1] or $bb_cou_pre != $bb_cou_aft){
				#can sort
				my %tmp;
				if($bb->{$pre}[1] != $bb->{$aft}[1]){
					$tmp{$pre}=$bb->{$pre}[1];$tmp{$aft}=$bb->{$aft}[1];
				}else{
					$tmp{$pre}=$bb_cou_pre;$tmp{$aft}=$bb_cou_aft;
				}
				my @t_key=sort {$tmp{$a}<=>$tmp{$b}} keys %tmp;
				my $pre_nn=$aa->{$t_key[0]}[0];
				my $aft_nn=$aa->{$t_key[1]}[0];
				$aa->{$pre}[0]=$pre_nn;
				$aa->{$aft}[0]=$aft_nn;
			}
		}
	}
	#block sort
	my %all;
	foreach my $ctg(@key){
		if($bb->{$ctg}){
			my $bbcm=$bb->{$ctg}[1];
			my $cou=@{$aa->{$ctg}};
			push @{$all{$aa->{$ctg}[1]}{$cou}{$bbcm}},$ctg;
		}
	}
	foreach my $aa_cm(keys %all){
		foreach my $cou(sort {$a<=>$b} keys $all{$aa_cm}){
			my %aa_pl;my %aa_re;
			my @bbcm;
			if($bb_s>0){
				@bbcm=sort {$a<=>$b} keys $all{$aa_cm}{$cou};
			}else{
				@bbcm=sort {$b<=>$a} keys $all{$aa_cm}{$cou};
			}
			foreach my $pctg(@bbcm){
				my @aa_ctg=@{$all{$aa_cm}{$cou}{$pctg}};
				foreach my $aa_ctg(@aa_ctg){
					$aa_pl{$aa_ctg}=$aa->{$aa_ctg}[0];
					$aa_re{$aa_ctg}=[@{$aa->{$aa_ctg}}];
				}
			}
			my @aa_pl_ctg=sort {$aa_pl{$a}<=>$aa_pl{$b}} keys %aa_pl;
			for(my $i=0;$i<@aa_pl_ctg-1;$i++){
				if($aa_pl{$aa_pl_ctg[$i+1]}-$aa_pl{$aa_pl_ctg[$i]}>1){
					print "$aa_pl_ctg[$i+1]\t$aa_pl{$aa_pl_ctg[$i+1]}\t$aa_pl_ctg[$i]\t$aa_pl{$aa_pl_ctg[$i]}\tuncorrect\n";
				}
			}
			my $count=0;
			foreach my $pctg(@bbcm){
				my @pctg=@{$all{$aa_cm}{$cou}{$pctg}};
				#my @aa_sort;
				#foreach(@pctg){push @aa_sort,$aa_pl{$_};}
				#@aa_sort=sort {$a<=>$b} @aa_sort;
				my %bb_sort;
				foreach(@pctg){$bb_sort{$_}=$bb->{$_}[0];}#print "pp\t$_\t$bb->{$_}[0]\n";}
				my @real_sort=sort {$bb_sort{$a}<=>$bb_sort{$b}} keys %bb_sort;
				#my @real_sort;
				#if($bb_s>0){@real_sort=sort {$bb_sort{$a}<=>$bb_sort{$b}} keys %bb_sort;}
				#else{@real_sort=sort {$bb_sort{$b}<=>$bb_sort{$a}} keys %bb_sort;}
				for(my $i=0;$i<@real_sort;$i++){
					my $new_sort=$aa_pl{$aa_pl_ctg[$count]};
					$count++;
					my $new_ctg=$real_sort[$i];
					#print "$new_ctg\n";
					my @ctg_inf=@{$aa_re{$new_ctg}};
					$ctg_inf[0]=$new_sort;
					$aa->{$new_ctg}=[@ctg_inf];
				}
			}
		}
	}

}

#buil paragraph and output the longest path
sub build{
	my ($aa,$bb,$aa_rv,$bb_rv,$used,$end1,$end2,$aa_st,$bb_st,$path,$shortest)=@_;
	#foreach my $key(sort {$aa->{$a}[0]<=>$aa->{$b}[0]} keys %{$aa}){
	#$used->{$aa_st}=1;
	$path="$path\t$aa_st";$used->{$aa_st}=1;my $cc=1;
	if($bb->{$aa_st}){
		#another map path
		my $index=$bb->{$aa_st}[0];my $next_index=$bb_rv->{$index}[1];
		if($next_index>0){
			my $ctg=$bb_rv->{$next_index}[0];
			if($ctg eq $end2){
				if($aa->{$ctg}){
					my $new_aa_st=$ctg;
					my %usd=%{$used};
					my $new_path=$path;
					build($aa,$bb,$aa_rv,$bb_rv,\%usd,$end1,$end2,$new_aa_st,$new_aa_st,$new_path,$shortest);
				}else{
					my @tmp_out=split /\s+/,$path;
					if(@tmp_out>=$shortest){print "normal end:\t$path\t$ctg\n";}
				}
			}elsif($used->{$ctg}){
				#print "circle end:\t$path\t$ctg\n";
			}else{
				my $new_aa_st=$ctg;
				my %usd=%{$used};
				my $new_path=$path;
				build($bb,$aa,$bb_rv,$aa_rv,\%usd,$end2,$end1,$new_aa_st,$new_aa_st,$new_path,$shortest);
			}
		}else{
			my @tmp_out=split /\s+/,$path;
			if(@tmp_out>=$shortest+1){print "normal a end:\t$path\n";}
			$cc=0;
		}
	}
	#
		#self map path
		my $index=$aa->{$aa_st}[0];my $next_index=$aa_rv->{$index}[1];
		if($next_index>0){
			my $ctg=$aa_rv->{$next_index}[0];
			if($ctg eq $end1){
				if($bb->{$ctg}){
					my $new_aa_st=$ctg;
					my %usd=%{$used};
					my $new_path=$path;
					build($bb,$aa,$bb_rv,$aa_rv,\%usd,$end2,$end1,$new_aa_st,$new_aa_st,$new_path,$shortest);
				}else{
					my @tmp_out=split /\s+/,$path;
					if(@tmp_out>=$shortest){print "normal end:\t$path\t$ctg\n";}
				}
			}elsif($used->{$ctg}){
				#print "circle end:\t$path\t$ctg\n";
			}else{
				my $new_aa_st=$ctg;
				my %usd=%{$used};
				my $new_path=$path;
				build($aa,$bb,$aa_rv,$bb_rv,\%usd,$end1,$end2,$new_aa_st,$new_aa_st,$new_path,$shortest);
			}
		}else{
			my @tmp_out=split /\s+/,$path;
			if(@tmp_out>=$shortest+1){
				if($cc){print "normal b end:\t$path\n";}
			}
		}
	#
}

sub del_link{
	my ($aa,$bb)=@_;
	my ($aa_st,$aa_ed,$bb_st,$bb_ed);
	($aa_st,$aa_ed)=st_ed($aa);
	($bb_st,$bb_ed)=st_ed($bb);
	print STDERR "$aa_st,$aa_ed,$bb_st,$bb_ed";
	my $aa_len=abs($aa_ed-$aa_st);
	my $bb_len=abs($bb_ed-$bb_st);
	my @check;
	foreach(keys %{$aa}){
		if($bb->{$_}){push @check,$_;}
	}
	foreach my $ctg(@check){
		my $aa_cm=$aa->{$ctg}[1];
		my $bb_cm=$bb->{$ctg}[1];
		my $tmp_aa_cm=abs($aa_cm-$aa_st);
		my $tmp_bb_cm=abs($bb_cm-$bb_st);
		my $dv_aa=$tmp_aa_cm/$aa_len*100;
		my $dv_bb=$tmp_bb_cm/$bb_len*100;
		my $dv_cha=abs($dv_bb-$dv_aa);
		#my $small=(sort {$a<=>$b} ($dv_aa,$dv_bb))[0];
		#my $dv_dv_cha;
		#if($small==0){$dv_dv_cha=$dv_cha;}else{$dv_dv_cha=$dv_cha/$small*100;}
		#if($dv_dv_cha>50){
		if($dv_cha>10){
			#wrong link, keep one.check frequency
			my @f_aa;my @f_bb;my $f_aa_large;my $f_bb_large;
			for(my $i=2;$i<@{$aa->{$ctg}};$i+=2){push @f_aa,$aa->{$ctg}[$i];$f_aa_large+=$aa->{$ctg}[$i];}
			for(my $i=2;$i<@{$bb->{$ctg}};$i+=2){push @f_bb,$bb->{$ctg}[$i];$f_bb_large+=$bb->{$ctg}[$i];}
			#my $f_aa_large=(sort {$b<=>$a} @f_aa)[0];
			#my $f_bb_large=(sort {$b<=>$a} @f_bb)[0];
			if($f_aa_large>$f_bb_large){
				#keep aa
				undef $bb->{$ctg};delete $bb->{$ctg};
				print STDERR "del\taa\t$ctg\n";
			}else{
				#keep bb
				undef $aa->{$ctg};delete $aa->{$ctg};
				print STDERR "del\tbb\t$ctg\n";
			}
		}
	}
}

sub st_ed{
	my ($hash)=@_;
	my @key=sort {$hash->{$a}[0]<=>$hash->{$b}[0]} keys %{$hash};
	return $hash->{$key[0]}[1],$hash->{$key[-1]}[1];
}

sub cM_link{
	my ($aa,$bb)=@_;
	my %ovl;
	foreach(keys %{$aa}){
		if($bb->{$_}){$ovl{$_}=1;}
	}
	my %tmp_aa;my %sort_aa;my $num=1;
	foreach my $ctg(keys %{$aa}){
		if($bb->{$ctg}){
			$tmp_aa{$aa->{$ctg}[1]}{$bb->{$ctg}[1]}++;
		}
		if($sort_aa{$aa->{$ctg}[1]}){next;}else{$sort_aa{$aa->{$ctg}[1]}=$num;$num++;}
	}
	my %tmp_bb;
	foreach my $ctg(keys %{$bb}){
		if($aa->{$ctg}){
			$tmp_bb{$bb->{$ctg}[1]}{$aa->{$ctg}[1]}++;
		}
	}
	#my @rel_aa;my @rel_bb;
	my %rel_aa;my %rel_bb;
	foreach my $aa_cm(keys %tmp_aa){
		my $bigest=(sort {$tmp_aa{$aa_cm}{$b}<=>$tmp_aa{$aa_cm}{$a}} keys $tmp_aa{$aa_cm})[0];
		#
		foreach my $bb_cm(keys $tmp_aa{$aa_cm}){
			if($tmp_aa{$aa_cm}{$bb_cm}==$bigest){$rel_aa{$aa_cm}{$bb_cm}=1;}
		}
	}
	foreach my $bb_cm(keys %tmp_bb){
		my $bigest=(sort {$tmp_bb{$bb_cm}{$b}<=>$tmp_bb{$bb_cm}{$a}} keys $tmp_bb{$bb_cm})[0];
		#my @rel;
		foreach my $aa_cm(keys $tmp_bb{$bb_cm}){
			if($tmp_bb{$bb_cm}{$aa_cm}==$bigest){$rel_bb{$bb_cm}{$aa_cm}=1;}
		}
	}
	#
}

sub split_block{
	my ($aa,$bb)=@_;
	my @check;my %sort_bb;
	my %check_sort;my $num=1;
	foreach(sort {$aa->{$a}[0]<=>$aa->{$b}[0]} keys %{$aa}){
		if($bb->{$_}){
			push @check,$_;$sort_bb{$bb->{$_}[0]}=$_;
			$check_sort{$_}=$num;$num++;
		}
	}
	my $j=0;my $end=0;
	my @split;
	for(my $i=0;$i<@check;$i++){
		my $cut_bb=$bb->{$check[$i]}[0];
		my $cut_aa=$aa->{$check[$i]}[0];
		if($i>0){
			my $flag=1;
			my $juc=$bb->{$check[$i]}[0]; 
			if($juc>$end){
				$end=$juc;
				for(;$j<$juc;$j++){
					if($sort_bb{$j} && $aa->{$sort_bb{$j}}[0]>$cut_aa){
						$flag=0;
						last;
					}
				}
			}else{
				$flag=0;
			}
			if($flag){
				#can cut block
				$j=$i+1;
				print "$check[$i] is split point\n";
				push @split,$check[$i];
			}else{
				$j=0;
			}
		}
	}
	my @real_split;
	for(my $i=0;$i<@split-1;$i++){
		my $fir=$split[$i];my $sec=$split[$i+1];
		my $fir_s=$check_sort{$fir};my $sec_s=$check_sort{$sec};
		if(abs($sec_s-$fir_s)==1){
			#neibor
			next;
		}else{
			push @real_split,$fir;
			print "$fir\treal split\n";
		}
	}
	return @real_split;
}

sub merge{
	my ($aa,$bb,$block)=@_;
	my @check;
	#my %sort_bb;my %check_sort;
	foreach(sort {$aa->{$a}[0]<=>$aa->{$b}[0]} keys %{$aa}){
		if($bb->{$_}){
			push @check,$_;
			#$sort_bb{$bb->{$_}[0]}=$_;
			#$check_sort{$_}=$num;$num++;
		}
	}
	my $num=1;my $flag_aft=0;
	for(my $i=1;$i<@check;$i++){
		my $pre=$check[$i-1];my $aft=$check[$i];
		if($aa->{$aft}[0]-$aa->{$pre}[0]==1 && $bb->{$aft}[0]-$bb->{$pre}[0]==1){
			push @{$block->{$num}},$pre;$flag_aft=1;
		}else{
			push @{$block->{$num}},$pre;$num++;$flag_aft=0;
		}
	}
	if($flag_aft){push @{$block->{$num}},$check[-1];}
	my %marker;
	foreach my $sort(sort {$a<=>$b} keys %{$block}){
		if(@{$block->{$sort}}==1){
			undef $block->{$sort};delete $block->{$sort};
		}else{
			foreach(@{$block->{$sort}}){$marker{$_}="block$sort";}
		}
	}
	#merge block and resort the hash
	my $new=1;my %get;
	foreach my $ctg(sort {$aa->{$a}[0]<=>$aa->{$b}[0]} keys %{$aa}){
		if($marker{$ctg}){
			if($get{$marker{$ctg}}){
				undef $aa->{$ctg};delete $aa->{$ctg};
			}
			else{
				$aa->{$marker{$ctg}}[0]=$new;$new++;$get{$marker{$ctg}}=1;
				my @tmp=@{$aa->{$ctg}};shift @tmp;push $aa->{$marker{$ctg}},@tmp;
				undef $aa->{$ctg};delete $aa->{$ctg};
			}
		}else{
			$aa->{$ctg}[0]=$new;$new++;
		}
	}
	$new=1;%get=();
	foreach my $ctg(sort {$bb->{$a}[0]<=>$bb->{$b}[0]} keys %{$bb}){
		if($marker{$ctg}){
			if($get{$marker{$ctg}}){
				undef $bb->{$ctg};delete $bb->{$ctg};
			}
			else{
				$bb->{$marker{$ctg}}[0]=$new;$new++;$get{$marker{$ctg}}=1;
				my @tmp=@{$bb->{$ctg}};shift @tmp;push $bb->{$marker{$ctg}},@tmp;
				undef $bb->{$ctg};delete $bb->{$ctg};
			}
		}else{
			$bb->{$ctg}[0]=$new;$new++;
		}
	}
}
