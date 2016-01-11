use strict;use warnings;

#huangguodong@genomics.cn 20140928

die "perl $0 <check.lst> <map1.check> <map2.check> <map1 sort:1/-1> <map2 sort:1/-1>" if(@ARGV==0);

my ($result,$map1,$map2,$map1sort,$map2sort)=@ARGV;
#firstly find out the consevus block and make them marker
open IN,$result or die;
my %check;
while(<IN>){
	chomp;
	/\d+/;
	my $id=$&;
	$check{$id}=$_;
}
close IN;
my %allmap;my $num=1;
foreach my $id(sort {$a<=>$b} keys %check){
	open IN,$check{$id} or die;
	while(<IN>){
		my @aa=split;
		foreach my $ele(@aa){
			if($ele=~/CTG/ or $ele=~/block/){
				$allmap{$ele}=$num;$num++;
			}
		}
	}
	close IN;
}



my %map1;my %map2;my %map1_ctg_sort;my %map2_ctg_sort;
read_map($map1,\%map1,\%map1_ctg_sort);
read_map($map2,\%map2,\%map2_ctg_sort);

#foreach my $ctg(sort {$allmap{$a}<=>$allmap{$b}} keys %allmap){
#	if($map1{$ctg}[0]){print "$map1{$ctg}[1]\t";}
#	if($map2{$ctg}[0]){print "$map2{$ctg}[1]\t";}
#	print "\n";
#}

#find the missing CTG and the pre CTG which on conservus map
my %miss1;my %miss2;
check_miss(\%map1,\%allmap,\%miss1);
check_miss(\%map2,\%allmap,\%miss2);

my %marker1;my %marker2;my %head;

if(keys %miss1>0){
	#1 have miss
	my @key=sort {$a<=>$b} keys %miss1;
	foreach my $cu(@key){
		my @ctg=sort {$miss1{$cu}{$a}[0]<=>$miss1{$cu}{$b}[0]} keys $miss1{$cu};
		my $marker_ctg_num=$map1{$ctg[0]}[0]-1;
		my $marker_ctg=$map1_ctg_sort{$marker_ctg_num};
		foreach my $ctg(@ctg){
			push @{$marker1{$marker_ctg}},$miss1{$cu}{$ctg}[1];
		}
	}
}
if(keys %miss2>0){
	#2 have miss
	my @key=sort {$a<=>$b} keys %miss2;
	foreach my $cu(@key){
		my @ctg=sort {$miss2{$cu}{$a}[0]<=>$miss2{$cu}{$b}[0]} keys $miss2{$cu};
		my $marker_ctg_num=$map2{$ctg[0]}[0]-1;
		if($marker_ctg_num<1){
			#only map2 can meet head missing
			my $new_marker_ctg_num=$map2{$ctg[-1]}[0]+1;
			my $marker_ctg=$map2_ctg_sort{$new_marker_ctg_num};
			foreach my $ctg(@ctg){
				push @{$head{$marker_ctg}},$miss2{$cu}{$ctg}[1];
			}
		}else{
			my $marker_ctg=$map2_ctg_sort{$marker_ctg_num};
			foreach my $ctg(@ctg){
				push @{$marker2{$marker_ctg}},$miss2{$cu}{$ctg}[1];
			}
		}
	}
}
#another special case the start ctg

#output the insert
my %have;my %ex;
foreach my $ctg(sort {$allmap{$a}<=>$allmap{$b}} keys %allmap){
	if($head{$ctg}[0]){
		#add map2 head
		foreach my $rec(@{$head{$ctg}}){
			print $rec,"\tadd2head\n";
			$ex{(split /\s+/,$rec)[0]}=1;
		}
	}
	if($map1{$ctg}[0]){print "$map1{$ctg}[1]\t";}
	if($map2{$ctg}[0]){print "$map2{$ctg}[1]\t";}
	print "\n";
	if($marker1{$ctg}[0]){
		if($marker2{$ctg}[0]){
			#both
			my %c1;
			my @aa=@{$marker1{$ctg}};
			my $blo=1;my $alo=1;my $pre_cm;
			foreach my $tmp(@aa){
				my @bb=split /\s+/,$tmp;
				unless($pre_cm){$pre_cm=$bb[1];}
				if($pre_cm != $bb[1]){$alo++;$pre_cm=$bb[1];}
				$c1{$bb[1]}{$alo}{$blo}=$tmp;
				$blo++;
			}
			my %c2;
			@aa=@{$marker2{$ctg}};
			$blo=1;$alo=1;$pre_cm=0;
			foreach my $tmp(@aa){
				my @bb=split /\s+/,$tmp;
				unless($pre_cm){$pre_cm=$bb[1];}
				if($pre_cm != $bb[1]){$alo++;$pre_cm=$bb[1];}
				$c2{$bb[1]}{$alo}{$blo}=$tmp;
				$blo++;
			}
			#output
			my @key_c1=sort {$a<=>$b} keys %c1;
			my @key_c2=sort {$a<=>$b} keys %c2;
			my @first;my @second;my %fir;my %sec;
			if(@key_c1>@key_c2){
				#output c1 first
				@first=@key_c1;@second=@key_c2;%fir=%c1;%sec=%c2;
			}else{
				@first=@key_c2;@second=@key_c1;%fir=%c2;%sec=%c1;
			}
			for(my $i=0;$i<@first;$i++){
				#outpurt c1
				my @cc=sort {$a<=>$b} keys $fir{$first[$i]};
				foreach my $alo(@cc){
					foreach my $blo(sort {$a<=>$b} keys $fir{$first[$i]}{$alo}){
						my @tmpp=split /\s+/,$fir{$first[$i]}{$alo}{$blo};
						if($ex{$tmpp[0]}){print STDERR "$tmpp[0] have\n";next;}
						else{
							print $fir{$first[$i]}{$alo}{$blo},"\tadd1\n";$ex{$tmpp[0]}=1;
						}
					}
				}
				if($i<@second){
					@cc=sort {$a<=>$b} keys $sec{$second[$i]};
					foreach my $alo(@cc){
						foreach my $blo(sort {$a<=>$b} keys $sec{$second[$i]}{$alo}){
							my @tmpp=split /\s+/,$sec{$second[$i]}{$alo}{$blo};
							if($ex{$tmpp[0]}){print STDERR "$tmpp[0] have\n";next;}
							else{
								print $sec{$second[$i]}{$alo}{$blo},"\tadd2\n";$ex{$tmpp[0]}=1;
							}
						}
					}
				}
			}
			$have{$ctg}=1;
		}else{
			#just 1
			foreach my $rec(@{$marker1{$ctg}}){
				my @tmpp=split /\s+/,$rec;
				if($ex{$tmpp[0]}){print STDERR "$tmpp[0] have\n";next;}
				else{print $rec,"\tadd1\n";$ex{$tmpp[0]}=1;}
			}
			$have{$ctg}=1;
		}
	}else{
		if($marker2{$ctg}[0]){
			#just 2
			foreach my $rec(@{$marker2{$ctg}}){
				my @tmpp=split /\s+/,$rec;
				if($ex{$tmpp[0]}){print STDERR "$tmpp[0] have\n";next;}
				else{print $rec,"\tadd2\n";$ex{$tmpp[0]}=1;}
			}
			$have{$ctg}=1;
		}
	}
}

foreach my $ctg(keys %have){
	if($marker1{$ctg} or $marker2{$ctg}){
		next;
	}else{
		print STDERR "miss\t$ctg\n";
	}
}


sub check_miss{
	my ($map,$allmap,$miss)=@_;
	my $cu=1;
	foreach my $ctg(sort {$map->{$a}[0]<=>$map->{$b}[0]} keys %{$map}){
		if($allmap->{$ctg}){
			#print "have\t$ctg\n";
			$cu++;
		}else{		
			$miss->{$cu}{$ctg}=[@{$map->{$ctg}}];
		}
	}
}


sub read_map{
	my ($map,$hash,$hash1)=@_;
	my $num=1;
	open IN,$map or die;
	while(<IN>){
		my @aa=split;chomp;
		$hash->{$aa[0]}=[$num,$_];
		$hash1->{$num}=$aa[0];
		$num++;
	}
	close IN;
}