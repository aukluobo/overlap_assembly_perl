use strict;
use warnings;

#huangguodong@genomics.cn 20141029

die "perl $0 <m8.best.nosame> <all_scf_len> <all_last_agp.lst> <rel>" if(@ARGV==0);

# this program si try to identify the redundant copy of the assembly sequence.
# the condition is the two sequence can match to each other totally
# the input m8 need to filter out column 1 eq column 2

my ($m8,$scf_len,$agp,$rel)=@ARGV;
my %rel;my %bacSort;my %BACinCTG;
read_rel($rel,\%rel,\%bacSort,\%BACinCTG);

my %scf_len;
open IN,$scf_len or die;
while(<IN>){
	my @aa=split;
	$scf_len{$aa[0]}=$aa[1];
}
close IN;

my %agp;
open IN,$agp or die;
while(<IN>){
	if(/bac/){
		my @aa=split;
		$agp{$aa[5]}=1;
	}
}
close IN;

open IN,$m8 or die;
my $cut=50;my %del;
while(<IN>){
	my @aa=split;my $rec=$_;chomp $rec;
	#bac_63_154_4    bac_53_359_1    401     1182    22771   23553   -       65
	my $id1=$aa[0];$id1=~s/_\d+$//;
	my $id2=$aa[1];$id2=~s/_\d+$//;
	#if($id1 eq $id2){
	#if(abs($bacSort{$id1}-$bacSort{$id2})<=3){
		#two case: contain and equal
		my ($id1_st,$id1_ed,$id2_st,$id2_ed);
		if($aa[2]<=$cut){$id1_st=1;}
		if($scf_len{$aa[0]}-$aa[3]<=$cut){$id1_ed=1;}
		if($aa[4]<=$cut){$id2_st=1;}
		if($scf_len{$aa[1]}-$aa[5]<=$cut){$id2_ed=1;}
		if($id1_st && $id1_ed){
			#id1 contain
			if($id2_st && $id2_ed){
				#id2 contain; one is redundant.
				#keep the one in agp.if not keep the long one
				my %keep;
				foreach(($aa[0],$aa[1])){
					if($agp{$_}){next;}
					else{$keep{$_}=1;}
				}
				my @keep=keys %keep;
				if(@keep == 1){
					print "delete $keep[0]\t$scf_len{$keep[0]}\n";$del{$keep[0]}=1;
				}
				elsif(@keep==2){
					if($scf_len{$aa[0]}>$scf_len{$aa[1]}){
						print "delete $aa[1]\t$scf_len{$aa[1]}\n";$del{$aa[1]}=1;
					}else{print "delete $aa[0]\t$scf_len{$aa[0]}\n";$del{$aa[0]}=1;}
				}else{
					print "$aa[0]\t$aa[1]\tboth in agp need check\n";$del{$aa[0]}=1;$del{$aa[1]}=1;
				}
			}else{
				#only id1 contain
				print "delete $aa[0] $scf_len{$aa[0]} only\n";$del{$aa[0]}=1;
			}
		}else{
			if($id2_st && $id2_ed){
				#only id2 contain
				print "delete $aa[1] $scf_len{$aa[1]} only\n";$del{$aa[1]}=1;
			}
		}
		if(($id1_st or $id1_ed) && ($id2_st or $id2_ed)){
			#may get overlap relationship
			if($del{$aa[0]} or $del{$aa[1]}){next;}
			if($aa[3]-$aa[2]<300 or $aa[5]-$aa[4]<300){next;}
			print STDERR "ovl\t$scf_len{$aa[0]}\t$scf_len{$aa[1]}\t";
			if($agp{$aa[0]}){print STDERR "1\t";}else{print STDERR "0\t";}
			if($agp{$aa[1]}){print STDERR "1\t";}else{print STDERR "0\t";}
			print STDERR "$rec\n";
		}

	#}
}
close IN;



sub read_rel{
  	my ($rel,$hash,$bacSort,$BACinCTG)=@_;
  	open IN,$rel or die;
  	while(<IN>){
    	if(/bac/){
     	my @aa=split;
      		if($aa[-1]=~/(CTG\d+)_(bac.*)_(\d+)$/){
		        my $ctg=$1;my $bac=$2;my $sort=$3;
				my $ex=0;
				#bac_id exists_or_not expect_len overlap    $aa[1] is expect total length; $aa[3] is expect overlap length
				$hash->{$ctg}{$sort}=[$bac,$ex,$aa[1],$aa[3]];
				$bacSort->{$bac}=$sort;
				$BACinCTG->{$bac}=$ctg;
      		}
    	}
  	}
  	close IN;
}







