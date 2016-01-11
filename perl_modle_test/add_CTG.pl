use strict;use warnings;

#huangguodong@genomics.cn 20141011

die "perl $0 <agp> <unused.fa.len> <rel> <data_dir> <out>" if(@ARGV==0);


my ($agp,$unused_scf,$rel,$data_dir,$out)=@ARGV;

my %rel;my %bacSort;my %BACinCTG;my %last_sort;
read_rel($rel,\%rel,$data_dir,\%bacSort,\%BACinCTG);

my %unused_scf_len;my %bacTolen;my %exCTG;my %add_longest;
read_scf_len($unused_scf,\%unused_scf_len,\%bacTolen,\%BACinCTG,\%exCTG,\%add_longest);

my %agpCTG;
open IN,$agp or die;
while(<IN>){
	my @aa=split;
	if($agpCTG{$aa[0]}){next;}else{$agpCTG{$aa[0]}=1;}
}
close IN;
my %thisBAC;
open OUT,">$out.addCTG.agp" or die;
foreach my $ctg(keys %rel){
	if($agpCTG{$ctg}){
		next;
	}
	if($exCTG{$ctg}){
		#missing CTG which need to supply
		my $st=1;my $ed=1;my $num=1;my $pre_sort;
		foreach my $sort(sort {$a<=>$b} keys $rel{$ctg}){
			unless($pre_sort){$pre_sort=$sort;}
			my $bac=$rel{$ctg}{$sort}[0];
			if($add_longest{$bac}){
				my $insert_scf=$add_longest{$bac}[1];
				my $insert_scf_len=$add_longest{$bac}[0];
				$thisBAC{$bac}=1;
				$ed=$st+$insert_scf_len-1;
				print OUT "$ctg\t$st\t$ed\t$num\tW\t$insert_scf\t1\t$insert_scf_len\t+\taddCTG\n";
				$st=$ed+1;$num++;
				my $gap=100;
				if($sort-$pre_sort>1){$gap=30000;}
				$ed=$st+$gap-1;
				print OUT "$ctg\t$st\t$ed\t$num\tN\t$gap\tscaffold\tno\tna\n";
				$st=$ed+1;$num++;
				$pre_sort=$sort;
			}
		}
	}
}
close OUT;
open IN,$unused_scf or die;
while(<IN>){
	my $abac=(split)[0];
	my $id=$abac;
	$abac=~s/_\d+$//;
	if($thisBAC{$abac}){
		print "$id\n";
	}
}
close IN;


my %block;my $sort=1;
open IN,"$out.addCTG.agp" or die;

while(<IN>){
	my @aa=split;my $re=$_;chomp $re;
	if($aa[4] eq "W"){
		push @{$block{$aa[0]}{$sort}},[@aa];
	}else{
		$sort++;
		if($aa[5]==30000){$sort++;}
	}
}
close IN;
open OUT,">$out.addCTG.agp.recal";
foreach my $ctg(keys %block){
	my $st=1;my $ed=1;my $num=1;
	my $pre_srt;
	foreach my $srt(sort {$a<=>$b} keys $block{$ctg}){
		unless($pre_srt){$pre_srt=$srt;}
		if($num>1){
			my $gap=100;
			if($srt-$pre_srt>1){$gap=30000;}
			$ed=$st+$gap-1;
			print OUT "$ctg\t$st\t$ed\t$num\tN\t$gap\tscaffold\tno\tna\n";
			$st=$ed+1;$num++;
		}
		foreach my $ele(@{$block{$ctg}{$srt}}){
			my @bb=@{$ele};
			my $len=$bb[2]-$bb[1]+1;
			$ed=$st+$len-1;
			$bb[1]=$st;$bb[2]=$ed;$bb[3]=$num;
			my $out=join "\t",@bb;
			print OUT $out,"\n";
			$st=$ed+1;$num++;
		}
		$pre_srt=$srt;
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





