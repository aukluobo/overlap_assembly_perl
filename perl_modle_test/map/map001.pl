use strict;
use warnings;

die "perl $0 <lst> <map>" if(@ARGV==0);

# find the most possible cM for CTG (if the neiborghood is the same then keep both record)

my $file=$ARGV[0];

my %all;
open IN,$file or die;
while(<IN>){
	my @aa=split;
	$all{$aa[2]}{$aa[0]}++;
}
close IN;

my $map=$ARGV[1];
my %sort;my $sort=1;
open IN,$map or die;
while(<IN>){
	my @aa=split;
	if($sort{$aa[-1]}){
		next;
	}else{
		$sort{$aa[-1]}=$sort;
		$sort++;
	}
}
close IN;


foreach my $ctg(keys %all){
	#find the largest cm. 
	my @keep;
	my $largest=(sort {$all{$ctg}{$a}<=>$all{$ctg}{$b}} keys $all{$ctg})[-1];
	my @cm;
	foreach my $cm(sort {$a<=>$b} keys $all{$ctg}){
		if($cm == $largest){push @cm,$cm;}
	}
	if(@cm==1){
		#only one  find neibor
		my $index=0;
		my @key=sort {$a<=>$b} keys $all{$ctg};
		foreach my $cm(@key){
			if($cm == $cm[0]){
				my $pre=$index-1;
				
				my $pre_sort=$sort{$key[$pre]};
				
				my $now=$sort{$key[$index]};
				my $pre_cha=$now-$pre_sort;
				if($pre_cha==1){
					#keep pre
					push @keep,$key[$pre];
				}
				push @keep,$cm;
				my $aft=$index+1;
				if($aft<@key){
					my $aft_sort=$sort{$key[$aft]};my $aft_cha=$aft_sort-$now;
					if($aft_cha==1){
						#keep aft
						push @keep,$key[$aft];
					}	
				}
				last;
			}
			$index++;
		}
		print $ctg,"\t";
		foreach my $out(@keep){
			print "$out\t$all{$ctg}{$out}\t";
		}
		print "\n";
	}else{
		#have much cm. output all.
		print $ctg,"\t";
		foreach my $ccm(sort {$a<=>$b} @cm){
			print "$ccm\tall{$ctg}{$ccm}\t"; 
		}
		print "\n";
	}
}


