use strict;use warnings;
#huangguodong@genomics.cn 20141226

#1\add 10kbp in between each agp
#2\change the N to be 5kbp
#3\change the scaffold name use the prefix of the agp file

die "perl $0 <input_chr:CHR> <agp.file> ...." if(@ARGV==0); 

my ($chrname,@file)=@ARGV;

my $st=1;my $ed=1;my $num=1;
foreach my $agpfile(@file){
	chomp $agpfile;
	open IN,$agpfile or die;
	my $prefix=$agpfile;$prefix=~s/\..*//;
	while(<IN>){
		chomp;
		my @aa=split;
		if($aa[4] eq "W"){
			my $len=$aa[2]-$aa[1]+1;
			$ed=$st+$len-1;
			$aa[5]="$prefix\_$aa[5]";
		}elsif($aa[4] eq "N"){
			my $len=5000;
			$ed=$st+$len-1;
			$aa[5]=5000;
		}
		$aa[2]=$ed;$aa[1]=$st;$aa[0]=$chrname;$aa[3]=$num;$num++;
		$st=$ed+1;
		my $out=join "\t",@aa;
		print $out,"\n";
	}
	close IN;
	$ed=$st+10000-1;
	#3AS     1225503 1225602 2       N       100     scaffold        no      map
	print "$chrname\t$st\t$ed\t$num\tN\t10000\tscaffold\tno\tmap\n";
	$st=$ed+1;$num++;
}







