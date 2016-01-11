use strict;use warnings;
package OVL_four;

#huangguodong@genomics.cn

sub check{
	my ($self,$st1,$ed1,$st2,$ed2)=@_;
	my %kk;
	$kk{"1"}=$st1;$kk{"2"}=$ed1;
	$kk{"3"}=$st2;$kk{"4"}=$ed2;
	my @key=sort {$kk{$a}<=>$kk{$b}} keys %kk;
	#foreach(@key){print "$kk{$_}\t";}print "\n";
	my $key=join "",@key;
	if($key=~/^12/ or $key=~/^34/){
		return 0;
	}else{
		return 1;
	}
}
1


