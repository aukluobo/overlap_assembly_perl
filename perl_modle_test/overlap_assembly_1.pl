use strict;
use warnings;
use lib "/ifs1/ST_PLANT/USER/huangguodong/perlfile/perl_modle_test/";
use Filter_Overlap_BAC_1;
use Filter_Overlap_BAC_2;
use Detect_Overlap;
use Paragraph_1;
#use Resolve_AGP_1;
use Resolve_AGP_1_change;
use Agp_rep_add;
die "perl $0 <relation.lst> <bac_data_dir> <make_blast:0/1> <scaffold_len.lst> <delete.scf.lst>" if(@ARGV==0);

my ($rel,$data_dir,$make_blast,$scf_len,$del)=@ARGV;

my %rel;my %bacSort;
read_rel($rel,\%rel,$data_dir,\%bacSort);
if($make_blast){make_blast($data_dir);die "Please run blast by hand\n";}
my %scf_len;my %bacTolen;
read_scf_len($scf_len,\%scf_len,\%bacTolen);

my %del;
read_del($del,\%del);

my %used_scf;
foreach my $ctg(keys %rel){
  print "checking CTG $ctg\n";
  my $m8="./data/$ctg/blast_out.m8";
  #print "$m8\n";
  my %m8;
  Filter_Overlap_BAC_1->filter($m8,\%m8,\%del);
  open OUT,">$m8.r";
  foreach my $t_key(keys %m8){
	foreach my $t_key1(keys $m8{$t_key}){
	  foreach my $t_key2(keys $m8{$t_key}{$t_key1}){
	    my $type=$m8{$t_key}{$t_key1}{$t_key2}[4];$type=~tr/+-/-+/;
		print OUT "$t_key1\t$t_key\t90\t0\t0\t0\t$m8{$t_key}{$t_key1}{$t_key2}[2]\t$m8{$t_key}{$t_key1}{$t_key2}[3]\t$m8{$t_key}{$t_key1}{$t_key2}[0]\t$m8{$t_key}{$t_key1}{$t_key2}[1]\t0\t500\t$type\n";
		#print "$t_key1\t$t_key\t90\t0\t0\t0\t$m8{$t_key}{$t_key1}{$t_key2}[2]\t$m8{$t_key}{$t_key1}{$t_key2}[3]\t$m8{$t_key}{$t_key1}{$t_key2}[0]\t$m8{$t_key}{$t_key1}{$t_key2}[1]\t0\t500\t$type\n";
	  }
	}
  }
  close OUT;
  my $m8r="$m8.r";my %m8r;
  Filter_Overlap_BAC_2->filter($m8r,\%m8r,\%del,\%scf_len,\%bacSort);
  foreach my $t_key(keys %m8r){
    foreach my $t_key1(keys $m8r{$t_key}){
	  foreach my $t_key2(sort {$m8r{$t_key}{$t_key1}{$a}[0]<=>$m8r{$t_key}{$t_key1}{$b}[0]} keys $m8r{$t_key}{$t_key1}){
	    print STDERR "$t_key\t$t_key1\t90\t0\t0\t0\t$m8r{$t_key}{$t_key1}{$t_key2}[0]\t$m8r{$t_key}{$t_key1}{$t_key2}[1]\t$m8r{$t_key}{$t_key1}{$t_key2}[2]\t$m8r{$t_key}{$t_key1}{$t_key2}[3]\t0\t500\t$m8r{$t_key}{$t_key1}{$t_key2}[4]\n";
	  }
	}
  }
  my %dec;
  Detect_Overlap->detect(\%m8r,\%scf_len,\%bacSort,\%dec);
  my %bac;
#  foreach my $bacbac(keys %dec){
#      foreach my $stat(keys $dec{$bacbac}){
#	  print "$bacbac\t$stat\t$dec{$bacbac}{$stat}[1]\n";
#      }
#  }
  my @sort=sort {$a<=>$b} keys $rel{$ctg};
  my $preBAC=shift @sort;
  my %used_bac;my %out;
  foreach my $sort(@sort){
    #unless(defined $preBAC){$preBAC=$sort;print "$preBAC\n";;next;}
    my @preBAC=@{$rel{$ctg}{$preBAC}};
    my @nowBAC=@{$rel{$ctg}{$sort}};
    my $bacID1=$preBAC[0];
    my $bacID2=$nowBAC[0];
    my $ex_ovl_len=$rel{$ctg}{$preBAC}[3];
    print "checking\t$bacID1\t$bacID2\t$preBAC\t$sort\n";
    my @result;my $flag=1;
    my %OVE;my $Onum=1;
    if($dec{$bacID1}{$bacID2}){
      my $contro=0;
      if($dec{$bacID1}{$bacID2}{"contain"}){
        print "$bacID1\t$bacID2\t$preBAC\t$sort\tcontain\n";
        foreach my $t_out(@{$dec{$bacID1}{$bacID2}{"contain"}}){
          if($contro==8){print "\n";$contro=0;$Onum++;}
          print "$t_out\t";$contro++;
          push @{$OVE{"contain"}{$Onum}},$t_out;
        }
        print "\n";
        $flag=0;
      }
      $contro=0;$Onum=1;
      if($dec{$bacID1}{$bacID2}{"overlap"}){
        print "$bacID1\t$bacID2\t$preBAC\t$sort\toverlap\n";
        foreach my $t_out(@{$dec{$bacID1}{$bacID2}{"overlap"}}){
          if($contro==8){print "\n";$contro=0;$Onum++;}print "$t_out\t";$contro++;
            push @{$OVE{"overlap"}{$Onum}},$t_out;
        }
        print "\n";
        $flag=0;
      }
      unless($flag){print "\n";}
    }
#   elsif($dec{$bacID2}{$bacID1}){
#	  my $contro=0;print "$bacID1\t$bacID2\t$preBAC\t$sort\tcontain\n";
#	  foreach my $t_out(@{$dec{$bacID2}{$bacID1}{"contain"}}){if($contro==2){print "\n";$contro=0;}print "$t_out\t";$contro++;}
#	  $contro=0;print "$bacID1\t$bacID2\t$preBAC\t$sort\toverlap\n";
#	  foreach my $t_out(@{$dec{$bacID2}{$bacID1}{"overlap"}}){if($contro==2){print "\n";$contro=0;}print "$t_out\t";$contro++;}
#      }
    if($flag){print "need hand check discard\n";$preBAC=$sort;next;}
    my $used_scf=Paragraph_1->para(\%OVE,$ex_ovl_len,\%out,$sort);
	  my %used_scf_para=%{$used_scf};
	  foreach my $key_used_scf(keys %used_scf_para){print "used\t$key_used_scf\n";}
	  #print %out to check 
    #if no then next
    if($used_scf_para{"no"}){
      print "no out";
    }else{
  	  foreach my $key(sort {$a<=>$b} keys $out{$sort}){
  	  	#print "$out{$sort}{$key}[0][0]\t$key\t$sort\n";
  		  if($out{$sort}{$key}[0]){
  		    print "out\t0\t$key\t$sort\t";
  		    foreach my $out_key(@{$out{$sort}{$key}[0]}){
            print "$out_key\t";
  		    }
  		    print "\n";
  		  }
  		  if($out{$sort}{$key}[1]){
          #print "$out{$sort}{$key}[1][0]\t$key\t$sort\n";
          print "out\t1\t$key\t$sort\t";
          foreach my $out_key(@{$out{$sort}{$key}[1]}){
  				  print "$out_key\t";
          }
          print "\n";
        }
      }
    }
    $preBAC=$sort;
  }
  #resolve order direction and discard duplication and output agp(add the scaffold contain)
  my $agp="./data/$ctg/agp";
  #Resolve_AGP_1->re_agp(\%out,\%used_scf,\%scf_len,\%bacTolen,\%{$rel{$ctg}},$ctg,$agp);
  Resolve_AGP_1_change->re_agp(\%out,\%scf_len,\%bacTolen,\%{$rel{$ctg}},$ctg,$agp);
  #check agp: the ovl repeat scff link and add the unused scaf;
  Agp_rep_add->good($agp,\%used_scf,\%scf_len,\%{$rel{$ctg}},$ctg,\%bacTolen,\%bacSort);
}
#output unused scaffold for check
open UN,">unused_scaffold_list.lst" or die;
foreach my $key_used(keys %used_scf){print "$key_used\n";}
foreach my $scf_key(keys %scf_len){
 if($used_scf{$scf_key}){
    #print UN "$scf_key\tused\n";
    next;
  }
  print UN "$scf_key\t$scf_len{$scf_key}\n"
}
close UN;





sub read_scf_len{
  my ($len,$hash,$tolen)=@_;
  open IN,$len or die;
  while(<IN>){
    my @aa=split;
    $hash->{$aa[0]}=$aa[1];
    $aa[0]=~s/_\d+$//;
    $tolen->{$aa[0]}+=$aa[1];
  }
  close IN;
}

sub make_blast{
  my ($data_dir)=@_;
  my $pwd=`pwd`;chomp $pwd;
  $data_dir=~s/\/$//;
  my @ctg_dir=glob("$data_dir/*");
  unless(@ctg_dir){die;}
  open CAT,">cat_data.sh" or die;
  open OUT,">blastall.sh" or die;
  mkdir("./data")unless(-d "./data");
  foreach my $ctg(@ctg_dir){
    chomp $ctg;#print "$ctg\n";die;
	$ctg=~s/$data_dir//;$ctg=~s/\///g;
	mkdir("./data/$ctg")unless(-d "./data/$ctg");
	print CAT "cat $data_dir/$ctg/*.fasta >$pwd/data/$ctg/$ctg.all.fa;\n";
	print OUT "cd $pwd/data/$ctg/;/opt/blc/genome/biosoft/blast-2.2.26/bin/formatdb -i $ctg.all.fa -p F;";
	print OUT "/opt/blc/genome/biosoft/blast-2.2.26/bin/megablast -d $pwd/data/$ctg/$ctg.all.fa -i $pwd/data/$ctg/$ctg.all.fa -e 1e-10 -m 8 -o $pwd/data/$ctg/blast_out.m8 -v 100 -b 100;\n";
  }
  close OUT;
  close CAT;
}


sub read_rel{
  my ($rel,$hash,$dir,$bacSort)=@_;
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
      }
    }
  }
  close IN;
}

sub read_del{
  my ($file,$hash)=@_;
  open IN,$file or die;
  while(<IN>){
    if(/#/ or /keep/){next;}
    chomp;$hash->{$_}=1;
  }
  close IN;
}




