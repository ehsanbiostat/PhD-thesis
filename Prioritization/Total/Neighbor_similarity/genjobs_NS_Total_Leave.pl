#!/usr/local/bin/perl -w

use strict;

for(my $i=141;$i<148;$i++){
  my $rfile="NS".$i."_Total_Leave.r";	
  my $jobfile="job_NS_L".$i.".sh";
  my $j=$i-2;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_NS_L$i.txt\n";
  close OUT;
  `qsub -q highmem $jobfile`;
}

