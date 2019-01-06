#!/usr/local/bin/perl -w

use strict;

for(my $i=137;$i<147;$i++){
  my $rfile="Total_Leave".$i.".r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-122;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -l hostname=blade$j $jobfile`;
}

