#!/usr/local/bin/perl -w

use strict;

for(my $i=121;$i<148;$i++){
  my $rfile="T".$i."_Total_Leave.r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-75;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -l hostname=hurricane2 $jobfile`;
}

