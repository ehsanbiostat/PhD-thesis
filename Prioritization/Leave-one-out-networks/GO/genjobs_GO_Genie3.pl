#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<2;$i++){
  my $rfile="GO".$i."_Genie3_Leave.r";	
  my $jobfile="job_Genie3".$i."_Leave.sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genie3$i_Leave.txt\n";
  close OUT;
  `qsub -q test $jobfile`;
}

