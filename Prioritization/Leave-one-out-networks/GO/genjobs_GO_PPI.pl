#!/usr/local/bin/perl -w

use strict;

for(my $i=29;$i<148;$i++){
  my $rfile="GO".$i."_PPI.r";	
  my $jobfile="job_PPI".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_PPI$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

