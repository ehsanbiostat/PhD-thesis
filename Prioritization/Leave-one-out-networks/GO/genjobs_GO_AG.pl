#!/usr/local/bin/perl -w

use strict;

for(my $i=131;$i<133;$i++){
  my $rfile="GO".$i."_AG.r";	
  my $jobfile="job_AG".$i.".sh";
  my $j=$i-107;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_AG$i.txt\n";
  close OUT;
  `qsub -l hostname=hurricane1 $jobfile`;
}

