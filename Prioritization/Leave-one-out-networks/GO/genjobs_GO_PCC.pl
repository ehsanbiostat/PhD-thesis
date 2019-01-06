#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $rfile="GO".$i."_PCC.r";	
  my $jobfile="job_PCC".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_PCC$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

