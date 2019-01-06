#!/usr/local/bin/perl -w

use strict;

for(my $i=40;$i<41;$i++){
  my $rfile="GO".$i."_text.r";	
  my $jobfile="job_text".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_text$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

