#!/usr/local/bin/perl -w

use strict;

for(my $i=5001; $i<6588; $i++){
  my $rfile="SA".$i.".r";	
  my $jobfile="job".$i.".sh";
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -l hostname=blizzard2 -l h_vmem=4.5G $jobfile`;
}

