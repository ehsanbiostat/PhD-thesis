#!/usr/local/bin/perl -w

use strict;

for(my $i=1; $i<701; $i++){
  my $rfile="SA".$i."_Tot.SD.r";	
  my $jobfile="job".$i."_T.sh";
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Total$i.txt\n";
  close OUT;
  `qsub -l h_vmem=4G $jobfile`;
}

