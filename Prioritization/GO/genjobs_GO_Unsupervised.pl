#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $rfile="GO".$i."_Unsupervised_Leave.r";	
  my $jobfile="job".$i."_Unsup.sh";
  my $j=$i-89;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -q lowmem $jobfile`;
}

