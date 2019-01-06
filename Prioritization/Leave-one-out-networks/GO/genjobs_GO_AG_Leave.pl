#!/usr/local/bin/perl -w

use strict;

for(my $i=64;$i<65;$i++){
  my $rfile="GO".$i."_AG_Leave.r";	
  my $jobfile="job".$i."_AG_Leave.sh";
  my $j=$i-114;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Leave_AG$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
}

