#!/usr/local/bin/perl -w

use strict;

for(my $i=141;$i<148;$i++){
  my $rfile="GO".$i."_PPI_Leave.r";	
  my $jobfile="job".$i."_PPI_Leave.sh";
  my $j=$i-126;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Leave_PPI$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  #`qsub -q test $jobfile`;
}

