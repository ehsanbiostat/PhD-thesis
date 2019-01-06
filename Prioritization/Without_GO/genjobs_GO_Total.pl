#!/usr/local/bin/perl -w

use strict;

for(my $i=101;$i<148;$i++){
  my $rfile="W_GO".$i."_Total_Leave.r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-5;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
}

