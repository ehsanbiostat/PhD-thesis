#!/usr/local/bin/perl -w

use strict;

for(my $i=135;$i<148;$i++){
  my $rfile="GO".$i."_Total_Leave.r";	
  my $jobfile="job_le".$i.".sh";
  my $j=$i-30;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_le$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
}

