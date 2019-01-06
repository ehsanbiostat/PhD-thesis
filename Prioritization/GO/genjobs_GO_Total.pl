#!/usr/local/bin/perl -w

use strict;

for(my $i=134;$i<147;$i++){
  my $rfile="GO".$i."_Total.r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-89;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
}

