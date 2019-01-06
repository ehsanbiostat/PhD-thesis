#!/usr/local/bin/perl -w

use strict;

for(my $i=49;$i<57;$i++){
  my $rfile="GO".$i."_IYG.r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-20;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -l hostname=tanith $jobfile`;
}

