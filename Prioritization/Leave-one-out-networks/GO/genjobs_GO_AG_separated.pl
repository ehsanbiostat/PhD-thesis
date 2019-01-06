#!/usr/local/bin/perl -w

use strict;

for(my $i=81;$i<148;$i++){
  my $rfile="GO".$i."_AG_separated.r";	
  my $jobfile="job".$i."_AG_separated.sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_AG_separated$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

