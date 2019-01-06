#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $rfile="Genie3_Leave".$i.".r";	
  my $jobfile="job_Genie3".$i.".sh";
  my $j=$i-122;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genie3$i.txt\n";
  close OUT;
  `qsub -q test $jobfile`;
}

