#!/usr/local/bin/perl -w

use strict;

for(my $i=138;$i<147;$i++){
  my $rfile="Total_Leave".$i.".r";	
  my $jobfile="job".$i.".sh";
  my $j=$i-89;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out$i.txt\n";
  close OUT;
  `qsub -l hostname=hurricane1 $jobfile`;
  #`qsub -l hostname=blade$j $jobfile`;
}

