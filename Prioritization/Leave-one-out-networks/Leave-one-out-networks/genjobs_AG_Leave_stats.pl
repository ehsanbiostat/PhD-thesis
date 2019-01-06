#!/usr/local/bin/perl -w

use strict;

for(my $i=13;$i<15;$i++){
  my $rfile="AG_Leave".$i.".r";	
  my $jobfile="job_AG".$i.".sh";
  my $j=$i-102;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_AG$i.txt\n";
  close OUT;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  `qsub -q lowmem -l hostname=hurricane2 $jobfile`;
}

