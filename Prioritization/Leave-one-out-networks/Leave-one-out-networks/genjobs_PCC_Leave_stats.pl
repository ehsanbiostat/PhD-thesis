#!/usr/local/bin/perl -w

use strict;

for(my $i=140;$i<148;$i++){
  my $rfile="PCC_Leave".$i.".r";	
  my $jobfile="job_PCC".$i.".sh";
  my $j=$i-85;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_PCC$i.txt\n";
  close OUT;
  #`qsub -q test $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  `qsub -q lowmem -l hostname=hurricane2 $jobfile`;
}

