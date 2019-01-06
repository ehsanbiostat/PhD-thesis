#!/usr/local/bin/perl -w

use strict;

for(my $i=129;$i<148;$i++){
  my $rfile="GO".$i."_PCC_Leave.r";	
  my $jobfile="job_PCC".$i."_Leave.sh";
  my $j=$i-119;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_PCC$i.txt\n";
  close OUT;
  #`qsub -q lowmem -l hostname=hurricane1 $jobfile`;
  `qsub -q lowmem -l hostname=blade$j $jobfile`;
  #`qsub -q test $jobfile`;
}

