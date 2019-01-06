#!/usr/local/bin/perl -w

use strict;

for(my $i=137;$i<148;$i++){
  my $rfile="GO".$i."_text_Leave.r";	
  my $jobfile="job".$i."_text_Leave.sh";
  my $j=$i-108;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Leave_text$i.txt\n";
  close OUT;
  #`qsub -q lowmem -l hostname=hurricane2 $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  `qsub -q test $jobfile`;
}

