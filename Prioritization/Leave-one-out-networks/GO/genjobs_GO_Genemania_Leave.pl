#!/usr/local/bin/perl -w

use strict;

for(my $i=128;$i<148;$i++){
  my $rfile="GO".$i."_Genemania_Leave.r";	
  my $jobfile="job".$i."_Genemania_Leave.sh";
  my $j=$i-80;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Leave_Genemania$i.txt\n";
  close OUT;
  `qsub -q lowmem -l hostname=hurricane1 $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  #`qsub -q test $jobfile`;
}

