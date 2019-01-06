#!/usr/local/bin/perl -w

use strict;

for(my $i=14;$i<19;$i++){
  my $rfile="GO".$i."_Genie3_Leave.r";	
  my $jobfile="job_Genie3".$i."_Leave.sh";
  my $j=$i-108;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genie3_Leave$i.txt\n";
  close OUT;
  #`qsub -q lowmem -l hostname=hurricane1 $jobfile`;
  #`qsub -q lowmem -l hostname=blade0$j $jobfile`;
  `qsub -q test $jobfile`;
}

