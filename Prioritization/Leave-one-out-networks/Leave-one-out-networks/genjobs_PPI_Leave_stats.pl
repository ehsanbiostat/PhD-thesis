#!/usr/local/bin/perl -w

use strict;

for(my $i=89;$i<93;$i++){
  my $rfile="PPI_Leave".$i.".r";	
  my $jobfile="job_PPI".$i.".sh";
  my $j=$i-74;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_PPI$i.txt\n";
  close OUT;
  #`qsub -q test $jobfile`;
  `qsub -q lowmem -l hostname=blade$j $jobfile`;
  #`qsub -q lowmem -l hostname=hurricane2 $jobfile`;
}

