#!/usr/local/bin/perl -w

use strict;

for(my $i=131;$i<148;$i++){
  my $rfile="NS".$i."_Total.r";	
  my $jobfile="job_NS".$i.".sh";
  my $j=$i-2;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_NS$i.txt\n";
  close OUT;
  `qsub -l hostname=hurricane2 $jobfile`;
}

