#!/usr/local/bin/perl -w

use strict;

for(my $i=139;$i<148;$i++){
  my $rfile="Mutant_Leave".$i.".r";	
  my $jobfile="job_Mutant".$i.".sh";
  my $j=$i-100;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Mutant$i.txt\n";
  close OUT;
  #`qsub -q test $jobfile`;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  `qsub -q lowmem -l hostname=hurricane2 $jobfile`;
}

