#!/usr/local/bin/perl -w

use strict;

for(my $i=108;$i<109;$i++){
  my $rfile="GO".$i."_Mutant.r";	
  my $jobfile="job_Muatnt".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Mutant3$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

