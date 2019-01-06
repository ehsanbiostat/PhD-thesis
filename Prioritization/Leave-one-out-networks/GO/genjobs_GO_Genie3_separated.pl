#!/usr/local/bin/perl -w

use strict;

for(my $i=138;$i<148;$i++){
  my $rfile="GO".$i."_Genie3_separated.r";	
  my $jobfile="job_Genie3_separated".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genie3_separated$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

