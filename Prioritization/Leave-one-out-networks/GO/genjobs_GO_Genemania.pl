#!/usr/local/bin/perl -w

use strict;

for(my $i=65;$i<88;$i++){
  my $rfile="GO".$i."_Genemania.r";	
  my $jobfile="job_Genemania".$i.".sh";
  my $j=$i-60;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genemania$i.txt\n";
  close OUT;
  `qsub $jobfile`;
}

