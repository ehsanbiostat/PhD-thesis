#!/usr/local/bin/perl -w

use strict;

for(my $i=101; $i<6311; $i++){
  my $rfile="SA_mod.SD".$i.".r";	
  my $jobfile="job_mod.SD".$i.".sh";
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_mod.SD$i.txt\n";
  close OUT;
  `qsub -l h_vmem=4G $jobfile`;
}

