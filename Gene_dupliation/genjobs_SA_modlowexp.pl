#!/usr/local/bin/perl -w

use strict;

for(my $i=6001; $i<6239; $i++){
  my $rfile="SA_mod".$i.".r";	
  my $jobfile="job_mod".$i.".sh";
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_mod$i.txt\n";
  close OUT;
  `qsub -l hostname=crunch1 -l h_vmem=4G $jobfile`;
}

