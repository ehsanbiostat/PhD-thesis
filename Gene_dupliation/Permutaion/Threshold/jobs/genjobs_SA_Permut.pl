#!/usr/local/bin/perl -w

use strict;

for(my $i=1; $i<16; $i++){
  my $rfile="SA".$i."_Permut.r";	
  my $jobfile="job".$i."_Permut.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Perm$i.txt\n";
  close OUT;
  # `qsub -l h_vmem=2G $jobfile`;
}

