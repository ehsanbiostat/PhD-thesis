#!/usr/local/bin/perl -w

use strict;

for(my $i=18001; $i<19826; $i++){
  my $rfile="SA".$i."_Permut.SD.r";	
  my $jobfile="job".$i."_T_Permut.sh";
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Perm$i.txt\n";
  close OUT;
  `qsub -l h_vmem=2G $jobfile`;
}

