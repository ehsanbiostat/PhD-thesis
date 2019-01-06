#!/usr/local/bin/perl -w

use strict;

for(my $i=2; $i<3; $i++){
  my $rfile="SA".$i."_Permut.SD_V2.r";	
  my $jobfile="job".$i."_T_Permut.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Perm$i.txt\n";
  close OUT;
  `qsub -pe serial 30 -l mem_free=150G,h_vmem=4G $jobfile`;
}

