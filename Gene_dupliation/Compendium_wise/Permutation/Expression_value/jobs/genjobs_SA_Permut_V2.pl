#!/usr/local/bin/perl -w

use strict;

for(my $i=16700; $i<19286; $i++){
  my $rfile="SA".$i."_Permut.SD_V2_Ind.r";	
  my $jobfile="job".$i."_T_Permut_Ind.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Perm$i.txt\n";
  close OUT;
  `qsub -l mem_free=1.1G $jobfile`;
}

