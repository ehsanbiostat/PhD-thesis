#!/usr/local/bin/perl -w

use strict;

for(my $i=4; $i<10; $i++){
  my $rfile="SA".$i."_Compen.r";	
  my $jobfile="job".$i."_T.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_$i.txt\n";
  close OUT;
  `qsub -pe serial 2 -l mem_free=4G,h_vmem=2G $jobfile`;
}

