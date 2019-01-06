#!/usr/local/bin/perl -w

use strict;

for(my $i=124;$i<148;$i++){
  my $rfile="Genemania_Leave".$i.".r";	
  my $jobfile="job_Genemania".$i.".sh";
  my $j=$i-99;
  open(OUT,">$jobfile");
  print OUT "module load R/\$ARCH/2.15.1\n";
  print OUT "R --vanilla < $rfile > out_Genemania$i.txt\n";
  close OUT;
  #`qsub -q lowmem -l hostname=blade$j $jobfile`;
  `qsub -q lowmem -l hostname=hurricane2 $jobfile`;
  #`qsub -q test $jobfile`;
}

