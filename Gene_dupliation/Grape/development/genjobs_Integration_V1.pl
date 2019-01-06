
#!/usr/local/bin/perl -w

use strict;

for(my $i= 1001; $i< 5001; $i++){
  my $rfile="Integ".$i.".r";	
  my $jobfile="job".$i."_Integ.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Integ$i.txt\n";
  close OUT;
  `qsub -S /bin/bash -l h_vmem=1.2G $jobfile`;
}


