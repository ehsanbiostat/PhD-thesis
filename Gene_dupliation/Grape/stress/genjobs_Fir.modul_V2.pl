#!/usr/local/bin/perl -w

use strict;


use strict;
print "module load gridengine\n";
my $i = 7752;
while ($i <= 7753){
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | awk '{print \$2}' > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | gawk '{print \$3}' | sed 's/[^0-9]*//g' >> /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
`awk '{ sum += \$1 } END { print sum }' /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt`;
open (MYFILE, '</ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running1.txt') or die "$!";
my $line = <MYFILE>;

if (not defined $line) {
	 $line = 0
}
if ($line < 200) {
	$i = $i + 1;

  my $rfile="Fir".$i.".modul.fast.grape.r";	
  my $jobfile="job".$i."_Fir.modul.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Fir.modul$i.txt\n";
  close OUT;
  `qsub -S /bin/bash -l h_vmem=0.7G $jobfile`;
  
}
# print $line;
print "Job number is $i\n";
}

