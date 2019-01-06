#!/usr/local/bin/perl -w

use strict;


use strict;
print "module load gridengine\n";
my $i = 1;
while ($i <= 19589){
`/group/biocomp/projects/group_tools/cj | gawk '/ehsab/' | awk '{print \$2}' > /ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt`;
open (MYFILE, '</ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/jobs_running.txt') or die "$!";
my $line = <MYFILE>;
if (not defined $line) {
	$line = 0
}
if ($line < 100) {
	$i = $i + 1;

  my $rfile="SA".$i."_Tot.SD.r";	
  my $jobfile="job".$i."_T.sh";
  open(OUT,">$jobfile");
  print OUT "module load R\n";
  print OUT "R --vanilla < $rfile > out_Total$i.txt\n";
  close OUT;
  `qsub -l h_vmem=1G $jobfile`;


}
# print $line;
print "Job number is $i\n";
}

